export SoilPlantHydrologyModel
"""
    struct SoilPlantHydrologyModel{
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        VM <: PlantHydraulics.AbstractVegetationModel{FT},
    } <: AbstractLandModel{FT}
        soil::SM
        vegetation::VM
    end

A concrete type of land model used for simulating systems with a 
vegetation and a soil component.
$(DocStringExtensions.FIELDS)
"""
struct SoilPlantHydrologyModel{
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    VM <: PlantHydraulics.AbstractVegetationModel{FT},
} <: AbstractLandModel{FT}
    "The soil model to be used"
    soil::SM
    "The vegetation model to be used"
    vegetation::VM
end


"""
    SoilPlantHydrologyModel{FT}(;
                      soil_model_type::Type{SM},
                      soil_args::NamedTuple = (;),
                      vegetation_model_type::Type{VM},
                      vegetation_args::NamedTuple = (;),
                      ) where {FT,
                               SM <: Soil.AbstractSoilModel{FT},
                               VM <: PlantHydraulics.AbstractVegetationModel{FT}}

A constructor for the `SoilPlantHydrologyModel`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `SoilPlantHydrologyModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function SoilPlantHydrologyModel{FT}(;
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    vegetation_model_type::Type{VM},
    vegetation_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    VM <: PlantHydraulics.AbstractVegetationModel{FT},
}
    # These may be passed in, or set, depending on use scenario
    boundary_fluxes = soil_args.boundary_conditions #FluxBC{FT}(FT(0.0), FT(0.0))
    transpiration = vegetation_args.transpiration #PrescribedTranspiration{FT}((t::FT) -> FT(0.0))

    ##These should always be set by the constructor.
    sources = (RootExtraction{FT}(),)
    root_extraction = PrognosticSoilPressure{FT}()

    soil = soil_model_type(;
        boundary_conditions = boundary_fluxes,
        sources = sources,
        soil_args...,
    )
    vegetation = vegetation_model_type(;
        root_extraction = root_extraction,
        transpiration = transpiration,
        vegetation_args...,
    )
    args = (soil, vegetation)
    return SoilPlantHydrologyModel{FT, typeof.(args)...}(args...)
end

interaction_vars(m::SoilPlantHydrologyModel) = (:root_extraction,)

interaction_types(m::SoilPlantHydrologyModel{FT}) where {FT} = (FT,)

interaction_domains(m::SoilPlantHydrologyModel) = (:subsurface,)

"""
    make_interactions_update_aux(
        land::SoilPlantHydrologyModel{FT, SM, RM},
    ) where {FT, SM <: Soil.RichardsModel{FT}, RM <: PlantHydraulics.PlantHydraulicsModel{FT}}

A method which makes a function; the returned function 
updates the auxiliary variable `p.root_extraction`, which
is needed for both the boundary condition for the soil model and the source
term (runoff) for the surface water model.

This function is called each ode function evaluation.
"""
function make_interactions_update_aux(#Do we want defaults, for land::AbstractLandModel?
    land::SoilPlantHydrologyModel{FT, SM, RM},
) where {
    FT,
    SM <: Soil.RichardsModel{FT},
    RM <: PlantHydraulics.PlantHydraulicsModel{FT},
}
    function update_aux!(p, Y, t)
        @unpack plant_vg_α, plant_vg_n, plant_vg_m, plant_ν, plant_S_s, plant_K_sat = land.vegetation.parameters
        #@show(land.vegetation.domain.root_depths)
        #@show(typeof(land.vegetation.domain.root_depths))
        #@show(land.soil.domain)
        #@show(land.soil.domain.zlim)
        #@show(land.vegetation.domain.root_depths)        
        #@show(land.vegetation.domain.compartment_midpoints[1])
        @show(p.soil.ψ)
        @show(Y.vegetation.ϑ_l)
        @show(Y.soil.ϑ_l)
        @show(PlantHydraulics.effective_saturation.(plant_ν,Y.vegetation.ϑ_l))
        @show(PlantHydraulics.water_retention_curve.(plant_vg_α,plant_vg_n,plant_vg_m,PlantHydraulics.effective_saturation.(plant_ν,Y.vegetation.ϑ_l),plant_ν,plant_S_s))               
        @show(plant_vg_α)
        @show(plant_vg_n)
        @show(plant_vg_m)
        @show(plant_ν)
        @show(plant_S_s)
        @show(plant_K_sat[:root])
        @show(plant_K_sat[land.vegetation.domain.compartment_labels[1]])
        @show(ClimaCore.Fields.coordinate_field(land.soil.domain.space).z)

        @. p.root_extraction = PlantHydraulics.flux.(ClimaCore.Fields.coordinate_field(land.soil.domain.space).z, #land.soil.domain.zlim[1], #, #soil.coordinates.z, #
                land.vegetation.domain.compartment_midpoints[1],
                p.soil.ψ,
                PlantHydraulics.water_retention_curve(plant_vg_α,plant_vg_n,plant_vg_m,PlantHydraulics.effective_saturation(plant_ν,Y.vegetation.ϑ_l[1]),plant_ν,plant_S_s),               
                plant_vg_α,
                plant_vg_n,
                plant_vg_m,
                plant_ν,
                plant_S_s,
                plant_K_sat[:root],
                plant_K_sat[land.vegetation.domain.compartment_labels[1]],
        ).* land.vegetation.parameters.root_distribution.(ClimaCore.Fields.coordinate_field(land.soil.domain.space).z) #land.vegetation.domain.root_depths land.soil.coordinates.z,) 
        
        #=
        @. p.root_extraction =
            Roots.ground_area_flux.(
                # we have to use coordinates here rather than root depth array to avoid broadcast error? 
                land.soil.coordinates.z,
                land.vegetation.domain.compartment_heights[1],
                p.soil.ψ .* FT(9800), # Pa: ρg * meters
                Roots.θ_to_p(Y.vegetation.θ[1]),
                a_root,
                b_root,
                K_max_root,
                land.vegetation.param_set.RAI,
            ) .*
            land.vegetation.param_set.root_distribution_function.(
                land.soil.coordinates.z,
            )
        =#
    end
    return update_aux!
end


## Extending methods of the Plant Hydraulics and Soil Model
## TBD if these should be linked more explicitly.
"""
   PrognosticSoilPressure{FT} <: PlantHydraulics.AbstractRootExtraction{FT}

Concrete type of PlantHydraulics.AbstractRootExtraction, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Soil.RootExtraction`:both 
are used at the same time,
ensuring that the water flux into the roots is extracted correctly
from the soil.
"""
struct PrognosticSoilPressure{FT} <: PlantHydraulics.AbstractRootExtraction{FT} end

"""
    PlantHydraulics.flux_out_roots(
        re::PrognosticSoilPressure{FT},
        model::PlantHydraulics.PlantHydraulicsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}

An extension of the `PlantHydraulics.flux_out_roots` function,
 which returns the
net flux of water between the
roots and the soil, when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flux of water between
roots and soil at each soil layer.
"""
function PlantHydraulics.flux_out_roots(
    re::PrognosticSoilPressure{FT},
    model::PlantHydraulics.PlantHydraulicsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
    return sum(p.root_extraction)
end


"""
   PrognosticRootFlux{FT} <: PlantHydraulics.AbstractRootExtraction{FT}

Concrete type of PlantHydraulics.AbstractRootExtraction, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Soil.RootExtraction`:both 
are used at the same time,
ensuring that the water flux into the roots is extracted correctly
from the soil.
"""
struct PrognosticRootFlux{FT} <: PlantHydraulics.AbstractRootExtraction{FT} end

"""
    PlantHydraulics.flux_out_roots(
        re::PrognosticRootFlux{FT},
        t::FT,
    )::FT where {FT}

An extension of the `PlantHydraulics.flux_out_roots` function,
 which returns the
net flux of water between the
roots and the soil, when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flux of water between
roots and soil at each soil layer.
"""
function PlantHydraulics.flux_out_roots(
    re::PrognosticRootFlux{FT},
)::FT where {FT}
    return re
end


"""
    RootExtraction{FT} <: Soil.AbstractSoilSource{FT}

Concrete type of Soil.AbstractSoilSource, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Vegetation.PrognosticSoilPressure`:both 
are used at the same time,
ensuring that the water flux into the roots is extracted correctly
from the soil.
"""
struct RootExtraction{FT} <: Soil.AbstractSoilSource{FT} end

"""
    Soil.source!(dY::ClimaCore.Fields.FieldVector,
                          src::RootExtraction{FT},
                          Y::ClimaCore.Fields.FieldVector,
                          p::ClimaCore.Fields.FieldVector)::ClimaCore.Fields.Field  where {FT}

An extension of the `Soil.source!` function,
 which computes source terms for the 
soil model; this method returns the water loss or gain due
to roots when a plant hydraulic prognostic model is included.
"""
function Soil.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::RootExtraction{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
)::ClimaCore.Fields.Field where {FT}
    return dY.soil.ϑ_l .+= p.root_extraction
end
