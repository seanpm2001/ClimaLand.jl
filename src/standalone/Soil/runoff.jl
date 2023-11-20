export soil_surface_infiltration

"""
    AbstractRunoffModel

The abstract type for soil runoff models to be used with
the following boundary condition types:
- `ClimaLSM.Soil.AtmosDrivenFluxBC`
- `ClimaLSM.Soil.RichardsAtmosDrivenFluxBC`,

and for these functions:
-`ClimaLSM.Soil.soil_surface_infiltration`
- `ClimaLSM.Soil.subsurface_runoff_source`
- `ClimaLSM.source!`.

Please see the documentation for these for more details.
The model should specify the subsurface runoff sink term as well
as the surface runoff implementation.
"""
abstract type AbstractRunoffModel end

"""
    NoRunoff <: AbstractRunoffModel

A concrete type of soil runoff model; the 
default choice which does not include the 
effects of runoff.
"""
struct NoRunoff <: AbstractRunoffModel end

"""
    soil_surface_infiltration(::NoRunoff, net_water_flux, _...)

A function which computes the infiltration into the soil
 for the default of `NoRunoff`.

If `net_water_flux = P+E`, where `P` is the precipitation and `E`
is the evaporation (both negative if towards the soil), 
this returns `P+E` as the water boundary flux for the soil.
"""
soil_surface_infiltration(::NoRunoff, net_water_flux, _...) = net_water_flux

"""
    subsurface_runoff_source(runoff::AbstractRunoffModel)::Union{Nothing, AbstractSoilSource} 

A function which returns the soil source for the runoff model 
`runoff`; the default returns nothing in which case no source is added.

"""
subsurface_runoff_source(
    runoff::AbstractRunoffModel,
)::Union{Nothing, AbstractSoilSource} = nothing

"""
    append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple

Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple = (srcs..., src)

"""
    append_source(src::Nothing , srcs::Tuple)::Tuple

Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::Nothing, srcs::Tuple)::Tuple = srcs


# Surface runoff at a point

"""
    PointSurfaceRunoff <: AbstractRunoffModel

A concrete type of soil runoff model, which only
models surface runoff at a point based on the soil
surface moisture.
"""
struct PointSurfaceRunoff <: AbstractRunoffModel end

"""
    soil_surface_infiltration(::PointSurfaceRunoff, net_water_flux, Y, p, _)

A function which computes the infiltration into the soil
in the case of `SF`.

If `net_water_flux = P+E`, where `P` is the precipitation and `E`
is the evaporation (both negative if towards the soil), 
this returns `P+E` as the water boundary flux for the soil.
"""
function soil_surface_infiltration(::PointSurfaceRunoff, net_water_flux, Y, p, _) where {FT <: AbstractFloat}
    i_c = soil_infiltration_capacity(Y, p) # field
    function infl_at_point(i_c, F)
        if i_c < 0
            # if i_c < 0 and smaller magnitude than net_water_flux<0, this returns i_c. Runoff = i_c - F
            # if i_c < 0 and larger magnitude than net_water_flux<0, this returns net_water_flux. Runoff = 0
            # if i_c < 0 and net_water_flux > 0, this returns net_water_flux. Runoff = 0
            return max(i_c,F)
        else
            if F < 0
                #if i_c > 0 and net_water_flux < 0, this returns i_c. Runoff = i_c - F
                return i_c
            else
                #if i_c >0 and net_water_flux > 0, this should return the sum. Runoff = i_c
                return i_c+F
            end
        end
    end
    return infl_at_point.(i_c, net_water_flux)
    #in all cases, runoff = boundary_flux - F
end

        


"""
    function soil_infiltration_capacity(
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
    )

Computes the infiltration capacity of the soil as
-K_sfc [ψ(ν) + δz- ψ(θ_sfc)]/δz. 

Note that if θ_sfc > ν,
the infiltration capacity is positive.
"""
function soil_infiltration_capacity(Y::ClimaCore.Fields.FieldVector, p::NamedTuple)
    FT = eltype(Y.soil.ϑ_l)
    z = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z
    Δz, _ = ClimaLSM.get_Δz(z)
    p_len = Spaces.nlevels(axes(p.soil.K))
    K_c = Fields.level(p.soil.K, p_len)
    ψ_c = Fields.level(p.soil.ψ, p_len)
    ψ_bc = FT(0)
    return ClimaCore.Fields.Field(ClimaLSM.diffusive_flux(K_c, ψ_bc .+ Δz, ψ_c, Δz), axes(Δz))
end
