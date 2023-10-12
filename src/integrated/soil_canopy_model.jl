export SoilCanopyModel
"""
    struct SoilCanopyModel{
        FT,
        SM <: Soil.EnergyHydrology{FT},
        VM <: Canopy.CanopyModel{FT},
    } <: AbstractLandModel{FT}
        soil::SM
        canopy::VM
    end

A concrete type of land model used for simulating systems with a 
canopy and a soil component.
$(DocStringExtensions.FIELDS)
"""
struct SoilCanopyModel{
    FT,
    SM <: Soil.EnergyHydrology{FT},
    VM <: Canopy.CanopyModel{FT},
} <: AbstractLandModel{FT}
    "The soil model to be used"
    soil::SM
    "The canopy model to be used"
    canopy::VM
end

"""
    CanopyRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT}

A struct used to compute radiative fluxes in land surface models,
indicating that 
canopy absorption and emission is taken into account when computing
radiation at the surface of the soil or snow.

The only other alternative at this stage is
ClimaLSM.PrescribedRadiativeFluxes, where the prescribed downwelling
short and longwave radiative fluxes are used directly,
without accounting for the canopy. There is a different method
of the function `soil_boundary_fluxes` in this case. 
"""
struct CanopyRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT} end

"""
    soil_boundary_fluxes(
        bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:CanopyRadiativeFluxes},
        boundary::ClimaLSM.TopBoundary,
        model::EnergyHydrology{FT},
        Y,
        Δz,
        p,
        t,
    ) where {FT}

A method of `ClimaLSM.Soil.soil_boundary_fluxes` which is used for
integrated land surface models; this computes and returns the net
energy and water flux at the surface of the soil for use as boundary
conditions.

The net radiative, sensible heat, latent heat, and evaporative fluxes 
are computed and stored in the auxiliary state of the integrated land 
surface models, and this function simply returns those. They are updated
each time step in `update_aux`.
"""
function soil_boundary_fluxes(
    bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:CanopyRadiativeFluxes},
    boundary::ClimaLSM.TopBoundary,
    model::EnergyHydrology{FT},
    Y,
    Δz,
    p,
    t,
) where {FT}
    infiltration = soil_surface_infiltration(
        bc.runoff,
        bc.atmos.liquid_precip(t) .+ p.soil_evap,
        Y,
        p,
        t,
        model.parameters,
    )
    G = @. -p.soil_LW_n - p.soil_SW_n + p.soil_lhf + p.soil_shf
    return infiltration, G
end

"""
    SoilCanopyModel{FT}(;
                         land_args::NamedTuple = (;),
                         soil_model_type::Type{SM},
                         soil_args::NamedTuple = (;),
                         canopy_component_types::NamedTuple = (;),
                         canopy_component_args::NamedTuple = (;),
                         canopy_model_args::NamedTuple = (;),
                         ) where {FT, SM <: Soil.EnergyHydrology{FT}}
A constructor for the `SoilCanopyModel`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `SoilCanopyModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function SoilCanopyModel{FT}(;
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    canopy_component_types::NamedTuple = (;),
    canopy_component_args::NamedTuple = (;),
    canopy_model_args::NamedTuple = (;),
) where {FT, SM <: Soil.EnergyHydrology{FT}}

    # These may be passed in, or set, depending on use scenario.
    (; atmos, radiation) = land_args
    # These should always be set by the constructor.
    Δz = minimum(
        ClimaCore.Fields.Δz_field(
            ClimaLSM.coordinates(soil_args.domain).subsurface,
        ),
    )
    sources = (RootExtraction{FT}(), Soil.PhaseChange{FT}(Δz))
    # add heat BC
    top_bc = ClimaLSM.Soil.AtmosDrivenFluxBC(atmos, CanopyRadiativeFluxes{FT}())
    zero_flux = FluxBC((p, t) -> eltype(t)(0.0))
    boundary_conditions = (;
        top = top_bc,
        bottom = (water = Soil.FreeDrainage(), heat = zero_flux),
    )
    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )

    transpiration = Canopy.PlantHydraulics.DiagnosticTranspiration{FT}()

    soil_driver = PrognosticSoil(;
        soil_α_PAR = soil.parameters.PAR_albedo,
        soil_α_NIR = soil.parameters.NIR_albedo,
    )

    canopy = Canopy.CanopyModel{FT}(;
        autotrophic_respiration = canopy_component_types.autotrophic_respiration(
            canopy_component_args.autotrophic_respiration...,
        ),
        radiative_transfer = canopy_component_types.radiative_transfer(
            canopy_component_args.radiative_transfer...,
        ),
        photosynthesis = canopy_component_types.photosynthesis(
            canopy_component_args.photosynthesis...,
        ),
        conductance = canopy_component_types.conductance(
            canopy_component_args.conductance...,
        ),
        hydraulics = canopy_component_types.hydraulics(;
            transpiration = transpiration,
            canopy_component_args.hydraulics...,
        ),
        soil_driver = soil_driver,
        atmos = atmos,
        radiation = radiation,
        canopy_model_args...,
    )

    return SoilCanopyModel{FT, typeof(soil), typeof(canopy)}(soil, canopy)
end

"""
    interaction_vars(m::SoilCanopyModel)

The names of the additional auxiliary variables that are 
included in the integrated Soil-Canopy model.
"""
interaction_vars(m::SoilCanopyModel) = (
    :root_extraction,
    :soil_LW_n,
    :soil_SW_n,
    :soil_evap,
    :soil_shf,
    :soil_lhf,
    :T_soil,
    :ρ_soil,
    :q_soil,
    :α_soil,
    :ϵ_soil,
    :q_canopy,
    :canopy_LW_n,
    :canopy_SW_n,
    :LW_out,
    :SW_out,
    :ustar,
    :T_canopy_air,
    :q_canopy_air,
    :shf,
    :lhf
)

"""
    interaction_types(m::SoilCanopyModel)

The types of the additional auxiliary variables that are 
included in the integrated Soil-Canopy model.
"""
interaction_types(m::SoilCanopyModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT)

"""
    interaction_domain_names(m::SoilCanopyModel)

The domain names of the additional auxiliary variables that are 
included in the integrated Soil-Canopy model.
"""
interaction_domain_names(m::SoilCanopyModel) = (
    :subsurface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
)

"""
    make_interactions_update_aux(
        land::SoilCanopyModel{FT, SM, RM},
    ) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Canopy.CanopyModel{FT}}

A method which makes a function; the returned function 
updates the auxiliary variable `p.root_extraction`, which
is needed for a sink term for the soil model and to create the
lower water boundary condition for the canopy model.
It also updates
the soil surface fluxes, which are affected by the presence of a canopy.

This function is called each ode function evaluation.
"""
function make_interactions_update_aux(
    land::SoilCanopyModel{FT, SM, RM},
) where {FT, SM <: Soil.EnergyHydrology{FT}, RM <: Canopy.CanopyModel{FT}}
    function update_aux!(p, Y, t)
        z =
            ClimaCore.Fields.coordinate_field(
                land.soil.domain.space.subsurface,
            ).z
        (; conductivity_model) = land.canopy.hydraulics.parameters
        area_index = p.canopy.hydraulics.area_index
        @. p.root_extraction =
            (
                area_index.root + getproperty(
                    area_index,
                    land.canopy.hydraulics.compartment_labels[1],
                )
            ) / 2 *
            PlantHydraulics.flux(
                z,
                land.canopy.hydraulics.compartment_midpoints[1],
                p.soil.ψ,
                p.canopy.hydraulics.ψ.:1,
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    p.soil.ψ,
                ),
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    p.canopy.hydraulics.ψ.:1,
                ),
            ) *
            (land.canopy.hydraulics.parameters.root_distribution(z))

        earth_param_set = land.soil.parameters.earth_param_set
        p.ϵ_soil .= surface_emissivity(land.soil, Y, p)
        p.T_soil .= ClimaLSM.surface_temperature(land.soil, Y, p, t)
        p.ρ_soil .= ClimaLSM.surface_air_density(land.soil.boundary_conditions.top.atmos, land.soil, Y, p, t, p.T_soil)
        p.q_soil .= ClimaLSM.surface_specific_humidity(land.soil, Y, p, p.T_soil, p.ρ_soil)

        # solve for T_canopy_air, q_canopy_air, ustar; update in place
        canopy_airspace_solve!(p, Y, t, earth_param_set, land, land.canopy.atmos)

        ustar = p.ustar

        _ρ_liq = LSMP.ρ_cloud_liq(earth_param_set)
        thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
        _LH_v0 = LSMP.LH_v0(earth_param_set)
        d_leaf = 0.04
        Cv = 0.01
        Cs_canopy = 0.004
        Cs_bare =@.  0.4/0.13*(0.01*ustar/(1.5e-5))^(-0.45)
        AI = @. p.canopy.hydraulics.area_index.leaf+p.canopy.hydraulics.area_index.stem
        W = @. exp(-AI)
        Cs = @. W*Cs_bare + (1-W)*Cs_canopy
        r_b = @. 1/Cv*(ustar/d_leaf)^-0.5/AI
        r_ah_cs = @. 1/Cs/ustar
        r_aw_cs = r_ah_cs
        r_soil = ClimaLSM.surface_resistance(land.soil, Y, p, t)
        r_canopy = ClimaLSM.surface_resistance(land.canopy, Y, p, t) .+ r_b
        
        ts_in = construct_atmos_ts(land.canopy.atmos, t, thermo_params)
        ρ_air = Thermodynamics.air_density(thermo_params, ts_in)
        cp_m = Thermodynamics.cp_m(thermo_params, ts_in)

        @. p.soil_shf = -ρ_air * cp_m * (p.T_canopy_air - p.T_soil) / r_ah_cs
        @. p.soil_evap = -ρ_air * (p.q_canopy_air - p.q_soil) / (r_aw_cs + r_soil) ./ _ρ_liq
        @. p.soil_lhf = p.soil_evap * _LH_v0 * _ρ_liq 

        # Compute transpiration using T_canopy
        p.q_canopy .=  ClimaLSM.surface_specific_humidity(land.canopy, Y, p, p.T_canopy_air, ρ_air)
        @. p.canopy.energy.shf  = -ρ_air * cp_m * (p.T_canopy_air - p.T_canopy_air) / r_b # equals zero, in this assumption
        @. p.canopy.conductance.transpiration = -ρ_air * (p.q_canopy_air - p.q_canopy) / r_canopy / _ρ_liq
        @. p.canopy.energy.lhf  = p.canopy.conductance.transpiration * _LH_v0 * _ρ_liq


        # Transpiration is per unit ground area, not leaf area (mult by LAI)
        fa = p.canopy.hydraulics.fa
        hydraulics = land.canopy.hydraulics
        n_stem = hydraulics.n_stem
        n_leaf = hydraulics.n_leaf
        i_end = n_stem + n_leaf

        fa.:($i_end) .= PlantHydraulics.transpiration_per_ground_area(
            hydraulics.transpiration,
            Y,
            p,
            t,
        )

        lsm_radiation_update!(
            p,
            land.canopy.radiative_transfer,
            land.canopy,
            land.soil,
            Y,
            t,
        )

    end
    return update_aux!
end


function canopy_airspace_solve!(p, Y, t, earth_param_set, land, atmos)
    FT = Float64
    T_soil = parent(p.T_soil)[1]
    ρ_soil = parent(p.ρ_soil)[1]
    q_soil = parent(p.q_soil)[1]


    earth_param_set = land.soil.parameters.earth_param_set
    _ρ_liq = LSMP.ρ_cloud_liq(earth_param_set)
    thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
    _LH_v0 = LSMP.LH_v0(earth_param_set)
    surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)
    
    u = atmos.u(t)
    h = atmos.h
    ts_in = construct_atmos_ts(atmos, t, thermo_params)
    ρ_air = Thermodynamics.air_density(thermo_params, ts_in)
    cp_m = Thermodynamics.cp_m(thermo_params, ts_in)

    h_sfc = FT(18.5)
    z_0m = FT(0.1)*h_sfc
    z_0b = FT(0.1)*z_0m
    d_sfc = FT(0.67)*h_sfc


    d_leaf = 0.04
    Cv = 0.01
    Cs_canopy = 0.004
    LAI =  parent(p.canopy.hydraulics.area_index.leaf)[1]
    AI = parent(p.canopy.hydraulics.area_index.leaf .+p.canopy.hydraulics.area_index.stem)[1]
    W = @. exp(-AI)
    r_soil = parent(ClimaLSM.surface_resistance(land.soil, Y, p, t))[1]
    r_stomata = parent(ClimaLSM.surface_resistance(land.canopy, Y, p, t))[1] # canopy upscaled
    initial_guess = [0.5 .*(atmos.T(t) .+ T_soil), 0.5 .*(atmos.q(t) .+ q_soil), u/2, 0.0, 0.0]
    function flux_equality(F, x)
        TS = x[1]
        qS = x[2]
        ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_air, TS, qS)
        state_sfc =
            SurfaceFluxes.SurfaceValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
        state_in =
            SurfaceFluxes.InteriorValues(h - d_sfc, SVector{2, FT}(u, 0), ts_in)
        
        sc = SurfaceFluxes.ValuesOnly{FT}(;
                                          state_in,
                                          state_sfc,
                                          z0m = z_0m,
                                          z0b = z_0b,
                                          gustiness = atmos.gustiness,
                                          )

        # Canopy airspace at z0+d and atmos at h
        conditions = SurfaceFluxes.surface_conditions(
            surface_flux_params,
            sc;
            tol_neutral = SFP.cp_d(surface_flux_params) / 100000,
        )
        r_a = 1 / (conditions.Ch * SurfaceFluxes.windspeed(sc))
        ustar = conditions.ustar
        x[3] = ustar
        x[4] = conditions.shf
        x[5] = conditions.lhf
        
        Cs_bare =@.  0.4/0.13*(0.01*ustar/(1.5e-5))^(-0.45)
        Cs = @. W*Cs_bare + (1-W)*Cs_canopy
        r_b = @. 1/Cv*(ustar/d_leaf)^-0.5 # not upscaled for the canopy
        r_ah_cs = @. 1/Cs/ustar
        
        soil_shf = -ρ_air * cp_m * (TS - T_soil) / r_ah_cs
        soil_evap = -ρ_air * (qS - q_soil) / (r_ah_cs + r_soil)

        Tc = TS
        q_canopy =  ClimaLSM.surface_specific_humidity(land.canopy, Y, p, Tc, ρ_air)
        canopy_shf  = -ρ_air * cp_m * (TS - Tc) / (r_b/AI)
        canopy_evap = -ρ_air * (qS - q_canopy) / (r_stomata + r_b/LAI)

        #shf = -ρ_air * cp_m * (atmos.T(t) - TS) / r_a
        #E = -ρ_air *(atmos.q(t) - qS) / r_a
        F[1] = (canopy_evap + soil_evap - conditions.evaporation) * _LH_v0
        F[2] = canopy_shf + soil_shf - conditions.shf
        F[3] = 0.0
        F[4] = 0.0
        F[5] = 0.0
    end
    soln = nlsolve(
        flux_equality,
        initial_guess,
        xtol = 1.0,
    )
    p.T_canopy_air .= soln.zero[1]
    p.q_canopy_air .= soln.zero[2]
    p.ustar .= soln.zero[3]
    p.shf .= soln.zero[4]
    p.lhf .= soln.zero[5]

end
    

    
    
"""
    net_radiation_at_ground(
        canopy_radiation::Canopy.AbstractRadiationModel{FT},
        canopy,
        ground_model::Soil.EnergyHydrology,
        Y,
        p,
        t,
    ) where {FT}


A function which computes the net radiation at the ground surface
give the canopy radiation model.

Returns the correct radiative fluxes for bare ground in the case
where the canopy LAI is zero.
"""
function lsm_radiation_update!(
    p,
    canopy_radiation::Canopy.AbstractRadiationModel{FT},
    canopy,
    ground_model::Soil.EnergyHydrology,
    Y,
    t,
) where {FT}
    radiation = canopy.radiation
    earth_param_set = canopy.parameters.earth_param_set
    _σ = LSMP.Stefan(earth_param_set)
    LW_d::FT = radiation.LW_d(t)
    SW_d::FT = radiation.SW_d(t)
    c = FT(LSMP.light_speed(earth_param_set))
    h = FT(LSMP.planck_constant(earth_param_set))
    N_a = FT(LSMP.avogadro_constant(earth_param_set))
    (; λ_γ_PAR, λ_γ_NIR, ϵ_canopy) = canopy_radiation.parameters
    energy_per_photon_PAR = h * c / λ_γ_PAR
    energy_per_photon_NIR = h * c / λ_γ_NIR
    T_canopy =
        ClimaLSM.Canopy.canopy_temperature(canopy.energy, canopy, Y, p, t)
    T_soil = p.T_soil
    ϵ_soil = p.ϵ_soil
    
    α_soil_PAR = Canopy.ground_albedo_PAR(canopy)
    α_soil_NIR = Canopy.ground_albedo_NIR(canopy)
    # in W/m^2
    PAR = p.canopy.radiative_transfer.par
    NIR = p.canopy.radiative_transfer.nir

    LW_net_canopy = p.canopy_LW_n
    SW_net_canopy = p.canopy_SW_n
    LW_net_soil = p.soil_LW_n
    SW_net_soil = p.soil_SW_n
    LW_out = p.LW_out
    SW_out = p.SW_out

    # in total: INC - OUT = CANOPY_ABS + (1-α_soil)*CANOPY_TRANS
    # SW out  = reflected par + reflected nir
    @. SW_out =
        energy_per_photon_NIR * N_a * p.canopy.radiative_transfer.rnir +
        energy_per_photon_PAR * N_a * p.canopy.radiative_transfer.rpar

    # net canopy = absorbed par + absorbed nir
    @. SW_net_canopy =
        energy_per_photon_NIR * N_a * p.canopy.radiative_transfer.anir +
        energy_per_photon_PAR * N_a * p.canopy.radiative_transfer.apar

    # net soil = (1-α)*trans for par and nir
    @. SW_net_soil =
        energy_per_photon_NIR *
        N_a *
        p.canopy.radiative_transfer.tnir *
        (1 - α_soil_NIR) +
        energy_per_photon_PAR *
        N_a *
        p.canopy.radiative_transfer.tpar *
        (1 - α_soil_PAR)

    LW_d_canopy = @. (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4 # double checked
    LW_u_soil = @. ϵ_soil * _σ * T_soil^4 + (1 - ϵ_soil) * LW_d_canopy # double checked
    @. LW_net_canopy =
        ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_soil # double checked
    @. LW_net_soil = ϵ_soil * LW_d_canopy - ϵ_soil * _σ * T_soil^4 # double checked
    @. LW_out = (1 - ϵ_canopy) * LW_u_soil + ϵ_canopy * _σ * T_canopy^4 # double checked
end

"""
    PlantHydraulics.root_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        s::PrognosticSoil,
        model::Canopy.PlantHydraulics.PlantHydraulicsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t::FT,
    ) where {FT}

An extension of the `PlantHydraulics.root_flux_per_ground_area!` function,
 which returns the
net flux of water between the
roots and the soil, per unit ground area, 
when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flux of water per ground area between
roots and soil at each soil layer.
"""
function PlantHydraulics.root_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    s::PrognosticSoil,
    model::Canopy.PlantHydraulics.PlantHydraulicsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t::FT,
) where {FT}
    fa .= sum(p.root_extraction)
end

"""
    RootExtraction{FT} <: Soil.AbstractSoilSource{FT}

Concrete type of Soil.AbstractSoilSource, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Canopy.PrognosticSoil`:both 
are used at the same time,
ensuring that the water flux into the roots is extracted correctly
from the soil.
"""
struct RootExtraction{FT} <: Soil.AbstractSoilSource{FT} end


"""
    ClimaLSM.source!(dY::ClimaCore.Fields.FieldVector,
                     src::RootExtraction,
                     Y::ClimaCore.Fields.FieldVector,
                     p::NamedTuple
                     model::EnergyHydrology)

An extension of the `ClimaLSM.source!` function,
 which computes source terms for the 
soil model; this method returns the water and energy loss/gain due
to root extraction.
"""
function ClimaLSM.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::RootExtraction,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::EnergyHydrology,
)
    @. dY.soil.ϑ_l += -1 * p.root_extraction
    @. dY.soil.ρe_int +=
        -1 *
        p.root_extraction *
        volumetric_internal_energy_liq(p.soil.T, model.parameters)
    # if flow is negative, towards soil -> soil water increases, add in sign here.
end
