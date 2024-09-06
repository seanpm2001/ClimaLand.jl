module Snow

using DocStringExtensions
using Dates
import ...Parameters as LP
using ClimaCore
using Thermodynamics
using ClimaLand
using ClimaLand:
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    turbulent_fluxes,
    net_radiation,
    AbstractModel,
    heaviside

import ClimaLand:
    make_update_aux,
    make_compute_exp_tendency,
    make_update_boundary_fluxes,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    prognostic_domain_names,
    auxiliary_types,
    auxiliary_domain_names,
    surface_temperature,
    surface_specific_humidity,
    surface_height,
    surface_albedo,
    surface_emissivity,
    get_drivers
export SnowParameters, SnowModel

"""
    AbstractSnowModel{FT} <: ClimaLand.AbstractExpModel{FT}

Defines a new type of abstract explicit model for snow modeling.
Currently, the only supported concrete example is called `SnowModel`
and is used as a bulk snow model.
"""
abstract type AbstractSnowModel{FT} <: ClimaLand.AbstractExpModel{FT} end

"""
    AbstractDensityModel{FT}

Defines the model type for density and depth parameterizations
for use within an `AbstractSnowModel` type. Current examples include the
`ConstantDensityModel` and the `Anderson1976` models.
"""
abstract type AbstractDensityModel{FT <: AbstractFloat} end

"""
    ConstantDensityModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}

Establishes the density parameterization where snow density
is always treated as a constant (type FT).
"""
struct ConstantDensityModel{FT} <: AbstractDensityModel{FT}
    ρ_snow::FT
end

"""
    ConstantDensityModel{FT}(ρ::FT)

An outer constructor for the `ConstantDensityModel` density parameterization for usage in a snow model.
"""
function ConstantDensityModel{FT}(ρ::FT) where {FT}
    return ConstantDensityModel{FT}(ρ)
end


"""
    ConstantDensityModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}

Establishes the density parameterization where snow density
compacts according to the seminal works of Anderson in the 1970s (and specifically the
numerical implementation documented in the Snow17 model).
"""
struct Anderson1976{FT} <: AbstractDensityModel{FT}
    c1::FT
    c2::FT
    c3::FT
    c4::FT
    c5::FT
    cx::FT
    ρ_d::FT
end

"""
    Anderson1976(FT; c1, c2, c3, c4, c5, ρ_d, cx)

An outer constructor for the `Anderson1976` density parameterization for usage in a snow model.
Uses the default constants defined by Anderson but allows for redefinition of the model constants.
These include:
- c1: fractional increase in density (0.026 cm/hr)
- c2: compression constant estimated by Kojima in 1967 (21 cm^3/g)
- c3: fractional settling rate at 0 degrees C for densities less than the critical density ρ_d (0.005 1/hr)
- c4: constant (0.1 1/degrees C)
- c5: scaling of the settling rate when no water is present in the snowpack (0, it is 2 when water is present)
- cx: destructive metamorphism decay factor for densities greater than the critical density ρ_d (23)
- ρ_d: the critical density (unitless as a fraction of the density of liquid water, 0.15)
"""
function Anderson1976(FT::DataType;
    c1 = 0.026,
    c2 = 21.0,
    c3 = 0.005,
    c4 = 0.10,
    c5 = 0.0, #Snow17 code says 1.0, paper says 0.0?
    ρ_d = 0.15,
    cx = 23.0,
)
return Anderson1976{FT}(FT(c1), FT(c2), FT(c3), FT(c4), FT(c5), FT(ρ_d), FT(cx))
end

"""
    SnowParameters{FT <: AbstractFloat, PSE}

A struct for storing parameters of the `SnowModel`.

Note that in our current implementation of runoff, a physical
timescale is required and computed using Ksat and the depth
of the snow. For shallow snowpacks, this will fall below the timestep
of the model. For that reason, we pass the timestep of the model as 
a parameter, and take the larger of the timestep and the physical timescale
as the value used in the model. Future implementations will revisit this.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SnowParameters{FT <: AbstractFloat, DM <: AbstractDensityModel, PSE}
    "Choice of parameterization for snow density"
    density::DM
    "Roughness length over snow for momentum (m)"
    z_0m::FT
    "Roughness length over snow for scalars (m)"
    z_0b::FT
    "Albedo of snow (unitless)"
    α_snow::FT
    "Emissivity of snow (unitless)"
    ϵ_snow::FT
    "Volumetric holding capacity of water in snow (unitless)"
    θ_r::FT
    "Hydraulic conductivity of wet snow (m/s)"
    Ksat::FT
    "Thermal conductivity of ice (W/m/K)"
    κ_ice::FT
    "Timestep of the model (s)"
    Δt::FT
    "Areal specific heat of ground interacting with snow (J/m^2/K)"
    ρcD_g::FT
    "Clima-wide parameters"
    earth_param_set::PSE
end

"""
   SnowParameters{FT}(Δt;
                      density = ConstantDensityModel(FT(200)),
                      z_0m = FT(0.0024),
                      z_0b = FT(0.00024),
                      α_snow = FT(0.8),
                      ϵ_snow = FT(0.99),
                      θ_r = FT(0.08),
                      Ksat = FT(1e-3),
                      κ_ice = FT(2.21),
                      ρcD_g = FT(3.553e5),
                      earth_param_set::PSE) where {FT, PSE}

An outer constructor for `SnowParameters` which supplies defaults for
all arguments but `earth_param_set`.
"""
function SnowParameters{FT}(
    Δt;
    density::DM = ConstantDensityModel(FT(200)),
    z_0m = FT(0.0024),
    z_0b = FT(0.00024),
    α_snow = FT(0.8),
    ϵ_snow = FT(0.99),
    θ_r = FT(0.08),
    Ksat = FT(1e-3),
    κ_ice = FT(2.21),
    ρcD_g = FT(3.553e5),
    earth_param_set::PSE,
) where {FT <: AbstractFloat, DM <: AbstractDensityModel, PSE}
    return SnowParameters{FT, DM, PSE}(
        density,
        z_0m,
        z_0b,
        α_snow,
        ϵ_snow,
        θ_r,
        Ksat,
        κ_ice,
        Δt,
        ρcD_g,
        earth_param_set,
    )
end

Base.broadcastable(ps::SnowParameters) = tuple(ps)

"""
    struct SnowModel{
        FT,
        PS <: SnowParameters{FT},
        ATM <: AbstractAtmosphericDrivers{FT},
        RAD <: AbstractRadiativeDrivers{FT},
        D,
    } <: AbstractSnowModel{FT}

A container/type for the bulk snow model, based on the UEB snow model
of Tarboton et al. (1995) and Tarboton and Luce (1996).
"""
struct SnowModel{
    FT,
    PS <: SnowParameters{FT},
    ATM <: AbstractAtmosphericDrivers{FT},
    RAD <: AbstractRadiativeDrivers{FT},
    D,
} <: AbstractSnowModel{FT}
    "Parameters required by the snow model"
    parameters::PS
    "The atmospheric drivers: Prescribed or Coupled"
    atmos::ATM
    "The radiation drivers: Prescribed or Coupled"
    radiation::RAD
    "The domain of the model"
    domain::D
end

function SnowModel(;
    parameters::SnowParameters{FT, DM, PSE},
    domain::ClimaLand.Domains.AbstractDomain,
    atmos::ATM,
    radiation::RAD,
) where {FT, DM, PSE, ATM, RAD}
    args = (parameters, atmos, radiation, domain)
    SnowModel{FT, typeof.(args)...}(args...)
end

"""
    prognostic_vars(::SnowModel)

Returns the prognostic variable names of the snow model.

For this model, we track the snow water equivalent S [m] and
the energy per unit area U [J/m^2] prognostically, as well as snow depth [m].
"""
prognostic_vars(::SnowModel) = (:S, :U, :Z)

"""
    prognostic_types(::SnowModel{FT})

Returns the prognostic variable types of the snow model;
both snow water equivalent and energy per unit area
are scalars.
"""
prognostic_types(::SnowModel{FT}) where {FT} = (FT, FT, FT)

"""
    prognostic_domain_names(::SnowModel)

Returns the prognostic variable domain names of the snow model;
both snow water equivalent and energy per unit area
are modeling only as a function of (x,y), and not as a function
of depth. Therefore their domain name is ":surface".
"""
prognostic_domain_names(::SnowModel) = (:surface, :surface, :surface)

"""
    auxiliary_vars(::SnowModel)

Returns the auxiliary variable names for the snow model. These
include the bulk snow density (`ρ_snow`, kg/m^3),
mass fraction in liquid water (`q_l`, unitless),
the bulk temperature (`T`, K), the surface temperature (`T_sfc`, K),
the SHF, LHF, and vapor flux (`turbulent_fluxes.shf`, etc),
the net radiation (`R_n, J/m^2/s)`, the energy flux in liquid water runoff
(`energy_runoff`, J/m^2/s), the water volume in runoff (`water_runoff`, m/s), and the total energy and water fluxes applied to the snowpack.

Since the snow can melt completely in one timestep, we clip the water and energy fluxes
such that SWE cannot become negative and U cannot become unphysical. The
clipped values are what are actually applied as boundary fluxes, and are stored in
`applied_` fluxes.
"""
auxiliary_vars(::SnowModel) = (
    :ρ_snow,
    :q_l,
    :T,
    :T_sfc,
    :turbulent_fluxes,
    :R_n,
    :energy_runoff,
    :water_runoff,
    :total_energy_flux,
    :total_water_flux,
    :applied_energy_flux,
    :applied_water_flux,
    :snow_cover_fraction,
)

auxiliary_types(::SnowModel{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    NamedTuple{(:lhf, :shf, :vapor_flux, :r_ae), Tuple{FT, FT, FT, FT}},
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
)

auxiliary_domain_names(::SnowModel) = (
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


ClimaLand.name(::SnowModel) = :snow

"""
    scf(x::FT; α = FT(1e-3))::FT where {FT}

Returns the snow cover fraction, assuming it is a heaviside
function at 1e-3 meters.

In the future we can play around with other forms, this makes it a bit more
explicit.
"""
function scf(x::FT; α = FT(1e-3))::FT where {FT}
    return heaviside(x - α)
end

function ClimaLand.make_update_aux(model::SnowModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        parameters = model.parameters

        p.snow.ρ_snow .= snow_density.(Ref(parameters.density), Y.snow.S, Y.snow.Z, parameters) #pass density model explicitly for multiple dispatch

        @. p.snow.q_l =
            snow_liquid_mass_fraction(Y.snow.U, Y.snow.S, parameters)

        @. p.snow.T =
            snow_bulk_temperature(Y.snow.U, Y.snow.S, p.snow.q_l, parameters)

        @. p.snow.T_sfc = snow_surface_temperature(p.snow.T)

        @. p.snow.water_runoff =
            compute_water_runoff(Y.snow.S, p.snow.q_l, p.snow.T, p.snow.ρ_snow, Y.snow.Z, parameters)

        @. p.snow.energy_runoff =
            p.snow.water_runoff * volumetric_internal_energy_liq(FT, parameters)
        @. p.snow.snow_cover_fraction = scf(Y.snow.S)
    end
end

function ClimaLand.make_update_boundary_fluxes(model::SnowModel{FT}) where {FT}
    function update_boundary_fluxes!(p, Y, t)
        p.snow.turbulent_fluxes .= turbulent_fluxes(model.atmos, model, Y, p, t)
        p.snow.R_n .= net_radiation(model.radiation, model, Y, p, t)
        # How does rain affect the below?
        P_snow = p.drivers.P_snow

        # Heaviside function because we cannot change the water content of the snow by
        # sublimation, evap, melting, or rain
        # unless there is already snow on the ground. 
        @. p.snow.total_water_flux =
            -P_snow +
            (-p.snow.turbulent_fluxes.vapor_flux + p.snow.water_runoff) *
            heaviside(Y.snow.S - eps(FT))

        # I think we want dU/dt to include energy of falling snow.
        # otherwise snow can fall but energy wont change
        # We are assuming that the sensible heat portion of snow is negligible.
        _LH_f0 = FT(LP.LH_f0(model.parameters.earth_param_set))
        _ρ_liq = FT(LP.ρ_cloud_liq(model.parameters.earth_param_set))
        ρe_falling_snow = -_LH_f0 * _ρ_liq # per unit vol of liquid water

        @. p.snow.total_energy_flux =
            -P_snow * ρe_falling_snow +
            (
                -p.snow.turbulent_fluxes.lhf - p.snow.turbulent_fluxes.shf -
                p.snow.R_n + p.snow.energy_runoff
            ) * heaviside(Y.snow.S - eps(FT))

        @. p.snow.applied_water_flux =
            clip_dSdt(Y.snow.S, p.snow.total_water_flux, model.parameters.Δt)
        @. p.snow.applied_energy_flux = clip_dUdt(
            Y.snow.U,
            Y.snow.S,
            p.snow.total_energy_flux,
            p.snow.total_water_flux,
            model.parameters.Δt,
        )
    end
end

function ClimaLand.make_compute_exp_tendency(model::SnowModel{FT}) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        @. dY.snow.S = p.snow.applied_water_flux
        @. dY.snow.U = p.snow.applied_energy_flux
        dY.snow.Z .= clip_dZdt.(
            Y.snow.S,
            Y.snow.Z,
            p.snow.applied_water_flux,
            dzdt(model.parameters.density, model, Y, p, t),
            model.parameters.Δt
        )
    end
    return compute_exp_tendency!
end

"""
    clip_dSdt(S, dSdt, Δt)

A helper function which clips the tendency of S such that 
S will not become negative.
"""
function clip_dSdt(S, dSdt, Δt)
    if S + dSdt * Δt < 0
        return -S / Δt
    else
        return dSdt
    end
end

"""
    clip_dUdt(U, S, dUdt, dSdt, Δt)

A helper function which clips the tendency of U such that 
U will not become positive, and which ensures that if S
goes to zero in a step, U will too.
"""
function clip_dUdt(U, S, dUdt, dSdt, Δt)
    if (U + dUdt * Δt) > 0
        return -U / Δt
    elseif S + dSdt * Δt < 0
        return -U / Δt
    else
        return dUdt
    end
end

"""
    clip_dZdt(U, S, dUdt, dSdt, Δt)

A helper function which clips the tendency of Z such that
its behavior is consistent with that of S.
"""
function clip_dZdt(S, Z, dSdt, dZdt, Δt)
    if (S + dSdt * Δt) < 0
        return -Z / Δt
    elseif (Z + dZdt * Δt) < (S + dSdt * Δt)
        #do we need to do something about setting q_l = 1 here?
        return (S-Z) / Δt + dSdt
    else
        return dZdt
    end
end

"""
    ClimaLand.get_drivers(model::SnowModel)

Returns the driver variable symbols for the SnowModel.
"""
function ClimaLand.get_drivers(model::SnowModel)
    return (model.atmos, model.radiation)
end

include("./snow_parameterizations.jl")
end
