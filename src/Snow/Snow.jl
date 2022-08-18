module Snow
using UnPack
using DocStringExtensions
using SurfaceFluxes
using SurfaceFluxes.UniversalFunctions
using Thermodynamics
using StaticArrays
using ClimaCore.Fields: coordinate_field
using ClimaLSM
using Flux
import ..Parameters as LSMP
import ClimaLSM.Domains: coordinates
import ClimaLSM:
    AbstractModel,
    make_update_aux,
    make_rhs,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types
export SnowModelParameters,
    SnowModel,
    surface_fluxes

abstract type AbstractSnowModel{FT} <: AbstractModel{FT} end
"""
     AbstractAtmosphericDrivers{FT <: AbstractFloat}

An abstract type of atmospheric drivers of the bucket model.
"""
abstract type AbstractAtmosphericDrivers{FT <: AbstractFloat} end

abstract type AlbedoModel{FT} end
abstract type EmissivityModel{FT} end
abstract type TauModel{FT} end
abstract type BucketModel{FT} end

ClimaLSM.name(::AbstractSnowModel) = :snow

struct ConstantAlbedo{FT} <: AlbedoModel{FT}
    α::FT
end

struct DynamicAlbedo{FT} <: AlbedoModel{FT}
    p1::FT
    p2::FT
    #etc. for other parameters
end

struct ConstantEmissivity{FT} <: EmissivityModel{FT}
    ϵ::FT
end

struct ConstantTau{FT} <: TauModel{FT}
    τ::FT
end

struct ConstantBucket{FT} <: BucketModel{FT}
    θ_max::FT
end

struct ModelParams{FT} #I imagine this will be edited to just lump CliMA parameters in
    σ::FT
    Lvap::FT
    Lfus::FT
    cp_w::FT
    cp_i::FT
    T_F::FT
    ρ_l::FT
end

mutable struct SnowState{FT}
    e_int::FT
    swe::FT
    z::FT
    has_snow::Bool
end

struct Forcing{FT}
    P::FT
    Gsol::FT
    T_air::FT
    T_soil::FT
    q_air::FT
    v_wind::FT
end

mutable struct DerivedParams{FT}
    ρ_s::FT
    ρ_a::FT
    cp_s::FT
    cp_a::FT
    T̄::FT
    Tsurf::FT
    g_aw::FT
    g_as::FT
    qw::FT
    qi::FT
    l::FT
    κ_s::FT
    κ_a::FT
    α::FT
    ϵ::FT
    τ::FT
    θ_max::FT
    imax::FT
end

mutable struct FluxTerms{FT}
    P̃::FT
    Ẽ::FT
    R̃::FT
    G::FT
    Qconv::FT
    Qlat::FT
    QR::FT
end

struct SnowModel{FT}
    α_model::AlbedoModel{FT}
    ϵ_model::EmissivityModel{FT}
    τ_model::TauModel{FT}
    θ_model::BucketModel{FT}
    state::SnowState{FT}
    params::ModelParams{FT}
    derived::DerivedParams{FT}
    flux::FluxTerms{FT}
    dt::FT
end

struct SnowModelParameters{
    FT <: AbstractFloat,
    AAM <: AlbedoModel,
    PSE,
}
    κ_snow::FT
    albedo::AAM
    z_0m::FT
    z_0b::FT
    "Earth Parameter set; physical constants, etc"
    earth_param_set::PSE
end

SnowModelParameters(
    κ_snow::FT,
    albedo::AAM,
    dz_dt_neural_net::NN
    z_0m::FT,
    z_0b::FT,
    earth_param_set::PSE,
) where {FT, AAM, PSE} = SnowModelParameters{FT, AAM, PSE}(
    κ_snow,
    albedo,
    z_0m,
    z_0b,
    earth_param_set,
)


struct PrescribedAtmosphere{FT, LP, TA, UA, QA, R} <:
       AbstractAtmosphericDrivers{FT}
    "Precipitation (m/s) function of time: positive by definition"
    liquid_precip::LP
    "Prescribed atmospheric temperature (function of time)  at the reference height (K)"
    T_atmos::TA
    "Prescribed wind speed (function of time)  at the reference height (m/s)"
    u_atmos::UA
    "Prescribed relative humidity (function of time)  at the reference height (_)"
    rel_q_atmos::QA
    "Reference height, relative to surface elevation(m)"
    h_atmos::FT
    "Surface air density (kg/m^3). Note that convergence errors result if ρ_sfc = ρ_atmos."
    ρ_sfc::FT # Eventually, computed from mean surface pressure (from Atmos) + land T_sfc
    "Downward radiation function of time (W/m^2): positive indicates towards surface"
    R_d::R
end

function PrescribedAtmosphere(
    precipitation,
    T_atmos,
    u_atmos,
    rel_q_atmos,
    ρ_atmos,
    h_atmos,
    ρ_sfc,
    R_d,
)
    args = (precipitation, T_atmos, u_atmos, rel_q_atmos, ρ_atmos, R_d)
    PrescribedAtmosphere{typeof(h_atmos), typeof.(args)...}(
        args...,
        h_atmos,
        ρ_sfc,
    )
end

struct SnowModel{
    FT,
    PS <: SnowModelParameters{FT},
    ATM <: AbstractAtmosphericDrivers{FT},
    D,
} <: AbstractSnowModel{FT}
    "Parameters required by the bucket model"
    parameters::PS
    "The atmospheric drivers: Prescribed or Coupled"
    atmos::ATM
    "The domain of the model"
    domain::D
end

function SnowModel(;
    parameters::SnowModelParameters{FT, PSE},
    domain::ClimaLSM.Domains.AbstractDomain,
    atmosphere::ATM,
) where {FT, PSE, ATM}
    args = (parameters, atmosphere, domain)
    SnowModel{FT, typeof.(args)...}(args...)
end

prognostic_types(::SnowModel{FT}) where {FT} = (FT, FT, FT)
prognostic_vars(::SnowModel) = (:swe, :z, :e_int)
auxiliary_types(::SnowModel{FT}) where {FT} = (FT, FT, FT, FT, FT)
auxiliary_vars(::SnowModel) = (:q_sfc, :E, :LHF, :SHF, :R_n)

"""
    surface_fluxes(Y,p,
                    t::FT,
                    parameters::P,
                    atmos::PA,
                    radiation::PR,
                    ) where {FT <: AbstractFloat, P <: BucketModelParameters{FT},  PA <: PrescribedAtmosphere{FT}, PR <: PrescribedRadiativeFluxes{FT}}

Computes the surface flux terms at the ground for a standalone simulation:
net radiation,  SHF,  LHF,
as well as the water vapor flux (in units of m^3/m^2/s of water).
Positive fluxes indicate flow from the ground to the atmosphere.

It solves for these given atmospheric conditions, stored in `atmos`,
 downwelling shortwave and longwave radiation (in `radiation`),
model parameters, and the surface temperature and surface specific
humidity.

Currently, we only support soil covered surfaces.
"""
function surface_fluxes(
    Y,
    p,
    t::FT,
    parameters::P,
    atmos::PA,
    radiation::PR,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PA <: PrescribedAtmosphere{FT},
    PR <: PrescribedRadiativeFluxes{FT},
}

    return surface_fluxes_at_a_point.(
        Y.bucket.T_sfc,
        p.bucket.q_sfc,
        Y.bucket.S,
        coordinate_field(Y.bucket.S),
        t,
        Ref(parameters),
        Ref(atmos),
        Ref(radiation),
    )
end

function surface_fluxes_at_a_point(
    T_sfc::FT,
    q_sfc::FT,
    S::FT,
    coords,
    t::FT,
    parameters::P,
    atmos::PA,
    radiation::PR,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PA <: PrescribedAtmosphere{FT},
    PR <: PrescribedRadiativeFluxes{FT},
}
    @unpack ρ_atmos, T_atmos, u_atmos, rel_q_atmos, h_atmos, ρ_sfc = atmos
    @unpack albedo, z_0m, z_0b, S_c, earth_param_set = parameters
    @unpack LW_d, SW_d = radiation
    _σ = LSMP.Stefan(earth_param_set)
    _ρ_liq = LSMP.ρ_cloud_liq(earth_param_set)

    thermo_params = LSMP.thermodynamic_parameters(earth_param_set)

    # call surface fluxes for E, SHF, LHF
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)
    T_air = T_atmos(t)
    P_atmos = some_function()
    q_atmos = rel_q_atmos(t) * q_sat(P_atmos, T_atmos)# look up this function

    ts_in = Thermodynamics.PhaseEquil_PTq(
        thermo_params,
        P_atmos,
        T_air,
        q_atmos,
    )

    # h_atmos is relative to surface height, so we can set surface height to zero.
    state_sfc = SurfaceFluxes.SurfaceValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
    state_in = SurfaceFluxes.InteriorValues(
        h_atmos,
        SVector{2, FT}(u_atmos(t), 0),
        ts_in,
    )

    # State containers
    sc = SurfaceFluxes.ValuesOnly{FT}(;
        state_in,
        state_sfc,
        z0m = z_0m,
        z0b = z_0b,
    )
    surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)
    conditions = SurfaceFluxes.surface_conditions(surface_flux_params, sc)

    α = surface_albedo(albedo)
    # Recall that the user passed the LW and SW downwelling radiation,
    # where positive values indicate toward surface, so we need a negative sign out front
    # in order to inidicate positive R_n  = towards atmos.
    R_n = -((FT(1) - α) * SW_d(t) + LW_d(t) - _σ * T_sfc^FT(4.0))
    # Land needs a volume flux of water, not mass flux
    E =
        SurfaceFluxes.evaporation(surface_flux_params, sc, conditions.Ch) /
        _ρ_liq
    return (R_n = R_n, LHF = conditions.lhf, SHF = conditions.shf, E = E)
end



"""
    make_rhs(model::BucketModel{FT}) where {FT}

Creates the rhs! function for the bucket model.
"""
function make_rhs(model::BucketModel{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        @unpack d_soil, T0, κ_soil, ρc_soil, S_c, W_f = model.parameters
        # Always positive
        @unpack SHF, LHF, R_n, E = p.snow
       
        dY.snow.z .= neural net # Equation (4a) of the text. 
        dY.snow.swe .=
        dY.snow.e_int .= 
    end
    return rhs!
end


"""
    make_update_aux(model::BucketModel{FT}) where {FT}

Creates the update_aux! function for the BucketModel.
"""
function make_update_aux(model::BucketModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        ρ_sfc = surface_air_density(p, model.atmos)
        p.snow.q_sfc .= # computed somehow

        fluxes = surface_fluxes(
            Y,
            p,
            t,
            model.parameters,
            model.atmos,
            model.radiation,
        )
        @. p.snow.LHF = fluxes.LHF
        @. p.snow.SHF = fluxes.SHF
        @. p.snow.R_n = fluxes.R_n
        @. p.snow.E = fluxes.E
    end
    return update_aux!
end

end