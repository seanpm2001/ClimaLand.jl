import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Plots
using Statistics
using Dates
using Insolation
using StatsBase
using StaticArrays
using Thermodynamics
using SurfaceFluxes
using NLsolve

using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
const FT = Float64
earth_param_set = create_lsm_parameters(FT)
climalsm_dir = pkgdir(ClimaLSM)
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/ozark/ozark_met_drivers_FLUXNET.jl",
    ),
)
include(joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_domain.jl"))
include(joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_domain.jl"))
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_parameters.jl"),
)
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_simulation.jl"),
)
P_air = atmos.P.(seconds)
q_air = atmos.q.(seconds)
T_air = atmos.T.(seconds)
u_air = atmos.u.(seconds)
h = atmos_h
h_sfc = FT(18.5)
z_0m = FT(0.1)*h_sfc
z_0b = FT(0.1)*z_0m
d_sfc = FT(0.67)*h_sfc
_σ = LSMP.Stefan(earth_param_set)
LW_d = radiation.LW_d.(seconds)
SW_d = radiation.SW_d.(seconds)

T_soil = TS
ϑ_soil = SWC
LAI = LAIfunction.(seconds)
d_leaf = 0.04
Ild = 5.0
Cv = 0.01
Cs_canopy = 0.004
_ρ_liq = LSMP.ρ_cloud_liq(earth_param_set)
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
_LH_v0 = LSMP.LH_v0(earth_param_set)
surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)
ref_time = atmos.ref_time


# soil params
soil_params = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = κ_dry,
    κ_sat_frozen = κ_sat_frozen,
    κ_sat_unfrozen = κ_sat_unfrozen,
    ρc_ds = ρc_ds,
    ν = soil_ν,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ν_ss_gravel = ν_ss_gravel,
    hydrology_cm = vanGenuchten(; α = soil_vg_α, n = soil_vg_n),
    K_sat = soil_K_sat,
    S_s = soil_S_s,
    θ_r = θ_r,
    earth_param_set = earth_param_set,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
);



# make canopy model
canopy_component_types = (;
    autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
    radiative_transfer = Canopy.TwoStreamModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
)
# Individual Component arguments
# Set up autotrophic respiration
autotrophic_respiration_args = (;
    parameters = AutotrophicRespirationParameters{FT}(;
        ne = ne,
        ηsl = ηsl,
        σl = σl,
        μr = μr,
        μs = μs,
        f1 = f1,
        f2 = f2,
    )
)
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = TwoStreamParameters{FT}(;
        Ω = Ω,
        ld = ld,
        α_PAR_leaf = α_PAR_leaf,
        λ_γ_PAR = λ_γ_PAR,
        λ_γ_NIR = λ_γ_NIR,
        τ_PAR_leaf = τ_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        τ_NIR_leaf = τ_NIR_leaf,
        n_layers = n_layers,
        ϵ_canopy = ϵ_canopy,
    )
)
# Set up conductance
conductance_args = (;
    parameters = MedlynConductanceParameters{FT}(;
        g1 = g1,
        Drel = Drel,
        g0 = g0,
    )
)
# Set up photosynthesis
photosynthesis_args = (;
    parameters = FarquharParameters{FT}(
        Canopy.C3();
        oi = oi,
        ϕ = ϕ,
        θj = θj,
        f = f,
        sc = sc,
        pc = pc,
        Vcmax25 = Vcmax25,
        Γstar25 = Γstar25,
        Kc25 = Kc25,
        Ko25 = Ko25,
        To = To,
        ΔHkc = ΔHkc,
        ΔHko = ΔHko,
        ΔHVcmax = ΔHVcmax,
        ΔHΓstar = ΔHΓstar,
        ΔHJmax = ΔHJmax,
        ΔHRd = ΔHRd,
    )
)
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    root_distribution = root_distribution,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
)
plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

# Canopy component args
canopy_component_args = (;
    autotrophic_respiration = autotrophic_respiration_args,
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
)
# Other info needed
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = canopy_domain)
transpiration = Canopy.PlantHydraulics.DiagnosticTranspiration{FT}()

soil_driver = PrognosticSoil(;
                             soil_α_PAR = soil_params.PAR_albedo,
                             soil_α_NIR = soil_params.NIR_albedo,
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


θs = canopy.radiation.θs.(seconds, Ref(radiation.orbital_data), ref_time)
c_co2_air = canopy.atmos.c_co2.(seconds)

(; Vcmax25, ΔHVcmax, f, ΔHRd, To, sc, pc) =
    canopy.photosynthesis.parameters
c = FT(LSMP.light_speed(earth_param_set))
planck_h = FT(LSMP.planck_constant(earth_param_set))
N_a = FT(LSMP.avogadro_constant(earth_param_set))
grav = FT(LSMP.grav(earth_param_set))
ρ_l = FT(LSMP.ρ_cloud_liq(earth_param_set))
(; ld, Ω, α_PAR_leaf, λ_γ_PAR, λ_γ_NIR) =
    canopy.radiative_transfer.parameters
energy_per_photon_PAR = planck_h * c / λ_γ_PAR
energy_per_photon_NIR = planck_h * c / λ_γ_NIR
R = FT(LSMP.gas_constant(earth_param_set))
thermo_params = canopy.parameters.earth_param_set.thermo_params
(; g1, g0, Drel) = canopy.conductance.parameters
   RT = canopy.radiative_transfer
K = @. extinction_coeff(ld, θs)
PAR = ClimaLSM.Canopy.compute_PAR.(Ref(RT), Ref(canopy.radiation), seconds)
NIR = ClimaLSM.Canopy.compute_NIR.(Ref(RT), Ref(canopy.radiation), seconds)

e_sat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        T_air,
        Ref(Thermodynamics.Liquid()),
    )
e =
    Thermodynamics.partial_pressure_vapor.(
        Ref(thermo_params),
        P_air,
        PhasePartition.(q_air),
    )
rel_hum = e ./ e_sat
DOY = @. Dates.dayofyear(ref_time + Dates.Second(floor(Int64, seconds)))
frac_diff = @. ClimaLSM.Canopy.diffuse_fraction(DOY, T_air, PAR + NIR, rel_hum, θs)
absorbance_tuple = ClimaLSM.Canopy.compute_absorbances.(
    Ref(RT),
    PAR / (energy_per_photon_PAR * N_a),
    NIR / (energy_per_photon_NIR * N_a),
    LAI,
    K,
    soil_α_PAR,
    soil_α_NIR,
    θs,
    frac_diff,
)
APAR = [absorbance_tuple[k].par[1].abs for k in 1:length(absorbance_tuple)]
ANIR = [absorbance_tuple[k].nir[1].abs  for k in 1:length(absorbance_tuple)]
TPAR = [absorbance_tuple[k].par[1].trans  for k in 1:length(absorbance_tuple)]
TNIR = [absorbance_tuple[k].nir[1].trans  for k in 1:length(absorbance_tuple)]
RPAR = [absorbance_tuple[k].par[1].refl for k in 1:length(absorbance_tuple)]
RNIR = [absorbance_tuple[k].nir[1].refl for k in 1:length(absorbance_tuple)]
medlyn_factor =
    medlyn_term.(g1, T_air, P_air, q_air, Ref(thermo_params))
β = 0.96
Rd = @. ClimaLSM.Canopy.dark_respiration(Vcmax25, β, f, ΔHRd, T_air, To, R)
An =
    ClimaLSM.Canopy.compute_photosynthesis.(
        Ref(canopy.photosynthesis),
        T_air,
        medlyn_factor,
        APAR,
        c_co2_air,
        β,
        R,
        Rd,
    )
GPP = @. ClimaLSM.Canopy.compute_GPP(An, K, LAI, Ω)
gs = @. ClimaLSM.Canopy.medlyn_conductance(g0, Drel, medlyn_factor, An, c_co2_air)
g_canopy = @. ClimaLSM.Canopy.upscale_leaf_conductance(gs, LAI, T_air, R, P_air)
r_canopy = 1 ./ g_canopy
g = LSMP.grav(earth_param_set)
R = LSMP.gas_constant(earth_param_set)
M_w = LSMP.molar_mass_water(earth_param_set)
hcm = ClimaLSM.Soil.vanGenuchten(; α = soil_vg_α, n = soil_vg_n)

ψ_soil = @. ClimaLSM.Soil.matric_potential(ClimaLSM.Soil.vanGenuchten(; α = soil_vg_α, n = soil_vg_n), (SWC - θ_r)/(0.52-θ_r))
ϵ_soil = 0.96
T_canopy_air = []
q_canopy_air = []
T_canopy = []
shf = []
lhf = []
soil_shf = []
soil_lhf = []
canopy_lhf = []
canopy_shf = []
frictionv = []
q_soil_ts =  []
LW_u = []
SW_u = energy_per_photon_NIR * N_a .* RNIR .+ energy_per_photon_PAR * N_a .* RPAR
LW_c = []
LW_soil = []
SW_c = energy_per_photon_NIR * N_a .* ANIR .+ energy_per_photon_PAR * N_a .* APAR
SW_soil = energy_per_photon_NIR * N_a .* TNIR .* (1-soil_α_NIR) .+ energy_per_photon_PAR * N_a .* TPAR .* (1-soil_α_PAR)
L_MO = []
steps = 120*48:1:180*48
for step in steps
    @info(step)
    ts_in = Thermodynamics.PhaseEquil_pTq(thermo_params, P_air[step], T_air[step], q_air[step])
    ρ_air = Thermodynamics.air_density(thermo_params, ts_in)
    cp_m = Thermodynamics.cp_m(thermo_params, ts_in)
    AI = SAI + LAI[step]
    W = exp(-AI)
    r_soil = ClimaLSM.Soil.soil_resistance(SWC[step], SWC[step], 0.0, soil_params)
    Rm_int = Thermodynamics.gas_constant_air(thermo_params, ts_in)
    ρ_soil = ρ_air *
        (T_soil[step] / T_air[step])^(Thermodynamics.cv_m(thermo_params, ts_in) / Rm_int)
     q_sat_soil =
         Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T_soil[step],
            ρ_soil,
             Thermodynamics.Liquid(),
         )
    q_soil  =  (q_sat_soil * exp(g * ψ_soil[step] * M_w / (R * T_soil[step])))
    # Ta, qa, Tc, ustar, shf, lhf
    initial_guess = [0.5 .*(T_air[step] .+ T_soil[step]), 0.5 .*(q_air[step] .+ q_soil), 0.5 .*(T_air[step] .+ T_soil[step]), u_air[step], 0.0, 0.0, 0.0]
    function flux_equality(F, x)
        #@info(x)
        #@info(F)
        TS = x[1]
        qS = x[2]
        ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_air, TS, qS)
        state_sfc =
            SurfaceFluxes.SurfaceValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
        state_in =
            SurfaceFluxes.InteriorValues(h - d_sfc, SVector{2, FT}(u_air[step], 0), ts_in)
        
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
            tol_neutral = SFP.cp_d(surface_flux_params) / 10000,
        )
        r_a = 1 / (conditions.Ch * SurfaceFluxes.windspeed(sc))
        ustar = conditions.ustar
        
        x[4] = ustar
        x[5] = conditions.shf
        x[6] = conditions.lhf
        x[7] = conditions.L_MO
        Cs_bare =@.  0.4/0.13*(0.01*ustar/(1.5e-5))^(-0.45)
        Cs = @. W*Cs_bare + (1-W)*Cs_canopy
        r_b = @. 1/Cv/sqrt(ustar)/Ild
        r_ah_cs = @. 1/Cs/ustar
        
        soil_shf = -ρ_air * cp_m * (TS - T_soil[step]) / r_ah_cs
        soil_evap = -ρ_air * (qS - q_soil) / (r_ah_cs + r_soil)
        
        T_c = x[3]
        q_canopy =  
            Thermodynamics.q_vap_saturation_generic(
                thermo_params,
                T_c,
                ρ_air,
                Thermodynamics.Liquid(),
            )
        canopy_shf  = -ρ_air * cp_m * (TS - T_c) / (r_b/AI)
        canopy_evap = -ρ_air * (qS - q_canopy) / (r_canopy[step] + r_b/LAI[step])
        LW_d_canopy = (1 - ϵ_canopy) * LW_d[step] + ϵ_canopy * _σ * T_c^4
        LW_u_soil = ϵ_soil * _σ * T_soil[step]^4 + (1 - ϵ_soil) * LW_d_canopy
        canopy_LW = ϵ_canopy * LW_d[step] - 2 * ϵ_canopy * _σ * T_c^4 + ϵ_canopy * LW_u_soil
        
        F[1] = canopy_shf + soil_shf - conditions.shf
        F[2] = (canopy_evap + soil_evap - conditions.evaporation) * _LH_v0
        F[3] = (-canopy_LW -SW_c[step] + canopy_evap * _LH_v0 + canopy_shf)
        F[4] = 0.0
        F[5] = 0.0
        F[6] = 0.0
        F[7] = 0.0
    end
    soln = nlsolve(
        flux_equality,
        initial_guess,
        ftol = 1e-8,
        iterations = 100
    )
    T_c = soln.zero[3]
    ust= soln.zero[4]
    T_ca = soln.zero[1]
    q_ca = soln.zero[2]
    Cs_bare =@.  0.4/0.13*(0.01*ust/(1.5e-5))^(-0.45)
    Cs = @. W*Cs_bare + (1-W)*Cs_canopy
    r_b = @. 1/Cv/sqrt(ust)/Ild
    r_ah_cs = @. 1/Cs/ust
     q_canopy =  
         Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T_c,
             ρ_air,
             Thermodynamics.Liquid(),
         )
    push!(L_MO, soln.zero[7])
    push!(q_soil_ts, q_soil)
    push!(T_canopy_air, T_ca)
    push!(q_canopy_air, q_ca)
    push!(frictionv, ust)
    push!(shf, soln.zero[5])
    push!(lhf, soln.zero[6])
    push!(soil_shf, -ρ_air * cp_m * (T_ca - T_soil[step]) / r_ah_cs)
    push!(soil_lhf, -ρ_air * _LH_v0 * (q_ca - q_soil) / (r_soil+r_ah_cs))
    push!(canopy_lhf,  -ρ_air * _LH_v0 * (q_ca - q_canopy) / (r_canopy[step]+r_b/LAI[step]))
    push!(canopy_shf,  -ρ_air * cp_m * (T_ca - T_c) / (r_b/AI))
    LW_d_canopy = (1 - ϵ_canopy) * LW_d[step] + ϵ_canopy * _σ * T_c^4
    LW_u_soil = ϵ_soil * _σ * T_soil[step]^4 + (1 - ϵ_soil) * LW_d_canopy
    canopy_LW = ϵ_canopy * LW_d[step] - 2 * ϵ_canopy * _σ * T_c^4 + ϵ_canopy * LW_u_soil
    push!(LW_c, canopy_LW)
    push!(LW_u, (1 - ϵ_canopy) * LW_u_soil + ϵ_canopy * _σ * T_c^4)
    push!(T_canopy, T_c)
    push!(LW_soil, ϵ_soil * LW_d_canopy - ϵ_soil * _σ * T_soil[step]^4)
end
    
idx_end = length(T_canopy_air)
days = seconds[steps[1:idx_end]] ./ 24 ./ 3600
Plots.plot(days, T_canopy_air, label = "T_ca")
Plots.plot!(days, T_air[steps[1:idx_end]], label = "Tair")
Plots.plot!(days, T_soil[steps[1:idx_end]], label = "Tsoil")
Plots.plot!(days, T_canopy, label = "Tcanopy")

Plots.plot(days, LW_u, label = "Model", title = "LW out")
Plots.plot!(days, LW_OUT[steps[1:idx_end]], label = "Data")

Plots.plot(days, SW_u[steps[1:idx_end]], label = "Model", title = "SW out")
Plots.plot!(days, SW_OUT[steps[1:idx_end]], label = "Data")


Plots.plot(days, q_canopy_air, label = "q_ca")
Plots.plot!(days, q_air[steps[1:idx_end]], label = "qair")
Plots.plot!(days, q_soil_ts, label = "qsoil")

# G still looks horrible
Plots.plot(days, -1 .* LW_soil .- SW_soil[steps[1:idx_end]] .+ soil_shf .+ soil_lhf, label = "Model")
Plots.plot!(days, G[steps[1:idx_end]], label = "Data")
# canopy SHF looks terrible
Plots.plot(days, soil_shf, label  = "soil shf")
Plots.plot!(days, canopy_shf, label  = "canopy shf")
Plots.plot!(days, H[steps[1:idx_end]], label = "data", xlim = [120,130])

Plots.plot(days, soil_lhf, label = "soil lhf")
Plots.plot!(days, canopy_lhf, label = "canopy lhf")
Plots.plot!(days, LE[steps[1:idx_end]], label = "data")
Plots.plot(days, frictionv)
