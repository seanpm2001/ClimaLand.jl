using Test
using Statistics
using DifferentialEquations
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, SSPRK33, RK4
using ClimaCore

if !("." in LOAD_PATH) # for ease of include
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots


FT = Float64

const ν = FT(0.495);
const Ksat = FT(0.0443 / 3600 / 100); # m/s
const S_s = FT(1e-3); #inverse meters
const vg_n = FT(2.0);
const vg_α = FT(2.6); # inverse meters
const vg_m = FT(1) - FT(1) / vg_n;
const θ_r = FT(0);
const zmax = FT(0);
const zmin = FT(-10);
const nelems = 50;

soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelems);
top_flux_bc = FT(-1e-5);
bot_flux_bc = FT(0.0);
boundary_fluxes = (top_flux_bc = top_flux_bc, bot_flux_bc = bot_flux_bc)
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r);

soil = Soil.RichardsModel{FT}(;
    param_set = params,
    domain = soil_domain,
    boundary_exchanges = boundary_fluxes,
)

Y, p, coords = initialize(soil)

# specify ICs
function init_soil!(Ysoil, z, params)
    function hydrostatic_profile(
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack ν, vg_α, vg_n, vg_m, θ_r = params
        #unsaturated zone only, assumes water table starts at z_∇
        z_∇ = FT(-10)# matches zmin
        S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
        ϑ_l = S * (ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
end

init_soil!(Y, coords, soil.param_set)

soil_ode! = make_ode_function(soil)

t0 = FT(0);
tf = FT(3600);
dt = FT(40);
# Try saving `p` differently - copy of integrator is required.
# Here we save at every timestep, and turn off adaptive stepping.
sv1 = SavedValues(FT, ClimaCore.Fields.FieldVector)
sv2 = SavedValues(FT, Vector{FT})                           
cb = SavingCallback((t, u, integrator) -> copy(integrator.p), sv1)
cb2 = SavingCallback((t, u, integrator) -> ([parent(integrator.p.soil.ψ)[end]]), sv2)
prob = ODEProblem(soil_ode!, Y, (t0, tf), p);
sol = solve(prob, RK4(); dt = dt, callback = CallbackSet(cb,cb2), save_every_step = true, adaptive = true);
N = length(sol.t)
ϑ = [parent(sol.u[k].soil.ϑ_l) for k in 1:1:N]
actual_S = [effective_saturation.(ν, ϑ[k], θ_r) for k in 1:1:N]
actual_K = [hydraulic_conductivity.(Ksat, vg_m, actual_S[k]) for k in 1:1:N]
actual_ψ = [pressure_head.(vg_α, vg_n, vg_m,θ_r, ϑ[k],ν, S_s) for k in 1:1:N]
saved_ψ = [parent(sv1.saveval[k].soil.ψ) for k in 1:1:N]
saved_K = [parent(sv1.saveval[k].soil.K) for k in 1:1:N]
# Saved values agree exactly with values computed offline from Y(t).
# This is confusing - when are the two saved? `p` is updated in the RHS to the current step,
# `Y` then updated to the next step, so there is a point when Y and p are at different steps.
Δψ = sum([mean(abs.(actual_ψ[k] .- saved_ψ[k])) for k in 2:1:N])
ΔK = sum([mean(abs.(actual_K[k] .- saved_K[k])) for k in 2:1:N])


# two saving methods agree
sum([parent(sv1.saveval[k].soil.ψ)[end] - sv2.saveval[k][1] for k in 1:1:N])


# note that saved and computed values do not agree if we use adaptive stepping, or interpolations,
# because that is applied to Y to get the value at t but not to p. But p should still be p(t).
# so ideally differences are small?


# What if we save at a non-interval of the timestep?

sv4 = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb4 = SavingCallback((t, u, integrator) -> copy(integrator.p), sv4,saveat = 0.0:250:tf)
prob = ODEProblem(soil_ode!, Y, (t0, tf), p);
sol4 = solve(prob, RK4(); dt = dt, callback = cb4, adaptive = false,saveat = 0.0:250:tf)
N =  length(sv4.t)
saved_ψ4 = [parent(sv4.saveval[k].soil.ψ) for k in 1:1:N]
saved_K4 = [parent(sv4.saveval[k].soil.K) for k in 1:1:N]


ϑ4 = [parent(sol4.u[k].soil.ϑ_l) for k in 1:1:N]
actual_S4 = [effective_saturation.(ν, ϑ[k], θ_r) for k in 1:1:N]
actual_K4 = [hydraulic_conductivity.(Ksat, vg_m, actual_S[k]) for k in 1:1:N]
actual_ψ4 = [pressure_head.(vg_α, vg_n, vg_m,θ_r, ϑ[k],ν, S_s) for k in 1:1:N]

Δψsa = sum([mean(abs.(saved_ψ4[k] .- actual_ψ4[k])) for k in 2:1:N])
ΔKsa = sum([mean(abs.(saved_K4[k] .- actual_K4[k])) for k in 2:1:N])
# now these dont agree, possibly because Y is interpolated to get to the save time differently
# than p?

# What if we use an interval of the timestep?
sv5 = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb5 = SavingCallback((t, u, integrator) -> copy(integrator.p), sv5,saveat = 0.0:320:tf)
prob = ODEProblem(soil_ode!, Y, (t0, tf), p);
sol5 = solve(prob, RK4(); dt = dt, callback = cb5, adaptive = false,saveat = 0.0:320:tf)
N =  length(sv5.t)
saved_ψ5 = [parent(sv5.saveval[k].soil.ψ) for k in 1:1:N]
saved_K5 = [parent(sv5.saveval[k].soil.K) for k in 1:1:N]


ϑ5 = [parent(sol5.u[k].soil.ϑ_l) for k in 1:1:N]
actual_S5 = [effective_saturation.(ν, ϑ[k], θ_r) for k in 1:1:N]
actual_K5 = [hydraulic_conductivity.(Ksat, vg_m, actual_S[k]) for k in 1:1:N]
actual_ψ5 = [pressure_head.(vg_α, vg_n, vg_m,θ_r, ϑ[k],ν, S_s) for k in 1:1:N]

Δψsa = sum([mean(abs.(saved_ψ5[k] .- actual_ψ5[k])) for k in 2:1:N])
ΔKsa = sum([mean(abs.(saved_K5[k] .- actual_K5[k])) for k in 2:1:N])
# no, do not agree.
