using OrdinaryDiffEq: ODEProblem, solve, RK4
using Plots
using Test
using ClimaCore
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains

# Import the functions we are extending for our model:
import ClimaLSM: name, make_rhs, prognostic_vars
import ClimaLSM.Domains: coordinates
include("carbon_model_utils.jl")


T_soil = (t) -> eltype(t)(25.0) # C
live_resp_model = Q10Respiration{Float64}(;Q10 = 2.0, base_turnover_time = 10.0) 
dead_resp_model = Q10Respiration{Float64}(;Q10 = 2.0, base_turnover_time = 100.0)
GPP = 60.0 # PgC/year
m_R = 0.2 # /year
R_g = 0.0 # PgC/year
carbon = BulkThreePools{Float64, typeof(live_resp_model), typeof(dead_resp_model)}(
    live_resp_model,
    dead_resp_model,
    GPP,
    R_g,
    m_R,
    T_soil
) # sets all the attributes in carbon model
Y, p, coords = initialize(carbon)

Y.carbon.NSC = [60.0]
Y.carbon.live = [600.0]
Y.carbon.dead = [1500.0]
ode_function! = make_ode_function(carbon);
#dY/dt = ode_function(dY,Y,p,t) # updates dY in place
dY = similar(Y)
ode_function!(dY, Y, p, 0.0)
#@test sum(abs.(dY.carbon.pool .- [0.0, 0.0])) < 1e-14


t0 = 0.0;
tf = 1000.#0; # years
dt = 1.0; # years
prob = ODEProblem(ode_function!, Y, (t0, tf), p);
sol = solve(prob, RK4(); dt = dt, reltol = 1e-6, abstol = 1e-6);

expected = [600.0, 1500.0]

NSC = [sol.u[k].carbon.NSC[1] for k in 1:1:length(sol.t)]
live = [sol.u[k].carbon.live[1] for k in 1:1:length(sol.t)] # time series
dead = [sol.u[k].carbon.dead[1] for k in 1:1:length(sol.t)]
plot(sol.t, live)
plot!(sol.t, dead)
plot!(sol.t,NSC)
# Follow up
# [X] Test steady state solution, look at time evolution plot
# time-varying NPP
# turnover times not constant, but parameterized
# time-varying soil moisture and temp 
# rhs in matrix form
# Hooking up with soil?

#=
function NPP(model::BulkTwoPools{FT},npp_model::Standalone)
    return model.NPP
end

function NPP(model::BulkTwoPools{FT}, npp_model::Coupled)
    p.NPP
end
=#
