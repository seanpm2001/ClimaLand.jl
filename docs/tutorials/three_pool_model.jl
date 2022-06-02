using OrdinaryDiffEq: ODEProblem, solve, RK4
using Plots
using Test
using ClimaCore
using LinearAlgebra
using Plots
using DifferentialEquations
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains

# Import the functions we are extending for our model:
import ClimaLSM: name, make_rhs, prognostic_vars
import ClimaLSM.Domains: coordinates

abstract type AbstractCarbonModel{FT} <: AbstractModel{FT} end

struct BulkThreePools{FT} <: AbstractCarbonModel{FT}
    live_turnover_rate::FT
    dead_turnover_rate::Function
    GPP::FT
    "Growth Respiration Term - function of NSC-> live biomass flux"
    R_g::FT
    "Maintenance - depends on climate"
    m_R::FT
    T_soil::FT
end

Base.@kwdef struct Q10Respiration{FT<:AbstractFloat}
    Q10::FT = 2.0
    base_turnover_time::FT = 40.0
    T_ref::FT = 25.0
end

function turnover_rate(T::FT, resp_model::Q10Respiration) where {FT}
    (; Q10, base_turnover_time, T_ref) = resp_model
    return Q10 ^((T-T_ref)/FT(10.0))/base_turnover_time
end


ClimaLSM.name(model::BulkThreePools) = :carbon

ClimaLSM.prognostic_vars(::BulkThreePools) = (:NSC, :live, :dead)

ClimaLSM.Domains.coordinates(model::BulkThreePools{FT}) where {FT} =
    FT.([0.0]);

function ClimaLSM.make_rhs(model::BulkThreePools{FT}) where {FT}
    function rhs!(dY, Y, p, t) # gets the entire Y
        # NPP = GPP - R_m - R_g
        # G = function (x1,x2, phen, environment)

        NPP = @.(model.GPP - model.m_R*Y.carbon.live - model.R_g)
        G = NSC_transfer.(model.GPP)
        @. dY.carbon.NSC = NPP - G
        @. dY.carbon.live = G - Y.carbon.live*live_turnover_rate
        @. dY.carbon.dead = Y.carbon.live*live_turnover_rate - Y.carbon.dead*model.turnover_rate(model.T_soil)
    end
    return rhs!
end

# stand in for future function
function NSC_transfer(GPP::FT) where{FT}
    return FT(0.1)*GPP
end

live_turnover_rate = 1.0/10.0 # years
dead_resp_model = Q10Respiration{Float64}(;Q10 = 2.0, base_turnover_time = 50.0)
dead_turnover_rate = (T) -> turnover_rate(T, dead_resp_model)
GPP = 60 # PgC/year
m_R = 0.2 # /year
R_g = 0.0 # PgC/year
T_soil = 25.0 # C
carbon = BulkThreePools{Float64}(
    live_turnover_rate,
    dead_turnover_rate,
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
@test sum(abs.(dY.carbon.pool .- [0.0, 0.0])) < 1e-14


t0 = 0.0;
tf = 1000.0; # years
dt = 1.0; # years
prob = ODEProblem(ode_function!, Y, (t0, tf), p);
sol = solve(prob, RK4(); dt = dt, reltol = 1e-6, abstol = 1e-6);

expected = [600.0, 1500.0]

pool1 = [sol.u[k].carbon.pool[1] for k in 1:1:length(sol.t)] # time series
pool2 = [sol.u[k].carbon.pool[2] for k in 1:1:length(sol.t)]
plot(sol.t, pool1)
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
