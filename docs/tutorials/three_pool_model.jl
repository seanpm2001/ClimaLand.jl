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
ClimaLSM.domain(::AbstractCarbonModel) = :surface
struct BulkThreePools{FT,D} <: AbstractCarbonModel{FT}
    live_turnover_time::FT
    dead_turnover_time::FT
    GPP::FT
    R_g::FT
    m_R::FT
    domain::D
end
function BulkThreePools{FT}(; live_turnover_time::FT, dead_turnover_time::FT, GPP::FT, R_g::FT, m_R::FT,
                            domain::ClimaLSM.Domains.AbstractDomain{FT} = ClimaLSM.Domains.Point(z_sfc=0.0),
                            ) where {FT}
    args = (live_turnover_time, dead_turnover_time, GPP, R_g, m_R)
    return BulkThreePools{FT, typeof(domain)}(args..., domain)
end

ClimaLSM.name(model::BulkThreePools) = :carbon

ClimaLSM.prognostic_vars(::BulkThreePools) = (:NSC, :live, :dead)
ClimaLSM.prognostic_types(::BulkThreePools{FT}) where{FT}= (FT, FT, FT)


function ClimaLSM.make_rhs(model::BulkThreePools{FT}) where {FT}
    function rhs!(dY, Y, p, t) 
        NPP = @.(model.GPP - model.m_R*Y.carbon.live - model.R_g)
        G = NSC_transfer.(model.GPP)
        @. dY.carbon.NSC = NPP - G
        @. dY.carbon.live = G - Y.carbon.live/live_turnover_time
        @. dY.carbon.dead = Y.carbon.live/live_turnover_time - Y.carbon.dead/dead_turnover_time
    end
    return rhs!
end

# stand in for future function
function NSC_transfer(GPP::FT) where{FT}
    return FT(0.1)*GPP
end


live_turnover_time = 10.0 # years
dead_turnover_time = 100.0 # years
GPP = 60.0 # PgC/year
m_R = 0.2 # /year
R_g = 0.0 # PgC/year

carbon = BulkThreePools{Float64}(; live_turnover_time = live_turnover_time,
                                 dead_turnover_time = dead_turnover_time,
                                 GPP = GPP,
                                 R_g = R_g,
                                 m_R = m_R)

Y, p, coords = initialize(carbon)

Y.carbon.live= 600.0
Y.carbon.dead = 1500.0
ode_function! = make_ode_function(carbon);
dY = similar(Y)
ode_function!(dY, Y, p, 0.0)
t0 = 0.0;
tf = 200.0;
dt = 1.0;
prob = ODEProblem(ode_function!, Y, (t0, tf), p);
sol = solve(prob, RK4(); dt = dt, reltol = 1e-6, abstol = 1e-6);

nsc_pool = [parent(sol.u[k].carbon.NSC)[1] for k in 1:1:length(sol.t)] # time series
live_pool = [parent(sol.u[k].carbon.live)[1] for k in 1:1:length(sol.t)] # time series
dead_pool = [parent(sol.u[k].carbon.dead)[1] for k in 1:1:length(sol.t)]
plot(sol.t, nsc_pool, label= "NSC")
plot!(sol.t, live_pool, label = "Live carbon")
plot!(sol.t,dead_pool, label = "Dead carbon")
