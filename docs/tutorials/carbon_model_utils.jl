abstract type AbstractCarbonModel{FT} <: AbstractModel{FT} end
abstract type AbstractRespirationModel{FT <: AbstractFloat} end
struct BulkThreePools{FT, LRM, DRM} <: AbstractCarbonModel{FT}
    live_carbon_respiration::LRM
    dead_carbon_respiration::DRM
    GPP::FT
    "Growth Respiration Term - function of NSC-> live biomass flux"
    R_g::FT
    "Maintenance - depends on climate"
    m_R::FT
    T_soil::Function
end

Base.@kwdef struct Q10Respiration{FT} <: AbstractRespirationModel{FT}
    Q10::FT = 2.0
    base_turnover_time::FT = 40.0
    T_ref::FT = 25.0
end

function turnover_rate(T::FT, resp_model::Q10Respiration{FT}) where {FT}
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
        G = NSC_transfer.(model, Y.carbon.NSC)
        if Y.carbon.NSC < G > dY.carbon.NSC:
            G = dY.carbon.NSC
        @. dY.carbon.NSC = NPP - G
        live_turnover_rate = turnover_rate(model.T_soil(t), model.live_carbon_respiration)
        dY.carbon.live = G .- Y.carbon.live .*live_turnover_rate
        dY.carbon.dead = Y.carbon.live .*live_turnover_rate .- Y.carbon.dead .*turnover_rate(model.T_soil(t), model.dead_carbon_respiration)
    end
    return rhs!
end

# stand in for future function
function NSC_transfer(model::BulkThreePools{FT}, NSC) where{FT}
    return FT(0.8)*GPP
end

#=

dY/dt = f(Y,parameters, environment, t) + g(Y,parameters, environment,...) = rhs1! + rhs2!
dY = Y(t+dt) - Y(t) = ∫ f(Y,....) dt = f(Y...)*dt + (dt^2)

if NSC - (NPP-G)*dt <0
    NSC = 0.0 
    remainder_flux = G - NSC/dt
else
    remainder_flux = G
end
dNSC/dt = input - output
dNSC/dt = NPP - G
dY/dt = -Y/τ -> Y(t) exponetially decays
dNSC_dt = NPP - live*parameter
dNSC/dt = min(dNSC_dt,)

dNSC/dt = NPP - min(NSC, live)*parameter 
dNSC/dt = NPP - min(NSC/turnover_rate, live*parameter) -> NSC/turnover_rate = gC / days
G = f(supply from NSC, demand by live biomass)

k(NSC, NSCcritical) = min( 1,  (1 - e^-alpha*NSC)/(1 - e^-alpha*NSCcritical) )
if NSC < NSCcritical   
    G = +NSC/τ
else
    G = live*parameter
    G = f(model, optimality, fixed allocation, source/sink)
end
=#