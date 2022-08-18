##### Constructors for required terms ######
abstract type AlbedoModel{FT} end
abstract type EmissivityModel{FT} end
abstract type TauModel{FT} end
abstract type BucketModel{FT} end

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

function makeConstantSnowModel(init::SnowState{FT}, a, e, t, b, params::ModelParams{FT}, dt)::SnowModel{FT}
    a_model = ConstantAlbedo{FT}(a)
    e_model = ConstantEmissivity{FT}(e)
    t_model = ConstantTau{FT}(t)
    b_model = ConstantBucket{FT}(b)
    init_derived = DerivedParams{FT}(zeros(FT, length(fieldnames(DerivedParams)))...)
    init_flux = FluxTerms{FT}(zeros(FT, length(fieldnames(FluxTerms)))...)
    return SnowModel{FT}(a_model, e_model, t_model, b_model, init, params, init_derived, init_flux, dt)
end

##### Functions for derived/tracer terms ######
function ρ_s(s::SnowState{FT}, p::ModelParams{FT})::FT
    return s.swe / s.z * p.ρ_l
end
#parametric type 
function l(e_int::FT, Lfus::FT)::FT where{FT}
    val = e_int /Lfus
    return max(min(val, one(FT)), zero(FT))
end

function cp_s(p::ModelParams{FT}, d::DerivedParams{FT})::FT
    l = d.l
    return (one(FT)-l)*p.cp_i+l*p.cp_w
end

function T̄(s::SnowState{FT}, d::DerivedParams{FT}, p::ModelParams{FT})::FT
    return (s.e_int-d.l*p.Lfus) / d.cp_s + p.T_F
end

function Tsurf(d::DerivedParams{FT}, f::Forcing{FT})::FT
    #assumes linear profile and T_bottom = T_soil
    return 2*d.T̄ - f.T_soil
end

function ρ_a(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
    #can be computed by CliMA Thermodynamics
end

function g_aw(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
end

function g_as(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
end

function qw(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
    #computed by CliMA Thermodynamics
end

function qi(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
    #computed by CliMA Thermodynamics
end

function cp_a(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
    #computed by CliMA Thermodynamics
end

function κ_s(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
    # make this into another model struct
    # or make a constant
end

function κ_a(m::SnowModel{FT}, f::Forcing{FT})::FT
    #something
end

function α(albedo_model::ConstantAlbedo{FT})::FT
    return albedo_model.α
end

#function α(albedo_model::DynamicAlbedo,)
#    #calculate something
#    a = albedo_model.p1 + albedo_model.p2
#    return a
#end

function ϵ(emis_model::ConstantEmissivity{FT})::FT
    return emis_model.ϵ
end

function τ(runoff_model::ConstantTau{FT})::FT
    return runoff_model.τ
end

function θ_max(bucket_model::ConstantBucket{FT})::FT
    return bucket_model.θ_max
end

##### Functions for flux terms ######
function El(d::DerivedParams{FT}, f::Forcing{FT})::FT
    return d.ρ_a * d.g_aw * (d.qw - f.q_air)
end

function Ei(d::DerivedParams{FT}, f::Forcing{FT})::FT
    return d.ρ_a * d.g_as * (d.qi - f.q_air)
end

function Ẽ(p::ModelParams{FT}, d::DerivedParams{FT}, f::Forcing{FT})::FT
    l = d.l
    el = El(d, f)
    ei = Ei(d, f)
    return (l*el + (one(FT)-l)*ei) / p.ρ_l
end

function G(p::ModelParams{FT}, d::DerivedParams{FT}, f::Forcing{FT})::FT
    return (one(FT)-d.α)*f.Gsol + d.ϵ*p.σ*(f.T_air^4 - d.Tsurf^4)
end

function Qconv(d::DerivedParams{FT}, f::Forcing{FT})::FT
    return d.ρ_a * d.cp_a * d.g_as * (f.T_air - d.Tsurf)
end

function Qlat(p::ModelParams{FT}, d::DerivedParams{FT}, f::Forcing{FT})::FT
    #confirm these two formulas!
    elg = p.Lvap
    esg = p.cp_i*(d.Tsurf - p.T_F) + p.Lfus + p.Lvap
    l = d.l
    ei = Ei(d, f)
    el = El(d, f)
    return elg*l*el + esg*(one(FT)-l)*ei
end

function imax(s::SnowState{FT}, d::DerivedParams{FT}, p::ModelParams{FT})::FT
    return s.z / s.swe * d.θ_max * p.Lfus
end

function QR(s::SnowState{FT}, d::DerivedParams{FT})::FT
    Δi = s.i - d.imax
    return  (Δi > zero(FT)) ? (Δi * d.ρ_s * s.z / d.τ) : zero(FT)
end

function R̃(p::ModelParams{FT}, d::DerivedParams{FT}, f::FluxTerms{FT})::FT
    return  f.QR / (p.ρ_l * p.Lfus)
end

function dIdt(f::FluxTerms{FT})::FT
    return f.G + f.Qcond - f.Qlat - f.QR
end

function dSWEdt(f::FluxTerms{FT})::FT
    return f.P̃ - f.Ẽ - f.R̃
end

function dzdt(f::Forcing{FT})
    #something
end

function zeromodel!(m::SnowModel{FT})
    m.state.e_int = 0.0
    m.state.swe = 0.0
    m.state.z = 0.0
    m.state.has_snow = false
    m.derived.ρ_s = 0.0
    m.derived.imax = Inf #?
end

function updateParams!(m::SnowModel{FT}, f::Forcing{FT})
    if !m.state.has_snow & f.P == zero(FT)
        return
    elseif m.state.e_int > m.p.Lfus
        zeromodel!(model)
        return
    end
    #is there a way to handle this without boolean to make it faster? (short-circuit or boolean-multiply the term-setting?)
    #handle in case we start with z = 0 or SWE = 0 and/or I = 0?
    #if internal energy is high enough to melt entire snowpack, do something else, otherwise:
    #or a case where z updates to a point where z < SWE, etc. think throguh all cases

    #order is purposeful
    m.derived.ρ_s = ρ_s(m.state,m.params)
    m.derived.l = l(m.state, m.params)
    m.derived.cp_s = cp_s(m.params, m.derived)
    m.derived.T̄ = T̄(m.state, m.derived, m.params)
    m.derived.Tsurf = Tsurf(model.derived, f) #move this one?

    m.derived.α = α(m.α_model)
    m.derived.ϵ = ϵ(m.ϵ_model)
    m.derived.τ = τ(m.τ_model)
    m.derived.θ = θ_max(m.θ_model)
    m.derived.imax = imax(m.state, m.derived, m.params)

    m.derived.ρ_a = ρ_a(m, f) #edit args if need be
    m.derived.cp_a = cp_a(m, f)
    m.derived.g_as = g_as(m, f)
    m.derived.g_aw = g_aw(m, f)
    m.derived.qi = qi(m, f)
    m.derived.qw = qf(m, f)
    m.derived.κ_a = κ_a(m, f)
    m.derived.κ_s = κ_s(m, f)

    m.flux.G = G(m.params, m.derived, f)
    m.flux.Qconv = Qconv(m.derived, f)
    m.flux.Qlat = Qlat(m.params, m.derived, f)
    m.flux.QR = QR(m.state, m.derived)
    
    m.flux.P̃ = f.P
    m.flux.Ẽ = Ẽ(m.params, m.derived, f)
    m.flux.R̃ = R̃(m.params, m.derived, m.flux)
end

function update!(model::SnowModel{FT}, f::Forcing{FT})
    #first update the derived parameters list
    updateParams!(model, f)
    δI = dIdt(model.flux)
    δSWE = dSWEdt(model.flux)
    δz = dzdt(f)
    #should i be using something that isn't Euler method?
    model.state.e_int += δI / (m.state.z * m.derived.ρ_s) * model.dt
    model.state.swe += δSWE * model.dt
    model.state.z += δz * model.dt
end

######################################

#for testing
alb_const = 0.4
emis_const = 0.99
tau_const = 2.0
bucket_const = 12.0
boltz = 5.67e-8
lvap = 1.0
lfus = 1.0
cpw = 1.0
cpi = 1.0
tf = 0.0
pl = 0.0
dt = convert(FT, 1.0)

alb_model = ConstantAlbedo{FT}(alb_const)
emis_model = ConstantEmissivity{FT}(emis_const)
tau_model = ConstantTau{FT}(tau_const)
bucket_model = ConstantBucket{FT}(bucket_const)
initial_state = SnowState{FT}(0.0, 0.0, 0.0, false)
mparams = ModelParams{FT}(boltz, lvap, lfus, cpw, cpi, tf, pl)
snow_model = makeConstantSnowModel(initial_state, alb_const, emis_const, tau_const, bucket_const, mparams, dt)

#import the CliMA stuff and use it (needs total humidity and pressure?)
#include your neural network into this model
#handle zeros/Inf case for swe = 0, z = 0, I = 0, I = enough to melt snowpack
#shift code across multiple files
#see if you can get a run to complete