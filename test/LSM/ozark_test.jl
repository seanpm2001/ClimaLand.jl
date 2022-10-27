using Test
using DiffEqCallbacks
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, RK4
using ClimaCore
import CLIMAParameters as CP
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using DelimitedFiles
using Dierckx
using Plots
using Statistics
using Dates
using Plots

using ClimaLSM
using ClimaLSM.Domains: Column, PlantHydraulicsDomain
using ClimaLSM.Soil
using ClimaLSM.PlantHydraulics
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))


FT = Float64

precip_θ_T =
    readdlm("test/LSM/precip.csv", ',', skipstart = 1)
et = readdlm("test/LSM/ET.csv", ',', skipstart = 1)
#et_spline = Spline1D(et[:, 1], max.(0.0, et[:, 2]))
p_spline = Spline1D(precip_θ_T[:, 1], -precip_θ_T[:, 2])
t = precip_θ_T[:, 1]

## Natan's data
data = readdlm("test/LSM/holtzman_clima_output_april1.csv", ',',skipstart = 1)
natan_et = data[2:end, 15] * 18 / 1e3 / 1e3 # convert to m^3/m^2/s
swc_column = data[2:end, 20]
swc_surface = data[2:end, 19]
swc_obs_surface = data[2:end, 12]
lwp = data[2:end,18]
dates = data[2:end, 1]
data = nothing # frees memory
dates_julia = tryparse.(DateTime, dates)
our_year = dates_julia[Dates.year.(dates_julia) .== 2005]
seconds = Dates.value.(our_year .- our_year[1]) ./ 1000
our_year_swc_column = FT.(swc_column[Dates.year.(dates_julia) .== 2005])
our_year_swc_surface = FT.(swc_surface[Dates.year.(dates_julia) .== 2005])
our_year_swc_obs_surface = FT.(swc_obs_surface[Dates.year.(dates_julia) .== 2005])
our_year_lwp = FT.(lwp[Dates.year.(dates_julia) .== 2005])
our_year_et = FT.(natan_et[Dates.year.(dates_julia) .== 2005])
et_spline = Spline1D(seconds, our_year_et)

# Capping T >0 and P <0. Not sure. P>0 seemed to be fine, but ET <0 led to some issues with
# the integration (numerical instability)
precip_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0
transpiration_function(t::FT) where {FT} =
    et_spline(t) > 0.0 ? et_spline(t) : 0.0
    #transpiration =
    #PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))

#@testset " Soil plant hydrology LSM integration test" begin
    # Plant hydraulics params
    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    earth_param_set = create_lsm_parameters(FT)
    SAI = FT(0.00242) # Basal area per ground area
    LAI = FT(4.2) # from Yujie's paper; checkout out yujiens source
    f_root_to_shoot = FT(1.0 / 5.0) # guess
    RAI = SAI * f_root_to_shoot # following CLM
    area_index = (root = RAI, stem = SAI, leaf = LAI)
    plant_K_sat = (root = FT(1.7658E-12), stem = FT(1.7658E-12), leaf = FT(1.7658E-12)) #m^3/m^2/s/Pa, from Natan (10 mol/s/m^2/MPa) 
    plant_vg_α = FT(0.06) # perhaps changing to linear in future
    plant_vg_n = FT(1.85) # vg curve very non-linear -- precision when close to saturation
    plant_vg_m = FT(1) - FT(1) / plant_vg_n
    plant_ν = FT(0.495)
    plant_S_s = FT(1e-3)
    # currently hardcoded to match the soil coordinates. this has to
    # be fixed eventually.
    z_root_depths = -Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0 # the root tips are at middle of soil elements, made calculating source term easier by avoiding to need to make extrapolation because source depends on soil pressure and vegetation water pressure
    # in previous Ozark test: z_root_depths = reverse(-Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0) # check order of soil z
    n_stem = Int64(1) 
    n_leaf = Int64(1)
    # need to adapt to make compartment sizes variable - check whether necessary actually
    Δz = FT(9.25) # height of compartments,  # height of trunk = 18.5, from Yujie's paper 
    plant_hydraulics_domain = 
        PlantHydraulicsDomain(z_root_depths, n_stem, n_leaf, Δz) # maybe pass in heights to make variable instead; if using clima core we wont need
    # 0.5 is from Natan
    function root_distribution(z::T) where {T}
        return T(1.0 / 0.5) * exp(z / T(0.5))
    end 
        
    plant_hydraulics_ps =
        PlantHydraulics.PlantHydraulicsParameters{FT, typeof(earth_param_set)}(
            area_index,
            plant_K_sat,
            plant_vg_α,
            plant_vg_n,
            plant_vg_m,
            plant_ν,
            plant_S_s,
            root_distribution,
            earth_param_set,
        )
    plant_hydraulics_args = (
        domain = plant_hydraulics_domain, 
        parameters = plant_hydraulics_ps, 
        transpiration = PrescribedTranspiration{FT}(transpiration_function)
    )

    # Soil params
    zmin = FT(-2.0)
    zmax = FT(0.0)
    nelements = 10
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelements)
    soil_ν = FT(0.55)
    soil_K_sat = FT(4e-7) # matches Natan, m/s
    soil_S_s = FT(1e-3) # inverse meters, guess
    soil_vg_n = FT(1.5)
    soil_vg_α = FT(0.10682); # inverse meters. From Natan (10.9/MPa)
    soil_vg_m = FT(1) - FT(1) / soil_vg_n
    θ_r = FT(0)
    soil_ps = Soil.RichardsParameters{FT}(soil_ν, soil_vg_α, soil_vg_n, soil_vg_m, soil_K_sat, soil_S_s, θ_r)
    #top_flux_bc = precip_function(t)
    #bot_flux_bc = FT(0.0)
    #boundary_fluxes = FluxBC{FT}(top_flux_bc, bot_flux_bc)
    #soil_args = (domain = soil_domain, parameters = soil_ps, boundary_conditions = boundary_fluxes,)   
    sources = ()
    soil_args = (
        parameters = soil_ps,
        domain = soil_domain,
        boundary_conditions = PrecipFreeDrainage{FT}(precip_function),
        sources = sources,
    )

    # Integrated plants and soil model
    # previously land_args =
    # (precipitation = precip_function, transpiration = transpiration_function)
    # land = RootSoilModel{FT}(;
    # land_args = land_args,
    land = SoilPlantHydrologyModel{FT}(;
        soil_model_type = Soil.RichardsModel{FT},
        soil_args = soil_args,
        vegetation_model_type = PlantHydraulics.PlantHydraulicsModel{FT},
        vegetation_args = plant_hydraulics_args,
    )
Y, p, cds = initialize(land)
ode! = make_ode_function(land)


# specify ICs
# previously specified from eq run: Y.vegetation.θ .=  [0.999596628794538, 0.9976596617465233]
function init_soil!(Ysoil, z, params)
    function hydrostatic_profile(
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack soil_ν, soil_vg_α, soil_vg_n, soil_vg_m, θ_r = params
        #unsaturated zone only, assumes water table starts at z_∇
        z_∇ = FT(-2)# matches zmin -- think
        S = FT((FT(1) + (soil_vg_α * (z - z_∇))^soil_vg_n)^(-soil_vg_m))
        ϑ_l = S * (soil_ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
    @show(Ysoil.soil.ϑ_l)
end
init_soil!(Y, cds.subsurface.z, land.soil.parameters)

## soil is at total ψ+z = -2.0 #m
# Want (ψ+z)_plant = (ψ+z)_soil 
p_stem_0 = (-2.0 - 0.0)
p_leaf_0 = (-2.0 - (n_stem+n_leaf)*Δz)

ϑ_l_stem_0 =
    inverse_water_retention_curve(plant_vg_α, plant_vg_n, plant_vg_m, p_stem_0, plant_ν, plant_S_s)
ϑ_l_leaf_0 =
    inverse_water_retention_curve(plant_vg_α, plant_vg_n, plant_vg_m, p_leaf_0, plant_ν, plant_S_s)
ϑ_l_0 = [ϑ_l_stem_0, ϑ_l_leaf_0]
@show(ϑ_l_0)

#vals = reverse(-Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)
#ic = [0.355686, 0.354263, 0.352855, 0.351462, 0.350083, 0.348718, 0.347368, 0.346031, 0.344708, 0.343398]

#ic_spline = Spline1D(vals, ic)
#Y.soil.ϑ_l = ic_spline.(cds.subsurface.z)
#update_aux! = make_update_aux(land) -- dont need, if we remove the returned p at t=0 might be weird
#update_aux!(p, Y, 0.0)

#sim
t0 = FT(0);
N_days = 60
tf = FT(3600 * 24 * N_days)
dt = FT(1);

sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
daily = Array(2:(3600 * 24):(N_days * 3600 * 24))
cb =
    SavingCallback((u, t, integrator) -> copy(integrator.p), sv; saveat = daily)
prob = ODEProblem(ode!, Y, (t0, tf), p);

#sol = solve(prob, RK4(), dt = dt) #, callback = cb); # 8 seconds for 180 days
update_interactions! = make_interactions_update_aux(land)
update_aux! = make_update_aux(land)
ode_function! = make_ode_function(land)

#=plots
ϕ_stem = [
    (Roots.θ_to_p.(sol.u[k].vegetation.θ[1]) .+ 0 * 9800) for
    k in 1:1:length(sol.t)
];
ϕ_leaf = [
    (Roots.θ_to_p.(sol.u[k].vegetation.θ[2]) .+ 18.5 .* 9800) for
    k in 1:1:length(sol.t)
];

lwp_leaf =
    [Roots.θ_to_p.(sol.u[k].vegetation.θ[2]) ./ 1e6 for k in 1:1:length(sol.t)];


    plot1 = plot(seconds ./ 3600 ./ 24, our_year_lwp, label = "LWP (Natan)")
    #plot!(sol.t ./ 3600 ./ 24, ϕ_stem ./ 1e6, label = "ϕ stem, MPa")
    plot!(sol.t ./ 3600 ./ 24, lwp_leaf, label = "LWP (clima) ")
    #plot!(    sol.t ./ 3600 ./ 24,ϕ_leaf ./ 1e6        label = "ϕ leaf, MPa",)
    plot!(
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
        ylim = [-3,0],
        xlabel = "t (days since Jan 1)",
        ylabel = "ϕ(MPa)",
    )


    plot!(legend = :bottomright)

    plot2 = plot(
        sol.t ./ 3600 ./ 24,
        [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)],
        label = "10cm",
        xtickfontsize = 5,
        ytickfontsize = 5,
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    )
    plot!(
        sol.t ./ 3600 ./ 24,
        [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
        label = "30cm",
    )
#    plot!(
#        sol.t ./ 3600 ./ 24,
#        [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
#        label = "50cm",
#    )
 #   plot!(
 #       sol.t ./ 3600 ./ 24,
 #       [parent(sol.u[k].soil.ϑ_l)[end - 3] for k in 1:1:length(sol.t)],
 #       label = "70cm",
 #   )
    plot!(
        sol.t ./ 3600 ./ 24,
        [mean(parent(sol.u[k].soil.ϑ_l)) for k in 1:1:length(sol.t)],
        label = "mean",
    )
    plot!(seconds ./ 3600 ./24, our_year_swc_surface, label = "Natan, surface")
    plot!(seconds ./ 3600 ./24, our_year_swc_obs_surface, label = "obs, surface")
    plot!(seconds ./ 3600 ./24, our_year_swc_column, label = "Natan, mean")
    plot!(legend = :topright, xlabel = "t (days)", ylabel = "θ soil")
    plot!(ylim = [0.3, 0.6])

    plot3 = plot(
        sol.t / 3600 ./ 24,
        precip_function.(sol.t) * 3600 * 39.37,
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
        label = "",
    )
    plot!(xlabel = "t (days)", ylabel = "Precip(in/hr)")
    plot4 = plot(
        sol.t / 3600 ./ 24,
        transpiration_function.(sol.t) * 3600 * 39.37,
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
        label = "",
    )
    plot!(xlabel = "t (days)", ylabel = "ET(in/hr)")
    plot(plot1, plot2)#, plot3, plot4, layout = 4)
    savefig("./test/LSM/ozark_2005.png")
=#
# Old code    
# Somewhat close to Weibull with C = 0.953 and B = 5.703 MPa
#const a_root = FT(0.1)
#const a_stem = a_root
#const b_root = FT(0.17 / 1e6) # Inverse Pa
#const b_stem = b_root
#const h_leaf = FT(0.005) #5mm, guess
#const h_stem = FT(18.5)# height of trunk, from Yujie's paper

# To get the equilibrium IC
#=
soil_args = (domain = soil_domain, param_set = soil_ps, boundary_conditions = FluxBC(0.0,0.0))
root_args = (domain = roots_domain, param_set = roots_ps)

land_args = (precipitation = (t) -> 0.0, transpiration = (t) -> 0.0)

land = RootSoilModel{FT}(;
                         land_args = land_args,
                         soil_model_type = Soil.RichardsModel{FT},
                         soil_args = soil_args,
                         vegetation_model_type = Roots.RootsModel{FT},
                         vegetation_args = root_args,
                         )
Y, p, cds = initialize(land)
ode! = make_ode_function(land)
p_stem_ini = -0.5e5
p_leaf_ini = -1e5
θ_stem_0 = Roots.p_to_θ(p_stem_ini)
θ_leaf_0 = Roots.p_to_θ(p_leaf_ini)
Y.vegetation.θ .= FT.([θ_stem_0, θ_leaf_0])
Y.soil.ϑ_l .= FT(0.35)
update_aux! = make_update_aux(land)
update_aux!(p,Y,0.0)

#sim
t0 = FT(0);
N_days = 300
tf = FT(3600*24*N_days)
dt = FT(1);

sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
daily = Array(2:3600*24:N_days*3600*24)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), sv; saveat = daily)
prob = ODEProblem(ode!, Y, (t0, tf), p);
sol = solve(prob, RK4(), dt = dt, callback = cb);

ϕ_soil =
    [parent(sv.saveval[k].soil.ψ .+ cds.soil.z) .* 9800 for k in 1:1:N_days]
pl = scatter()
function plot_each_day(pl, indices)
    for i in indices
        plot!(pl, ϕ_soil[i][:] ./ 1e6, parent(cds.soil.z), label = "")
    end
    plot!(
        pl,
        xlabel = "ϕ_soil (MPa)",
        ylabel = "z(m)",
        legend = :bottomright,
    )
end
plot_each_day(pl, 1:1:N_days)
=#