# Uses data and computes fluxes using a parallel flux model
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
using CurveFit

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
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/ozark/ozark_met_drivers_FLUXNET.jl",
    ),
)
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_domain.jl"),
)
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_parameters.jl"),
)
include(joinpath(climalsm_dir, "experiments/integrated/flux_helper_functions.jl"))
include(joinpath(climalsm_dir, "experiments/integrated/ozark/set_up_timeseries.jl"))

# Variables we will save
T_canopy_model = []
soil_shf = []
soil_lhf = []
canopy_lhf = []
canopy_shf = []
LW_out_model = []
LW_canopy = []
LW_soil = []
diagnostic = false
# step through the data and compute surface fluxes according to the model
steps = 120*48:1:180*48
for step in steps
    @info(step)
    if diagnostic
        # Diagnostic T canopy
        initial_guess =  [0.5 .*(T_air[step] .+ T_soil[step])]
        function flux_equality(F, x)
            T_canopy = x[1]
            q_canopy =  
                Thermodynamics.q_vap_saturation_generic(
                    thermo_params,
                    T_canopy,
                    ρ_air[step],
                    Thermodynamics.Liquid(),
                )
            ts_canopy = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_air[step], T_canopy, q_canopy)
            canopy_fluxes =  atmos_sfc_fluxes(ts_canopy, ts_in[step], h, d_sfc, u_air[step],z_0m, z_0b, surface_flux_params, thermo_params, r_canopy[step])
            lw_fluxes = land_lw_fluxes(LW_d[step], T_canopy, T_soil[step], ϵ_soil, ϵ_canopy, _σ)
            # Functions to find the roots of:
            F[1] = (-lw_fluxes.LW_canopy - SW_canopy[step] + canopy_fluxes.evaporation * _LH_v0 + canopy_fluxes.shf)
        end
        soln = nlsolve(
            flux_equality,
            initial_guess,
            ftol = 1e-8,
            iterations = 100
        )
        T_canopy = soln.zero[1]
        push!(T_canopy_model, T_canopy)
    else
        push!(T_canopy_model, T_air[step])
    end
    T_canopy = T_canopy_model[end] # valuefor this step
    q_canopy =  
        Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T_canopy,
            ρ_air[step],
            Thermodynamics.Liquid(),
            )
    ts_canopy = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_air[step], T_canopy, q_canopy)
    canopy_fluxes =  atmos_sfc_fluxes(ts_canopy, ts_in[step], h, d_sfc, u_air[step],z_0m, z_0b, surface_flux_params, thermo_params, r_canopy[step],cp_m[step])
    lw_fluxes = land_lw_fluxes(LW_d[step], T_canopy, T_soil[step], ϵ_soil, ϵ_canopy, _σ)
    ts_soil =  Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_air[step], T_soil[step], q_soil[step])
    soil_fluxes = atmos_sfc_fluxes(ts_soil, ts_in[step], h, 0.0, u_air[step],0.01, 0.001, surface_flux_params, thermo_params, r_soil[step], cp_m[step])

    # Save output
    push!(soil_shf, soil_fluxes.shf)
    push!(soil_lhf, soil_fluxes.evaporation * _LH_v0)
    push!(canopy_lhf, canopy_fluxes.evaporation * _LH_v0)
    push!(canopy_shf,  canopy_fluxes.shf)
    push!(LW_canopy, lw_fluxes.LW_canopy)
    push!(LW_out_model,lw_fluxes.LW_out)
    push!(LW_soil, lw_fluxes.LW_soil)
              
end
        
# Plotting
idx_end = length(soil_shf)
days = seconds[steps[1:idx_end]] ./ 24 ./ 3600
#Model - We do not expect conservation because T_canopy set to T_air, but the fluxes dont balance
# That is, R_n(canopy) ≠ LH(canopy)
G_model = LW_soil .+ SW_soil[steps[1:idx_end]] .- soil_shf .- soil_lhf;
net_R_model = @.( LW_d[steps[1:idx_end]] .- LW_out_model) + (SW_d -SW_out_model)[steps[1:idx_end]];
RminusG_model = net_R_model .- G_model;
HplusL_model = soil_shf .+ canopy_lhf .+ soil_lhf;

#Data
net_R = @. (LW_d-LW_OUT) + (SW_d - SW_OUT);
RminusG = net_R .- G;
HplusL = LE .+ H;

#=
# model balance
Plots.scatter(RminusG_model, HplusL_model, label = "Model", xlabel = "R-G", ylabel = "H+L", title= "Flux Conservation")
(intercept, slope) = CurveFit.linear_fit(RminusG_model, HplusL_model)
x = Array(minimum(RminusG_model):10:maximum(RminusG_model))
y = x .* slope .+intercept
Plots.plot!(x, y, label = "Fit; slope = $slope")

# Energy balance at the flux tower site:
Plots.scatter(RminusG, HplusL, label = "Data", xlabel = "R-G", ylabel = "H+L", title= "Flux Conservation")
(intercept, slope) = CurveFit.linear_fit(RminusG, HplusL)
x = Array(minimum(RminusG):10:maximum(RminusG))
y = x .* slope .+intercept
Plots.plot!(x, y, label = "Fit; slope = $slope")
=#
# Compare net rad
plt_scatter = Plots.scatter(net_R_model, net_R[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "Net Radiation")
(intercept, slope) = CurveFit.linear_fit(net_R_model, net_R[steps[1:idx_end]])
x = Array(minimum(net_R_model):10:maximum(net_R_model))
y = x .* slope .+intercept
Plots.plot!(plt_scatter, x, y, label = "Fit; slope = $slope")

# Not perfect but no immediately wrong
plt_lw = Plots.plot(0.5:0.5:24, diurnal_avg(LW_out_model, 48), label = "Model", title = "LW out", xlabel = "Hour of day")
Plots.plot!(plt_lw, 0.5:0.5:24, diurnal_avg(LW_OUT[steps[1:idx_end]], 48), label = "Data")

plt_sw = Plots.plot(0.5:0.5:24, diurnal_avg(SW_out_model[steps[1:idx_end]], 48), label = "Model", title = "SW out", xlabel = "Hour of day")
Plots.plot!(plt_sw, 0.5:0.5:24, diurnal_avg(SW_OUT[steps[1:idx_end]], 48), label = "Data")

Plots.plot(plt_scatter, plt_lw, plt_sw, layout = (1,3))
savefig("parallel_fluxes_radiation.png")


# Compare LHF - not insane
lhf = canopy_lhf .+ soil_lhf
plt_scatter = Plots.scatter(lhf, LE[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "LHF")
(intercept, slope) = CurveFit.linear_fit(lhf, LE[steps[1:idx_end]])
x = Array(minimum(lhf):10:maximum(lhf))
y = x .* slope .+intercept
Plots.plot!(plt_scatter, x, y, label = "Fit; slope = $slope")

plot_lhf = Plots.plot(0.5:0.5:24, diurnal_avg(soil_lhf, 48), label  = "Soil-Atmos", title = "LHF", xlabel = "Hour of day")
Plots.plot!(plot_lhf, 0.5:0.5:24, diurnal_avg(canopy_lhf, 48), label  = "Canopy-Atmos")
Plots.plot!(plot_lhf, 0.5:0.5:24, diurnal_avg(canopy_lhf .+ soil_lhf, 48), label  = "Land-Atmos")
Plots.plot!(plot_lhf, 0.5:0.5:24, diurnal_avg(LE[steps[1:idx_end]], 48), label = "data")
Plots.plot(plt_scatter, plot_lhf, layout =(1,2))
savefig("parallel_fluxes_lhf.png")

# G
plt_scatter = Plots.scatter(G_model, G[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "G")
(intercept, slope) = CurveFit.linear_fit(G_model, G[steps[1:idx_end]])
x = Array(minimum(G_model):10:maximum(G_model))
y = x .* slope .+intercept
Plots.plot!(plt_scatter, x, y, label = "Fit; slope = $slope")

plot_G = Plots.plot(0.5:0.5:24, diurnal_avg(G_model, 48), label = "Net soil heat flux", title = "G", xlabel = "Hour of day")
Plots.plot!(plot_G, 0.5:0.5:24, diurnal_avg(G[steps[1:idx_end]], 48), label = "Data")
Plots.plot(plt_scatter, plot_G, layout =(1,2))
savefig("parallel_fluxes_G.png")


# Looks wrong all around
shf = soil_shf
plt_scatter = Plots.scatter(shf, H[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "SHF")
(intercept, slope) = CurveFit.linear_fit(shf, H[steps[1:idx_end]])
x = Array(minimum(shf):10:maximum(shf))
y = x .* slope .+intercept
Plots.plot!(plt_scatter, x, y, label = "Fit; slope = $slope")

plot_shf = Plots.plot(0.5:0.5:24, diurnal_avg(soil_shf, 48), label  = "Soil-Atmos", title = "SHF", xlabel = "Hour of day")
Plots.plot!(plot_shf, 0.5:0.5:24, diurnal_avg(H[steps[1:idx_end]], 48), label = "data",)
Plots.plot(plt_scatter, plot_shf, layout =(1,2))
savefig("parallel_fluxes_H.png")
