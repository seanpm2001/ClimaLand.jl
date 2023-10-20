# Uses data and computes fluxes using a CLM-style model
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
T_airspace_model = []
q_airspace_model = []
T_canopy_model = []
shf = []
lhf = []
soil_shf = []
soil_lhf = []
canopy_lhf = []
canopy_shf = []
ustar_model = []
LW_out_model = []
LW_canopy = []
LW_soil = []
L_MO = []
converged = []

# step through the data and compute surface fluxes according to the model
steps = 120*48:1:180*48
# Errors on near step 7613
for step in 120*48:1:5911#7612
    @info(step)
    # Tairspace, qairspace, Tcanopy, ustar, shf, lhf, LMO. Only the first three are optimized, so we set F[4:7] = 0.
    initial_guess = [0.5 .*(T_air[step] .+ T_soil[step]), 0.5 .*(q_air[step] .+ q_soil[step]), 0.5 .*(T_air[step] .+ T_soil[step]), u_air[step], 0.0, 0.0, 0.0]
    function flux_equality(F, x)
        # Compute the surface fluxes between the airspace and the atmosphere. CLM 2.5.100 and CLM2.5.87. See discussion right after 2.5.92.
        # The surface corresponds to the canopy airspace, with values at height z_0+d
        T_airspace = x[1]
        q_airspace = x[2]
        T_canopy = x[3]
        
        ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_air[step], T_airspace, q_airspace)
        airspace_atmos_fluxes =  atmos_airspace_fluxes(ts_sfc, ts_in[step], h, d_sfc, u_air[step],z_0m, z_0b, surface_flux_params, cp_m[step])
        ustar = airspace_atmos_fluxes.ustar
        soil_fluxes = airspace_soil_fluxes(ustar, T_airspace, q_airspace, W[step], T_soil[step], q_soil[step], r_soil[step], ρ_air[step], cp_m[step])
        canopy_fluxes = airspace_canopy_fluxes(ustar, T_airspace, q_airspace, T_canopy, ρ_air[step], cp_m[step], LAI[step], AI[step], r_canopy[step])
        lw_fluxes = land_lw_fluxes(LW_d[step], T_canopy, T_soil[step], ϵ_soil, ϵ_canopy, _σ)
            
        # Functions to find the roots of:
        F[1] = canopy_fluxes.shf + soil_fluxes.shf - airspace_atmos_fluxes.shf
        F[2] = (canopy_fluxes.transpiration + soil_fluxes.evaporation - airspace_atmos_fluxes.evaporation) * _LH_v0
        F[3] = (-lw_fluxes.LW_canopy - SW_canopy[step] + canopy_fluxes.transpiration * _LH_v0 + canopy_fluxes.shf)
        F[4] = 0.0
        F[5] = 0.0
        F[6] = 0.0
        F[7] = 0.0

        # Store latest values from this iteration.
        x[4] = ustar
        x[5] = airspace_atmos_fluxes.shf
        x[6] = airspace_atmos_fluxes.lhf
        x[7] = airspace_atmos_fluxes.L_MO
    end
    soln = nlsolve(
        flux_equality,
        initial_guess,
        ftol = 1e-8,
        iterations = 100
    )
    # If we converge according to either metric (F vs x)
    if ~soln.f_converged & ~soln.x_converged
        push!(converged, true)
    else
        push!(converged, false)
    end
    
    T_airspace = soln.zero[1]
    q_airspace = soln.zero[2]
    T_canopy = soln.zero[3]
    ustar= soln.zero[4]
    
    soil_fluxes = airspace_soil_fluxes(ustar, T_airspace, q_airspace, W[step], T_soil[step], q_soil[step], r_soil[step], ρ_air[step], cp_m[step])
    canopy_fluxes = airspace_canopy_fluxes(ustar, T_airspace, q_airspace, T_canopy, ρ_air[step], cp_m[step], LAI[step], AI[step], r_canopy[step])
    lw_fluxes = land_lw_fluxes(LW_d[step], T_canopy, T_soil[step], ϵ_soil, ϵ_canopy, _σ)
    
    
    # Save output
    push!(T_airspace_model, T_airspace)
    push!(q_airspace_model, q_airspace)
    push!(ustar_model, ustar)
    push!(shf, soln.zero[5])
    push!(lhf, soln.zero[6])
    push!(L_MO, soln.zero[7])
    
    push!(soil_shf, soil_fluxes.shf)
    push!(soil_lhf, soil_fluxes.evaporation * _LH_v0)
    push!(canopy_lhf, canopy_fluxes.transpiration * _LH_v0)
    push!(canopy_shf,  canopy_fluxes.shf)
    push!(LW_canopy, lw_fluxes.LW_canopy)
    push!(T_canopy_model, T_canopy)
    push!(LW_out_model,lw_fluxes.LW_out)
    push!(LW_soil, lw_fluxes.LW_soil)
end

# Plotting
idx_end = length(T_airspace_model)
days = seconds[steps[1:idx_end]] ./ 24 ./ 3600

#Model
G_model = LW_soil .+ SW_soil[steps[1:idx_end]] .- soil_shf .- soil_lhf
net_R_model = @.( LW_d[steps[1:idx_end]] .- LW_out_model) + (SW_d -SW_out_model)[steps[1:idx_end]]
RminusG_model = net_R_model .- G_model
HplusL_model = shf .+ lhf

#Data
net_R = @. (LW_d-LW_OUT) + (SW_d - SW_OUT)
RminusG = net_R .- G
HplusL = LE .+ H


# model balance
plt_scatter_model = Plots.scatter(RminusG_model, HplusL_model, label = "Model", xlabel = "R-G", ylabel = "H+L", title= "Flux Conservation")
(intercept, slope) = CurveFit.linear_fit(RminusG_model, HplusL_model)
x = Array(minimum(RminusG_model):10:maximum(RminusG_model))
y = x .* slope .+intercept
Plots.plot!(plt_scatter_model, x, y, label = "Fit; slope = $slope")

# Energy balance at the flux tower site:
plt_scatter_data = Plots.scatter(RminusG, HplusL, label = "Data", xlabel = "R-G", ylabel = "H+L", title= "Flux Conservation")
(intercept, slope) = CurveFit.linear_fit(RminusG, HplusL)
x = Array(minimum(RminusG):10:maximum(RminusG))
y = x .* slope .+intercept
Plots.plot!(plt_scatter_data, x, y, label = "Fit; slope = $slope")
Plots.plot(plt_scatter_model, plt_scatter_data, layout = (1,2))
savefig("CLM_conservation.png")


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
savefig("CLM_fluxes_radiation.png")


# Compare LHF - not insane
plt_scatter = Plots.scatter(lhf, LE[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "LHF")
(intercept, slope) = CurveFit.linear_fit(lhf, LE[steps[1:idx_end]])
x = Array(minimum(lhf):10:maximum(lhf))
y = x .* slope .+intercept
Plots.plot!(plt_scatter, x, y, label = "Fit; slope = $slope")

plot_lhf = Plots.plot(0.5:0.5:24, diurnal_avg(soil_lhf, 48), label  = "Soil-Airspace", title = "LHF", xlabel = "Hour of day")
Plots.plot!(plot_lhf, 0.5:0.5:24, diurnal_avg(canopy_lhf, 48), label  = "Canopy-Airspace")
Plots.plot!(plot_lhf, 0.5:0.5:24, diurnal_avg(canopy_lhf .+ soil_lhf, 48), label  = "Airspace-Atmos")
Plots.plot!(plot_lhf, 0.5:0.5:24, diurnal_avg(LE[steps[1:idx_end]], 48), label = "data", legend = :topleft)
Plots.plot(plt_scatter, plot_lhf, layout =(1,2))
savefig("CLM_fluxes_lhf.png")


# Compare R-L - still looks good
#Plots.scatter(net_R_model .- lhf, (net_R .- LE)[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "Rn- LHF")
#(intercept, slope) = CurveFit.linear_fit(net_R_model .- lhf, (net_R .- LE)[steps[1:idx_end]])
#x = Array(minimum(net_R_model .- lhf):10:maximum(net_R_model .- lhf))
#y = x .* slope .+intercept
#Plots.plot!(x, y, label = "Fit; slope = $slope")



# G
plt_scatter = Plots.scatter(G_model, G[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "G")
(intercept, slope) = CurveFit.linear_fit(G_model, G[steps[1:idx_end]])
x = Array(minimum(G_model):10:maximum(G_model))
y = x .* slope .+intercept
Plots.plot!(plt_scatter, x, y, label = "Fit; slope = $slope")

plot_G = Plots.plot(0.5:0.5:24, diurnal_avg(G_model, 48), label = "Net soil heat flux", title = "G", xlabel = "Hour of day")
Plots.plot!(plot_G, 0.5:0.5:24, diurnal_avg(G[steps[1:idx_end]], 48), label = "Data")
Plots.plot(plt_scatter, plot_G, layout =(1,2))
savefig("CLM_fluxes_G.png")


# Looks wrong all around
plt_scatter = Plots.scatter(shf, H[steps[1:idx_end]], xlabel = "Model", ylabel = "Data", title= "SHF")
(intercept, slope) = CurveFit.linear_fit(shf, H[steps[1:idx_end]])
x = Array(minimum(shf):10:maximum(shf))
y = x .* slope .+intercept
Plots.plot!(plt_scatter, x, y, label = "Fit; slope = $slope")

plot_shf = Plots.plot(0.5:0.5:24, diurnal_avg(soil_shf, 48), label  = "Soil-Airspace", title = "SHF", xlabel = "Hour of day")
Plots.plot!(plot_shf, 0.5:0.5:24, diurnal_avg(canopy_shf, 48), label  = "Canopy-Airspace")
Plots.plot!(plot_shf, 0.5:0.5:24, diurnal_avg(canopy_shf .+ soil_shf, 48), label  = "Airspace-Atmos")
Plots.plot!(plot_shf, 0.5:0.5:24, diurnal_avg(H[steps[1:idx_end]], 48), label = "data",)
Plots.plot(plt_scatter, plot_shf, layout =(1,2))
savefig("CLM_fluxes_H.png")


# Weird spikes but not wrong immediately. due to MO solve?
plt_temp = Plots.plot( 0.5:0.5:24, diurnal_avg(T_airspace_model, 48), label = "Airspace", title = "Temperature")
Plots.plot!(plt_temp,  0.5:0.5:24, diurnal_avg(T_air[steps[1:idx_end]], 48), label = "Atmos")
Plots.plot!(plt_temp,  0.5:0.5:24, diurnal_avg(T_soil[steps[1:idx_end]], 48), label = "Soil")
Plots.plot!(plt_temp,  0.5:0.5:24, diurnal_avg(T_canopy_model, 48), label = "Canopy")
savefig("CLM_temperatures.png")

# same spikes
plot_q = Plots.plot( 0.5:0.5:24, diurnal_avg(q_airspace_model, 48), label = "Airspace", title = "Specific Humidity")
Plots.plot!(plot_q, 0.5:0.5:24, diurnal_avg(q_air[steps[1:idx_end]], 48), label = "Atmos")
Plots.plot!(plot_q,  0.5:0.5:24, diurnal_avg(q_soil[steps[1:idx_end]], 48), label = "Soil")
savefig("CLM_q.png")

# Spikes when ustar -> 0?
plot_ustar = Plots.plot(days, ustar_model, title = "Ustar")
