# # Global run of soil model

# The code sets up and runs the soil model for on a spherical domain,
# using ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 1in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 365 d
# Timestep: 900 s
# Timestepper: ARS343
# Fixed number of iterations: 1
# Jacobian update: every new timestep
# Atmos forcing update: every 3 hours
import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts

using ClimaDiagnostics
using ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs: LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
import GeoMakie
using CairoMakie
using Dates
import NCDatasets

using Poppler_jll: pdfunite

const FT = Float64;
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
context = ClimaComms.context()
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "soil_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_prob(t0, tf, Δt; outdir = outdir, nelements = (101, 15))
    earth_param_set = LP.LandParameters(FT)
    radius = FT(6378.1e3)
    depth = FT(50)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = nelements,
        npolynomial = 1,
        dz_tuple = FT.((10.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    start_date = DateTime(2021)
    # Forcing data
    era5_artifact_path =
        ClimaLand.Artifacts.era5_land_forcing_data2021_folder_path(; context)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )
    spatially_varying_soil_params =
        ClimaLand.default_spatially_varying_soil_parameters(
            subsurface_space,
            surface_space,
            FT,
        )
    (;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo,
        NIR_albedo,
        f_max,
    ) = spatially_varying_soil_params
    soil_params = Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo = PAR_albedo,
        NIR_albedo = NIR_albedo,
    )

    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    )

    soil_args = (domain = domain, parameters = soil_params)
    soil_model_type = Soil.EnergyHydrology{FT}
    sources = (Soil.PhaseChange{FT}(),)# sublimation and subsurface runoff are added automatically
    top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(
        atmos,
        radiation,
        runoff_model,
        (:soil,),
    )
    zero_flux = Soil.HeatFluxBC((p, t) -> 0.0)
    boundary_conditions = (;
        top = top_bc,
        bottom = Soil.WaterHeatBC(;
            water = Soil.FreeDrainage(),
            heat = zero_flux,
        ),
    )
    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )


    Y, p, cds = initialize(soil)
    z = ClimaCore.Fields.coordinate_field(cds.subsurface).z
    lat = ClimaCore.Fields.coordinate_field(cds.subsurface).lat
    # This function approximates the hydrostatic equilibrium solution in
    # the vadose and unsaturated regimes by solving for ∂(ψ+z)/∂z = 0,
    # assuming the transition between the two is at a coordinate of z_∇.

    # The approximation arises because the porosity, residual water fraction,
    # and van Genuchtem parameters are spatially varying but treated constant
    # in solving for equilibrium. Finally, we make a plausible but total guess
    # for the water table depth using the TOPMODEL parameters. 
    function hydrostatic_profile(
        lat::FT,
        z::FT,
        ν::FT,
        θ_r::FT,
        α::FT,
        n::FT,
        S_s::FT,
        fmax,
    ) where {FT}
        m = 1 - 1 / n
        zmin = FT(-50.0)
        zmax = FT(0.0)

        z_∇ = FT(zmin / 5.0 + (zmax - zmin) / 2.5 * (fmax - 0.35) / 0.7)
        if z > z_∇
            S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
            ϑ_l = S * (ν - θ_r) + θ_r
        else
            ϑ_l = -S_s * (z - z_∇) + ν
        end
        return FT(ϑ_l)
    end
    vg_α = hydrology_cm.α
    vg_n = hydrology_cm.n
    Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_n, S_s, f_max)
    Y.soil.θ_i .= FT(0.0)
    T = FT(276.85)
    ρc_s =
        Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            soil_params.ρc_ds,
            soil_params.earth_param_set,
        )
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            soil_params.earth_param_set,
        )

    set_initial_cache! = make_set_initial_cache(soil)
    exp_tendency! = make_exp_tendency(soil)
    imp_tendency! = ClimaLand.make_imp_tendency(soil)
    jacobian! = ClimaLand.make_jacobian(soil)
    set_initial_cache!(p, Y, t0)

    # set up jacobian info
    jac_kwargs =
        (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )

    updateat = Array(t0:(3600 * 3):tf)
    drivers = ClimaLand.get_drivers(soil)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    # ClimaDiagnostics

    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        subsurface_space,
        outdir;
        start_date,
    )

    diags = ClimaLand.default_diagnostics(
        soil,
        start_date;
        output_writer = nc_writer,
        average_period = :monthly,
    )

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    return prob, SciMLBase.CallbackSet(driver_cb, diag_cb)
end

function setup_and_solve_problem(; greet = false)

    t0 = 0.0
    tf = 60 * 60.0 * 24 * 365
    Δt = 900.0
    nelements = (101, 15)
    if greet
        @info "Run: Global Soil Model"
        @info "Resolution: $nelements"
        @info "Timestep: $Δt s"
        @info "Duration: $(tf - t0) s"
    end

    prob, cb = setup_prob(t0, tf, Δt; nelements)

    # Define timestepper and ODE algorithm
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 3,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
    return nothing
end

setup_and_solve_problem(; greet = true);
# read in diagnostics and make some plots!
#### ClimaAnalysis ####
simdir = ClimaAnalysis.SimDir(outdir)
short_names = ["swc", "si", "sie"]
mktempdir(root_path) do tmpdir
    for short_name in short_names
        var = get(simdir; short_name)
        times = [ClimaAnalysis.times(var)[1], ClimaAnalysis.times(var)[end]]
        for t in times
            var = get(simdir; short_name)
            fig = CairoMakie.Figure(size = (600, 400))
            kwargs = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict()
            viz.heatmap2D_on_globe!(
                fig,
                ClimaAnalysis.slice(var, time = t; kwargs...),
                mask = viz.oceanmask(),
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                    :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
                ),
            )
            CairoMakie.save(joinpath(tmpdir, "$(short_name)_$t.pdf"), fig)
        end
    end
    figures = readdir(tmpdir, join = true)
    pdfunite() do unite
        run(Cmd([unite, figures..., joinpath(root_path, "figures.pdf")]))
    end
end
