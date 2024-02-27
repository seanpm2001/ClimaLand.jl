using CairoMakie
using Statistics
using ArtifactWrappers
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil

import ClimaLand
import ClimaLand.Parameters as LP
context = ClimaComms.context()
outdir = joinpath(pkgdir(ClimaLand), "experiments/standalone/Soil/artifacts")
!ispath(outdir) && mkpath(outdir)
FT = Float64
radius = FT(6378.1e3);
depth = FT(50)
domain = ClimaLand.Domains.SphericalShell(;
                                          radius = radius,
                                          depth = depth,
                                          nelements = (101, 15),
                                          npolynomial = 1,
                                          dz_tuple = FT.((10.0, 0.1)),
                                          );
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface
# Read in f_max data and land sea mask
infile_path = "/Users/katherinedeck/Desktop/code/ClimaLand.jl/means_2.5_new.nc"
regrid_dirpath = "regrid"
outfile_root = "static_data_cgll"
ClimaLand.Regridder.hdwrite_regridfile_rll_to_cgll(
    FT,
    regrid_dirpath,
    infile_path,
    ["fmax", "landsea_mask"],
    surface_space,
    outfile_root;
    mono = true,
)

file_info =
    ClimaLand.FileReader.FileInfo(infile_path, regrid_dirpath, ["fmax", "landsea_mask"], outfile_root, [], [])
data =  ClimaLand.FileReader.PrescribedDataStatic{typeof(file_info)}(file_info)
f_max = ClimaLand.FileReader.get_data_at_date(data,
                                             surface_space,
                                             "fmax",
                                              )
mask = ClimaLand.FileReader.get_data_at_date(data,
                                             surface_space,
                                             "landsea_mask",
                                             )
oceans_to_zero(field, mask) = mask > 0.5 ? field : 0
f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4/1000) # m/s
runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;f_over=f_over, f_max=f_max, R_sb=R_sb)
vg_α = ClimaCore.Fields.ones(subsurface_space) .* FT(0.2)
hydrology_cm = map(vg_α) do (α)
    FT = typeof(α)
    ClimaLand.Soil.vanGenuchten{FT}(;α = α, n = α +FT(2))
end
θ_r = ClimaCore.Fields.zeros(subsurface_space)
ν =  ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.5)
K_sat =  ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-6)
S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)

soil_params = ClimaLand.Soil.RichardsParameters(; hydrology_cm = hydrology_cm,
                                                ν = ν,
                                                K_sat = K_sat,
                                                S_s = S_s,
                                                θ_r = θ_r)

function precip_function(t)
    return -1e-6
end
precip = ClimaLand.TimeVaryingInput(precip_function)
atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip)
bottom_bc = ClimaLand.Soil.FluxBC((p,t) -> 0.0)
bc = (; top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(atmos, runoff_model), bottom =  bottom_bc)
model = ClimaLand.Soil.RichardsModel{FT}(; parameters = soil_params, domain = domain,
                                         boundary_conditions = bc, sources = (), lateral_flow=false)
Y,p,t = initialize(model)
z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
function hydrostatic_profile(lat::FT, z::FT, ν::FT, θ_r::FT, α::FT, n::FT, S_s::FT, fmax) where {FT}
    m = 1-1/n
    zmin = FT(-50.0)
    zmax = FT(0.0)
    
    z_∇ = FT( zmin / 5.0+(zmax - zmin) / 2.5 * (fmax - 0.35)/0.7)
    if z > z_∇
        S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
        ϑ_l = S * (ν - θ_r) + θ_r
    else
        ϑ_l = -S_s * (z - z_∇) + ν
    end
    return FT(ϑ_l)
end
t0 = 0.0
tf = 360.0*10
dt = 360.0*2
Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_α .+2, S_s, f_max)
@. Y.soil.ϑ_l = oceans_to_zero(Y.soil.ϑ_l, mask)
set_initial_cache! = make_set_initial_cache(model)
exp_tendency! = make_exp_tendency(model);
imp_tendency! = ClimaLand.make_imp_tendency(model);
update_jacobian! = ClimaLand.make_update_jacobian(model);

set_initial_cache!(p,Y,t0)
stepper = CTS.ARS111()
norm_condition = CTS.MaximumError(FT(1e-8))
conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 50,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

# set up jacobian info
jac_kwargs = (;
              jac_prototype = RichardsTridiagonalW(Y),
              Wfact = update_jacobian!,
              )

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
saveat = Array(t0:dt:tf)
sv = (;
      t = Array{Float64}(undef, length(saveat)),
      saveval = Array{NamedTuple}(undef, length(saveat)),
      )
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
updateat = deepcopy(saveat)
updatefunc = ClimaLand.make_update_drivers(atmos, nothing)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
@time sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = dt, callback = cb)

longpts = range(-180.0, 180.0, 101)
latpts = range(-90.0, 90.0, 101)
hcoords = [ClimaCore.Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]
remapper = ClimaCore.Remapping.Remapper(hcoords, nothing, surface_space)

h∇ = [ClimaCore.Remapping.interpolate(remapper,sv.saveval[k].soil.h∇) for k in 1:length(sol.t)]
R_s = [ClimaCore.Remapping.interpolate(remapper,sv.saveval[k].soil.R_s) for k in 1:length(sol.t)]
R_ss = [ClimaCore.Remapping.interpolate(remapper,sv.saveval[k].soil.R_ss) for k in 1:length(sol.t)]
θ_sfc = [ClimaCore.Remapping.interpolate(remapper,ClimaLand.Soil.get_top_surface_field(sol.u[k].soil.ϑ_l, surface_space)) for k in 1:length(sol.t)]
int_ϑ = [similar(sv.saveval[k].soil.h∇) for k in 1:length(sol.t)]
[ClimaCore.Operators.column_integral_definite!(int_ϑ[k], sol.u[k].soil.ϑ_l) for k in 1:length(sol.t)];
int_ϑ = [ClimaCore.Remapping.interpolate(remapper,int_ϑ[k]) for k in 1:length(sol.t)]

for (i, (field, field_name)) in enumerate(
    zip(
        [h∇, R_s, R_ss, θ_sfc, int_ϑ],
        ["Water table thickness(m)", "Surface Runoff (m/s)", "Subsurface Runoff (m/s)", "θ_sfc (m³/m³)", "∫ϑdz (m)"],
    ),
)
    fig = Figure(size = (2000, 1000))
    ax = Axis(
        fig[1, 1],
        xlabel = "Longitude",
        ylabel = "Latitude",
        title = field_name,
    )
    a2 = Axis(
        fig[2, 1],
        xlabel = "Longitude",
        title = "Δ",
    )
    CairoMakie.heatmap!(ax,
                        longpts,
                        latpts,
                        field[end])
    
    CairoMakie.heatmap!(
        ax2
        longpts,
        latpts,
        field[end]-field[1],
    )
    outfile  = joinpath(
        outdir,
        string("heatmap_", field_name, ".png"),
    )
    CairoMakie.save(outfile, fig)
end
