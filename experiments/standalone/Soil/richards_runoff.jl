using Plots
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

FT = Float64
radius = FT(6378.1e3);
depth = FT(50)
domain = ClimaLand.Domains.SphericalShell(;
                                          radius = radius,
                                          depth = depth,
                                          nelements = (101, 15),
                                          npolynomial = 1,
                                          dz_tuple = FT.((5.0, 0.5)),
                                          );
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface
# Read in f_max data and land sea mask
infile_path = "/Users/katherinedeck/Desktop/code/ClimaLand.jl/means_2.5.nc"
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
    FileInfo(infile_path, regrid_dirpath, varnames, outfile_root, [], [])
data =  PrescribedDataStatic{typeof(file_info)}(file_info)



f_max = ClimaCore.Fields.ones(surface_space).*FT(0.5)
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
    return 1e-6
end
precip = ClimaLand.TimeVaryingInput(precip_function)
bottom_bc = ClimaLand.Soil.FluxBC((p,t) -> 0.0)
bc = (; top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(precip, runoff_model), bottom =  bottom_bc)
model = ClimaLand.Soil.RichardsModel{FT}(; parameters = soil_params, domain = domain,
                                         boundary_conditions = bc, sources = ())
Y,p,t = initialize(model)
z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
function hydrostatic_profile(lat::FT, z::FT, ν::FT, θ_r::FT, α::FT, n::FT, S_s::FT) where {FT}
    m = 1-1/n
    zmin = FT(-50.0)
    zmax = FT(0.0)
    
    z_∇ = FT(
        zmin / 2.0 +
        (zmax - zmin) / 20.0 * sin(π * 2 * lat / 180.0),
    )
    if z > z_∇
        S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
        ϑ_l = S * (ν - θ_r) + θ_r
    else
        ϑ_l = -S_s * (z - z_∇) + ν
    end
    return FT(ϑ_l)
end
t0 = 0.0
tf = 720.0*2
dt = 720.0
Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_α .+2, S_s)
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

sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = dt, callback = saving_cb)

h∇ = [sv.saveval[k].soil.h∇ for k in 1:length(sol.t)]
R_s = [sv.saveval[k].soil.R_s for k in 1:length(sol.t)]
R_ss = [sv.saveval[k].soil.R_ss for k in 1:length(sol.t)]
θ_sfc = [ClimaLand.Soil.get_top_surface_field(sol.u[k].soil.ϑ_l, domain.space.surface) for k in 1:length(sol.t)]
int_ϑ = [similar(θ_sfc[k]) for k in 1:length(sol.t)]
[ClimaCore.Operators.column_integral_definite!(int_ϑ[k], sol.u[k].soil.ϑ_l) for k in 1:length(sol.t)];
