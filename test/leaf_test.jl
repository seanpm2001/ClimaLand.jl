using Photosynthesis, PlantHydraulics, StomataModels, WaterPhysics

using Statistics
using DifferentialEquations
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
using Plots 

FT = Float64 

# Create a dummy BallBerry C3 Leaf model:
params = ClimaLSM.Leaf.LandLeafModel(
    C3CLM(FT),
    ESMBallBerry{FT}(g0=0.0),
    LeafHydraulics{FT}()
);

leaf = ClimaLSM.Leaf.LeafModel{FT}(;
    param_set = params,
    clayer    = CanopyLayer{FT}(n_leaf = 1, g_min = -10.0, τ_esm = 600),
    environ   = AirLayer{FT}(),
);

# initialize manually:
Y = initialize_prognostic(leaf,[1.0])
p = initialize_auxiliary(leaf,[1.0])

# Set g_sw to 0 at the beginning:
Y.leaf.g_sw .= 0.001;

leaf_ode! = make_ode_function(leaf)

saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values)

t0 = FT(0);tf = FT(15000);
dt = FT(1);
prob = ODEProblem(leaf_ode!, Y, (t0, tf), p);
sol = solve(prob, Euler(), dt = dt,callback = cb);

# Get gsw out of solution again:
gsw = [(sol.u[k].leaf.g_sw[1]) for k in 1:1:length(sol.t)]
plot(sol.t, gsw)

plot(saved_values.t, [saved_values.saveval[i].leaf.Ci[1]   for i in eachindex(saved_values.t)])
plot(saved_values.t, [saved_values.saveval[i].leaf.An[1]   for i in eachindex(saved_values.t)])
plot(saved_values.t, [saved_values.saveval[i].leaf.NPQ[1]  for i in eachindex(saved_values.t)])
plot(saved_values.t, [saved_values.saveval[i].leaf.APAR[1] for i in eachindex(saved_values.t)])