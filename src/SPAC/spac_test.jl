
using Land.EarthSurface: solar_zenith_angle
using Land.ClimaCache: BetaFunction, BetaParameterG1, BetaParameterPsoil, GCO₂Mode, MedlynSM, MonoMLTreeSPAC, Soil, VanGenuchten
using Land.CanopyRadiativeTransfer: MODIS_BLUE, MODIS_EVI, MODIS_NDVI, MODIS_NIR, MODIS_NIRv, MODIS_NIRvR, MODIS_RED, TROPOMI_SIF683, TROPOMI_SIF740, OCO2_SIF759, OCO2_SIF770
using Land.SoilPlantAirContinuum: BETA, CNPP, GPP, PPAR, T_VEG, initialize!, soil_plant_air_continuum!, update!
using Land.SoilPlantAirContinuum

using OrdinaryDiffEq: ODEProblem, solve, Midpoint
using DiffEqCallbacks
# We use [ClimaCore](https://github.com/CliMA/ClimaCore.jl)
# for setting up the domain/coordinate points. While
# this infrastructure isn't really necessary for standalone simulations,
# adhering to it makes setting up coupled simulations very easy. It also
# is nice to rely on ClimaCore utilities because they have been designed
# in advance for running distributed simulations.
using ClimaCore

# We also use CLIMAParameters, which strives to ensure a common set of
# parameters across all Clima models, and to make parameter estimation
# more seamless.
import CLIMAParameters as CP

# Lastly, let's bring in the bucket model types (from ClimaLSM) that we
# will need access to.
using ClimaLSM.SPAC:
    SpacModel,
    SpacModelParameters,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    BulkAlbedo
using ClimaLSM.Domains: coordinates, LSMSingleColumnDomain
using ClimaLSM: initialize, make_update_aux, make_ode_function
# We also want to plot the solution
using Plots

FT = Float64;

# As mentioned we use CLIMAParameters for earth parameters that are
# required across models (e.g. the density of water and ice, the latent
# heat of fusion at a reference temperature, etc). The land model requires
# additional parameters as described in the text above. These two sets
# are combined in the object `BucketModelParameters` as follows:
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"));
earth_param_set = create_lsm_parameters(FT);

# Define our `BulkAlbedo` model using a constant soil and snow albedo:
# The soil albedo is a function of coordinates, which would be
# (x,y) on a plane, and (lat,lon) on a sphere. In the future, we
# will support other albedo models.
α_soil = (coordinate_point) -> FT(0.2);
α_snow = FT(0.8);
albedo = BulkAlbedo{FT}(α_snow, α_soil);
# The critical snow level setting the scale for when we interpolate between
# snow and soil albedo
σS_c = FT(0.2);
# The field capacity of the soil
W_f = FT(0.15);
# Roughness lengths (meters)
z_0m = FT(1e-2);
z_0b = FT(1e-3);
# Thermal parameters of soil
κ_soil = FT(0.7);
ρc_soil = FT(2e6);
Δt = FT(3600.0);

spac_parameters = SpacModelParameters(
    κ_soil,
    ρc_soil,
    albedo,
    σS_c,
    W_f,
    z_0m,
    z_0b,
    Δt,
    earth_param_set,
);

# Set up the model domain. At every surface coordinate point, we'll solve
# an ODE for `W` and `Ws`, and for every subsurface point, we solve for `T`.
# In coupled simulations run at the same
# resolution as the atmosphere, the bucket horizontal resolution would match the
# horizontal resolution at the lowest level of the atmosphere model. In general, however, the two
# resolutions do not need to match. Here we just set up something
# simple - an LSMSingleColumnDomain, consisting of a single column.

soil_depth = FT(3.5);
spac_domain =
    LSMSingleColumnDomain(; zlim = (-soil_depth, 0.0), nelements = 10);


# To drive the system in standalone mode,
# the user must provide
# prescribed functions of time for the water volume flux in precipitation,
#  for the net downward shortwave and longwave
# radiative energy fluxes (`SW↓, LW↓`, W/m^2),
# for the atmospheric temperature `T_a`,
# wind speed `u_a` (m/s), specific humidity `q_a`, and air density
# `ρ_a` (kg/m^3) at a reference height `h_a` (m),
# as well as for the air density `ρ_sfc` (kg/m^3)
# at the surface of the earth.

# Here we define the model drivers, starting with downward radiation.
SW_d = (t) -> eltype(t)(300);
LW_d = (t) -> eltype(t)(300);
spac_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d);

# Prescribed atmospheric variables

# Stochastic precipitation:
precip = (t) -> eltype(t)(0);
snow_precip = (t) -> eltype(t)(5e-7 * (t > 3 * 86400) * (t < 4 * 86400));
# Diurnal temperature variations:
T_atmos = (t) -> eltype(t)(275.0 + 5.0 * sin(2.0 * π * t / 86400 + 7200));
# Constant otherwise:
u_atmos = (t) -> eltype(t)(3.0);
q_atmos = (t) -> eltype(t)(0.005);
h_atmos = FT(2);
ρ_atmos = (t) -> eltype(t)(1.13);
ρ_sfc = FT(1.15);
spac_atmos = PrescribedAtmosphere(
    precip,
    snow_precip,
    T_atmos,
    u_atmos,
    q_atmos,
    ρ_atmos,
    h_atmos,
    ρ_sfc,
);

# Create a SPAC model:
spac = MonoMLTreeSPAC{FT}();

# Then, we create the model object, which contains the drivers, parameters,
# domain, and is associated with the correct differential equations
# for the bucket model:
model = SpacModel(
    spac = spac,
    parameters = spac_parameters,
    domain = spac_domain,
    atmosphere = spac_atmos,
    radiation = spac_rad,
);

# Note the holder structs for the radiation and atmosphere functions: they
# are named `Prescribed`. In coupled simulations, we would use a different
# type and rely on multiple dispatch to obtain the atmospheric and radiative
# quantitites from the coupler.

# Like all ClimaLSM models, we set up the state vector using `initialize`:
Y, p, coords = initialize(model);

# We can inspect the prognostic and auxiliary variables of the model:
ClimaLSM.prognostic_vars(model)
Y.bucket |> propertynames
# The auxiliary variables in this case are the surface temperature, the turbulent fluxes, the
# net radiation, and the surface specific humidity.
ClimaLSM.auxiliary_vars(model)
p.bucket |> propertynames


# Next is to set initial conditions. 
Y.bucket.T .= FT(270);
Y.bucket.W .= FT(0.05);
Y.bucket.Ws .= FT(0.0);
Y.bucket.σS .= FT(0.08);

# Then to create the entire right hand side function for the system
# of ordinary differential equations:
ode_function! = make_ode_function(model);

# Then set up the simulation
t0 = FT(0.0);
tf = FT(7 * 86400);
prob = ODEProblem(ode_function!, Y, (t0, tf), p);
# We need a callback to get and store the auxiliary fields, as they
# are not stored by default.
saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector);
cb = SavingCallback(
    (u, t, integrator) -> copy(integrator.p),
    saved_values;
    saveat = 0:Δt:tf,
);

sol = solve(prob, Midpoint(); dt = Δt, saveat = 0:Δt:tf, callback = cb);

# Extracting the solution from what is returned by the ODE.jl commands
# is a bit clunky right now, but we are working on hiding some of this.
# `parent` extracts the underlying data from the fields stored in
# the ClimaCore.Fields.FieldVector,
# and we loop over the solution `sol` because of how the data is stored
# within solutions returned by ODE.jl - indexed by timestep.
W = [parent(sol.u[k].bucket.W)[1] for k in 1:length(sol.t)];
Ws = [parent(sol.u[k].bucket.Ws)[1] for k in 1:length(sol.t)];
σS = [parent(sol.u[k].bucket.σS)[1] for k in 1:length(sol.t)];
T_sfc =
    [parent(saved_values.saveval[k].bucket.T_sfc)[1] for k in 1:length(sol.t)];
evaporation = [
    parent(saved_values.saveval[k].bucket.evaporation)[1] for
    k in 1:length(sol.t)
];
R_n = [parent(saved_values.saveval[k].bucket.R_n)[1] for k in 1:length(sol.t)];
# The turbulent energy flux is the sum of latent and sensible heat fluxes.
turbulent_energy_flux = [
    parent(saved_values.saveval[k].bucket.turbulent_energy_flux)[1] for
    k in 1:length(sol.t)
];