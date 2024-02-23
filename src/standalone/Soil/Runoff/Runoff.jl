module Runoff
using ClimaLSM
using ClimaLSM.Regridder: regrid_netcdf_to_field
import ClimaLSM: source!
using ..ClimaLSM.Soil: AbstractSoilSource, AbstractSoilModel
import ..ClimaLSM.Soil: set_initial_parameter_field!
export soil_surface_infiltration,
    TOPMODELRunoff,
    AbstractRunoffModel,
    TOPMODELSubsurfaceRunoff,
    subsurface_runoff_source

"""
    AbstractRunoffModel
The abstract type for soil runoff models to be used with
the following boundary condition types:
- `ClimaLSM.Soil.AtmosDrivenFluxBC`
- `ClimaLSM.Soil.RichardsAtmosDrivenFluxBC`,
and for these functions:
-`ClimaLSM.Soil.soil_surface_infiltration`
- `ClimaLSM.Soil.subsurface_runoff_source`
- `ClimaLSM.source!`.
Please see the documentation for these for more details.
The model should specify the subsurface runoff sink term as well
as the surface runoff implementation.
"""
abstract type AbstractRunoffModel end

"""
    NoRunoff <: AbstractRunoffModel
A concrete type of soil runoff model; the 
default choice which does not include the 
effects of runoff.
"""
struct NoRunoff <: AbstractRunoffModel end

"""
    soil_surface_infiltration(::NoRunoff, net_water_flux, _...)
A function which computes the infiltration into the soil
 for the default of `NoRunoff`.
If `net_water_flux = P+E`, where `P` is the precipitation and `E`
is the evaporation (both negative if towards the soil), 
this returns `P+E` as the water boundary flux for the soil.
"""
soil_surface_infiltration(::NoRunoff, net_water_flux, _...) = net_water_flux

"""
    subsurface_runoff_source(runoff::AbstractRunoffModel)::Union{Nothing, AbstractSoilSource} 
A function which returns the soil source for the runoff model 
`runoff`; the default returns nothing in which case no source is added.
"""
subsurface_runoff_source(
    runoff::AbstractRunoffModel,
)::Union{Nothing, AbstractSoilSource} = nothing

# TOPMODEL

struct TOPMODELSubsurfaceRunoff{FT} <: AbstractSoilSource{FT}
    R_sb::FT
end

struct TOPMODELRunoff{FT <: AbstractFloat} <: AbstractRunoffModel
    f_over::FT
    f_max::F
    subsurface_source::TOPMODELSubsurfaceRunoff{FT}
end

"""
TOPMODEL infiltration
"""
function soil_surface_infiltration(
    runoffmodel::TOPMODELRunoff,
    net_water_flux,
    Y,
    p,
    soil_parameters,
    surface_space
)

    # we need this at the surface only
    surface_space = axes(Δz_top)
    infiltration_capacity = get_top_surface_field(-soil_parameters.K_sat, surface_space)
    z∇ = sum(ClimaCore.Fields.ones(axes(Y.soil.ϑ_l))) - sum(heaviside.(Y.soil.ϑ_l .- ν .-0.001))
    # net_water_flux is negative if towards the soil; take smaller in magnitude -> max
    # net_water_flux is positive if away from soil -> use as BC.
    return @. (1 - runoffmodel.f_max*exp(-1/2*runoffmodel.f_over*z∇)) * max(infiltration_capacity, net_water_flux)
    # check limit - should be zero if z∇ is equal to soil depth!
end

"""
    TOPMODEL sink term for baseflow - standin
"""
function ClimaLSM.source!(
    dY,
    src::TOPMODELSubsurfaceRunoff,
    Y,
    p,
    model::AbstractSoilModel,
)
    ν = model.parameters.ν
    h∇ = sum(heaviside.(Y.soil.ϑ_l .- ν .-0.001))# water table thickness
    z∇ = model.domain.depth - h∇
    @. dY.soil.ϑ_l -= src.R_sb * exp(-f_over*z∇) / h∇ * heaviside.(Y.soil.ϑ_l .- ν- 0.001) # apply only to saturated layers

end

subsurface_runoff_source(runoff::TOPMODELRunoff) = runoff.subsurface_source
end
