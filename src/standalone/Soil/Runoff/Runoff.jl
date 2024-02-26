module Runoff
using ClimaCore.Operators: column_integral_definite!
using ClimaCore.Fields: Field
using ClimaLand
import ClimaLand: source!
using ..ClimaLand.Soil: AbstractSoilSource, AbstractSoilModel
export soil_surface_infiltration,
    TOPMODELRunoff,
    AbstractRunoffModel,
    TOPMODELSubsurfaceRunoff,
    subsurface_runoff_source

"""
    AbstractRunoffModel
The abstract type for soil runoff models to be used with
the following boundary condition types:
- `ClimaLand.Soil.AtmosDrivenFluxBC`
- `ClimaLand.Soil.RichardsAtmosDrivenFluxBC`,
and for these functions:
-`ClimaLand.Soil.soil_surface_infiltration`
- `ClimaLand.Soil.subsurface_runoff_source`
- `ClimaLand.source!`.
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
    soil_surface_infiltration(::NoRunoff, precip, _...)
A function which computes the infiltration into the soil
 for the default of `NoRunoff`.

Returns the precipitation.
"""
soil_surface_infiltration(::NoRunoff, precip, _...) = precip

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

struct TOPMODELRunoff{FT <: AbstractFloat, F <: Field} <: AbstractRunoffModel
    f_over::FT
    f_max::F
    subsurface_source::TOPMODELSubsurfaceRunoff{FT}
end

"""
TOPMODEL infiltration
"""
function soil_surface_infiltration(
    runoff::TOPMODELRunoff,
    precip,
    Y,
    p,
    surface_space,
    model::AbstractSoilModel,
)

    # we need this at the surface only
    surface_space = axes(Δz_top)
    (; hydrology_cm, K_sat, ν, θ_r) = model.parameters
    K_eff = p.soil.K ./  ClimaLand.Soil.hydraulic_conductivity(
                hydrology_cm,
                K_sat,
                effective_saturation(ν, Y.soil.ϑ_l, θ_r),
            )
    infiltration_capacity = get_top_surface_field(-K_eff, surface_space)
    #needs to be by column
    h∇ = ClimaCore.Fields.zeros(surface_space)
    water_content = Y.soil.ϑ_l .+ Y.soil.θ_i
    column_integral_definite!(h∇, heaviside.(water_content .- ν .-0.001))# water table thickness
    fsat = @. runoff.f_max*exp(-1/2*runoff.f_over*(model.domain.depth - h∇)).*heaviside(h∇) # only apply where there is a water table
    return @. (1 - fsat) * max(infiltration_capacity, precip)
end

"""
    TOPMODEL sink term for baseflow
"""
function ClimaLand.source!(
    dY,
    src::TOPMODELSubsurfaceRunoff,
    Y,
    p,
    model::AbstractSoilModel,
)
    ν = model.parameters.ν
    #needs to be by column
    h∇ = ClimaCore.Fields.zeros(surface_space)
    water_content = Y.soil.ϑ_l .+ Y.soil.θ_i
    h∇ = column_integral_definite!(h∇, sum(heaviside.(water_content .- ν .-0.001)))# water table thickness
    @. dY.soil.ϑ_l -= src.R_sb * exp(-f_over*(model.domain.depth - h∇)) / h∇ * heaviside.(water_content .- ν- 0.001) # apply only to saturated layers

end

subsurface_runoff_source(runoff::TOPMODELRunoff) = runoff.subsurface_source
end
