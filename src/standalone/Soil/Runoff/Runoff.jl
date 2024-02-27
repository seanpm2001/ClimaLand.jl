module Runoff
using ClimaCore
using ClimaLand
import ClimaLand: source!
using ..ClimaLand.Soil: AbstractSoilSource, AbstractSoilModel, RichardsModel, EnergyHydrology
export soil_surface_infiltration,
    TOPMODELRunoff,
    NoRunoff,
    AbstractRunoffModel,
    TOPMODELSubsurfaceRunoff,
    subsurface_runoff_source,
    soil_surface_infiltration,
    topmodel_ss_flux
    

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
    f_over::FT
end

struct TOPMODELRunoff{FT <: AbstractFloat, F <: ClimaCore.Fields.Field} <: AbstractRunoffModel
    f_over::FT
    f_max::F
    subsurface_source::TOPMODELSubsurfaceRunoff{FT}
end

function TOPMODELRunoff{FT}(; f_over::FT, f_max::F, R_sb::FT) where {FT, F}
    subsurface_source = TOPMODELSubsurfaceRunoff{FT}(R_sb, f_over)
    return TOPMODELRunoff{FT, F}(f_over, f_max, subsurface_source)
end


"""
TOPMODEL infiltration
"""
function soil_surface_infiltration(
    runoff::TOPMODELRunoff,
    precip,
    Y,
    p,
    model,
)
    flux_ic = soil_infiltration_capacity_flux(model, Y, p)
    return @. topmodel_surface_infiltration(p.soil.h∇, runoff.f_max, runoff.f_over, model.domain.depth - p.soil.h∇, flux_ic, precip)
end
function topmodel_surface_infiltration(h∇, f_max, f_over, z∇, f_ic, precip)
    f_sat = f_max*exp(-f_over/2*z∇) # This will be extremely small if the depth is ~50m and h∇ = 0
    return (1-f_sat)*max(f_ic, precip)
end

function soil_infiltration_capacity_flux(model::RichardsModel, Y, p)
    return ClimaLand.Soil.get_top_surface_field(-1 .* model.parameters.K_sat, model.domain.space.surface)
end

function soil_infiltration_capacity_flux(model::EnergyHydrology, Y, p)
    (; Ksat, θ_r, Ω, γ, γT_ref) = model.parameters
    K_eff = @. Ksat * impedance_factor(Y.soil.θ_i / (p.soil.θ_l + Y.soil.θ_i - θ_r), Ω) * viscosity_factor(p.soil.T, γ, γT_ref)
    return ClimaLand.Soil.get_top_surface_field(-1 .* Keff, model.domain.space.surface)
end


"""
    TOPMODEL sink term for baseflow
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::TOPMODELSubsurfaceRunoff,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::RichardsModel,
)
    FT = eltype(Y.soil.ϑ_l)
    ν = model.parameters.ν
    h∇ = p.soil.h∇
    @. dY.soil.ϑ_l -= p.soil.R_ss / max(h∇, eps(FT)) * ClimaLand.heaviside(Y.soil.ϑ_l  - ν) # apply only to saturated layers
end

"""
    TOPMODEL sink term for baseflow
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::TOPMODELSubsurfaceRunoff,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::EnergyHydrology,
)
    FT = eltype(Y.soil.ϑ_l)
    ν = model.parameters.ν
    h∇ = p.soil.h∇
    @. dY.soil.ϑ_l -= p.soil.R_ss / max(h∇, eps(FT)) * ClimaLand.heaviside(Y.soil.ϑ_l + Y.soil.θ_i  - ν) # apply only to saturated layers
end

function topmodel_ss_flux(R_sb::FT, f_over::FT, z∇::FT) where {FT}
    return R_sb * exp(-f_over*z∇)
end

subsurface_runoff_source(runoff::TOPMODELRunoff) = runoff.subsurface_source
end
