"""
    IntegratedAtmosDrivenFluxBC{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
        R <: AbstractRunoffModel,
        C <: Tuple
    } <: AbstractAtmosDrivenFluxBC


$(DocStringExtensions.FIELDS)
"""
struct IntegratedAtmosDrivenFluxBC{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
    R <: ClimaLand.Soil.AbstractRunoffModel,
    C,
} <: AbstractAtmosDrivenFluxBC
    "The atmospheric conditions driving the model"
    atmos::A
    "The radiative fluxes driving the model"
    radiation::B
    "The runoff model. The default is no runoff."
    runoff::R
    "Component list"
    components::C
end


function ClimaLand.Soil.sublimation_source(
    bc::IntegratedAtmosDrivenFluxBC{AbstractAtmosphericDrivers{FT}},# dispatch off of whether :snow and :soil are in components?
) where {FT}
    if :snow ∈ bc.components 
        return SoilSublimationwithSnow{FT}()
    else
        return SoilSublimation{FT}()
    end
    
end

"""
    SoilSublimationwithSnow{FT} <: AbstractSoilSource{FT}

Soil Sublimation source type. Used to defined a method
of `ClimaLand.source!` for soil sublimation with snow present.
"""
struct SoilSublimationwithSnow{FT} <: ClimaLand.Soil.AbstractSoilSource{FT} end

"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::SoilSublimationwithSnow{FT},
             Y::ClimaCore.Fields.FieldVector,
             p::NamedTuple,
             model
             )

Updates dY.soil.θ_i in place with a term due to sublimation; this only affects
the surface layer of soil.

"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::SoilSublimationwithSnow{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
) where {FT}
    _ρ_i = FT(LP.ρ_cloud_ice(model.parameters.earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(model.parameters.earth_param_set))
    z = model.domain.fields.z
    Δz_top = model.domain.fields.Δz_top # this returns the center-face distance, not layer thickness
    @. dY.soil.θ_i +=
        -p.soil.turbulent_fluxes.vapor_flux_ice *
        (1 - p.snow.snow_cover_fraction) *
        _ρ_l / _ρ_i * heaviside(z + 2 * Δz_top) # only apply to top layer, recall that z is negative
end

function soil_boundary_fluxes!(
    bc::IntegratedAtmosDrivenFluxBC,
    boundary::ClimaLand.TopBoundary,
    soil::EnergyHydrology{FT},
    Δz,
    Y,
    p,
    t,
) where {FT}
    soil_boundary_fluxes!(bc,
                          Val(bc.components),
                          soil,
                          Y,
                          p,
                          t)
end

                          
