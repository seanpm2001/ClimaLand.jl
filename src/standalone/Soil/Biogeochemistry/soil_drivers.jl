"""
    AbstractSoilDriver

An abstract type for drivers of soil CO2 production and diffusion.
These are soil temperature, soil moisture, soil properties,
root carbon, soil organic matter and microbe carbon, and atmospheric pressure.
Soil temperature and moisture, as well as soc, vary in space (horizontally and vertically) and time.
Atmospheric pressure vary in time (defined at the surface only, not with depth).
"""
abstract type AbstractSoilDriver end

"""
    SoilCO2Drivers

A container which passes in the soil drivers to the biogeochemistry
model. These drivers are either of type Prescribed (for standalone mode)
or Prognostic (for running with a prognostic model for soil temp and moisture).

$(DocStringExtensions.FIELDS)
"""
struct SoilCO2Drivers{
    FT,
    S <: AbstractSoilDriver,
    ATM <: PrescribedAtmosphere{FT},
}
    "Soil temperature and moisture drivers - Prescribed or Prognostic"
    soil::S
    "Prescribed atmospheric variables"
    atmos::ATM
end

"""
    PrescribedSoil <: AbstractSoilDriver

A container which holds the prescribed functions for soil temperature
and moisture.

This is meant for use when running the biogeochemistry model in standalone mode,
without a prognostic soil model.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedSoil{FT, F1 <: Function, F2 <: Function, F3 <: Function, F <:Union{AbstractFloat, ClimaCore.Fields.Field}} <: AbstractSoilDriver
    "The temperature of the soil, of the form f(z::FT,t) where FT <: AbstractFloat"
    temperature::F1
    "Soil moisture, of the form f(z::FT,t) FT <: AbstractFloat"
    volumetric_liquid_fraction::F2
    "Carbon content of soil organic matter, of the form f(z::FT, t) where FT <: AbstractFloat"
    soil_organic_carbon::F3
    "Soil porosity (m³ m⁻³)"
    ν::F
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)"
    θ_a100::F
    "Absolute value of the slope of the line relating log(ψ) versus log(θ) (unitless)"
    b::F

end

function PrescribedSoil{FT}(
    temperature::Function,
    volumetric_liquid_fraction::Function,
    soil_organic_carbon::Function,
    ν::F,
    θ_a100::F,
    b::F
) where {FT <: AbstractFloat, F}
    return PrescribedMet{
        FT,
        typeof(temperature),
        typeof(volumetric_liquid_fraction),
        typeof(soil_organic_carbon),
        F}(
            temperature,
            volumetric_liquid_fraction,
            soil_organic_carbon,
            ν
            θ_a100,
            b
        )
end


"""
    soil_temperature(driver::PrescribedSoil, p, Y, t, z)

Returns the soil temperature at location (z) and time (t) for the prescribed
soil case.
"""
function soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.temperature.(z, t)
end

"""
    soil_moisture(driver::PrescribedSoil, p, Y, t, z)

Returns the soil moisture at location (z) and time (t) for the prescribed
soil case.
"""
function soil_moisture(driver::PrescribedSoil, p, Y, t, z)
    return driver.volumetric_liquid_fraction.(z, t)
end

"""
    soil_som_C(driver::PrescribedSoil, p, Y, t, z)

Returns the carbon soil organic matter (SOM) at location (z) and time (t) for the prescribed
soil case.
"""
function soil_SOM_C(driver::PrescribedSoil, p, Y, t, z)
    return driver.soil_organic_carbon.(z, t)
end
