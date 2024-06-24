export PrognosticSoil, AbstractSoilDriver

"""
    AbstractSoilDriver

An abstract type for time varying and/or spatially varying 
soil properties (temperature, water content, organic carbon, etc)
which are either prescribed or prognostic.
"""
abstract type AbstractSoilDriver{FT <:AbstractFloat} <: ClimaLand.AbstractClimaLandDrivers{FT} end

struct PrognosticSoil{FT,
                      TVI <: Union{Nothing, AbstractTimeVaryingInput},
                      F <:Union{Nothing, AbstractFloat, ClimaCore.Fields.Field},
                      F2 <:Union{Nothing, AbstractFloat, ClimaCore.Fields.Field}
                      } <: AbstractSoilDriver{FT}
    "Carbon content of soil organic matter, of the form f(z::FT, t) where FT <: AbstractFloat"
    soil_organic_carbon::TVI
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)"
    θ_a100::F
    "Absolute value of the slope of the line relating log(ψ) versus log(θ) (unitless)"
    b::F
    α_PAR::F2
    α_NIR::F2
end
"""
    soil_temperature(driver::PrognosticSoil, p, Y)
Returns the prognostic soil temperature.
"""
function soil_temperature(driver::PrognosticSoil, p, Y)
    return p.soil.T
end

"""
    soil_moisture(driver::PrognosticSoil, p, Y, t, z)

Returns the volumetric liquid fraction, computed by the soil
model from the prognostic liquid and ice fractions.
"""
function soil_moisture(driver::PrognosticSoil, p, Y)
    return p.soil.θ_l
end

function Canopy.ground_albedo_PAR(soil_driver::PrognosticSoil, Y, p, t)
    return soil_driver.α_PAR
end

function Canopy.ground_albedo_NIR(soil_driver::PrognosticSoil, Y, p, t)
    return soil_driver.α_NIR
end
