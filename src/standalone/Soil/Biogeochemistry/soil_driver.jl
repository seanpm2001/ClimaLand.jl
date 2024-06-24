
"""
    AbstractSoilDriver{FT}

An abstract type for time varying and/or spatially varying 
soil properties (temperature, water content, organic carbon, etc)
which are either prescribed or prognostic.
"""
abstract type AbstractSoilDriver{F} <: AbstractClimaLandDrivers {FT} end

"""
    PrescribedSoil <: AbstractSoilDriver

A container which holds the prescribed functions for soil temperature, moisture,
and carbon.

This is meant for use when running the biogeochemistry model in standalone mode,
without a prognostic soil model.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedSoil{FT, F1 <: AbstractTimeVaryingInput, F2 <: AbstractTimeVaryingInput, F3 <: AbstractTimeVaryingInput, F <:Union{AbstractFloat, ClimaCore.Fields.Field}} <: AbstractSoilDriver
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

function ClimaLand.initialize_drivers(s::PrescribedSoil{FT}, coords) where {FT}
    keys = (:T,:θ_l, :soc, :ν, :θa_100, :b)
    types = (FT,)
    domain_names = (:subsurface,)
    model_name = :drivers
    # intialize_vars packages the variables as a named tuple,
    # as part of a named tuple with `model_name` as the key.
    # Here we just want the variable named tuple itself
    vars =
        ClimaLand.initialize_vars(keys, types, domain_names, coords, model_name)
    return vars.drivers
end

function PrescribedSoil{FT}(
    T_soil,
    θ_l_soil,
    soil_organic_carbon,
    ν::F,
    θa_100::F,
    b::F
) where {FT <: AbstractFloat}
    parametric_types = {FT, typeof(T_soil), typeof(θ_l_soil), typeof(soil_organic_carbon), typeof(ν)}
    return PrescribedSoil{parametric_types...}(T_soil, θ_l_soil, soil_organic_carbon, ν, θa_100, b)
end

"""
    soil_temperature(driver::PrescribedSoil, p, Y, t, z)

Returns the soil temperature for the prescribed
soil case.
"""
function soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return p.drivers.T_soil
end

"""
    soil_moisture(driver::PrescribedSoil, p, Y)

Returns the soil moisture for the prescribed
soil case.
"""
function soil_moisture(driver::PrescribedSoil, p, Y)
    return p.drivers.θ_l_soil
end

"""
    soil_som_C(driver::PrescribedSoil, p, Y)

Returns the carbon soil organic matter (SOM) for the prescribed
soil case.
"""
function soil_SOM_C(driver::PrescribedSoil, p, Y)
    return p.drivers.soil_organic_carbon
end

"""
    make_update_drivers(a::PrescribedSoil{FT}) where {FT}

Creates and returns a function which updates the driver variables
in the case of a PrescribedSoil.
"""
function make_update_drivers(d::PrescribedSoil{FT}) where {FT}
    function update_drivers!(p, t)
        evaluate!(p.drivers.T_soil, d.T_soil, t)
        evaluate!(p.drivers.θ_l_soil, d.θ_l_soil, t)
        evaluate!(p.drivers.soil_organize_carbon, d.soil_organic_carbon, t)
    end
    return update_drivers!
end
