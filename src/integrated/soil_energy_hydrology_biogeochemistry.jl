export LandSoilBiogeochemistry, PrognosticSoilforBiogeochemistry

"""
    struct LandSoilBiogeochemistry{
        FT,
        SEH <: Soil.EnergyHydrology{FT},
        SB <: Soil.Biogeochemistry.SoilCO2Model{FT},
    } <: AbstractLandModel{FT}

A concrete type of land model used for simulating systems with a
soil energy, hydrology, and biogeochemistry component.
$(DocStringExtensions.FIELDS)"""
struct LandSoilBiogeochemistry{
    FT,
    SEH <: Soil.Soil.EnergyHydrology{FT},
    SB <: Soil.Biogeochemistry.SoilCO2Model{FT},
} <: AbstractLandModel{FT}
    "The soil model"
    soil::SEH
    "The biochemistry model"
    soilco2::SB
end

"""
    LandSoilBiogeochemistry{FT}(;
        soil_args::NamedTuple = (;),
        biogeochemistry_args::NamedTuple = (;),
    ) where {FT}
A constructor for the `LandSoilBiogeochemistry` model, which takes in
the required arguments for each component, constructs those models,
and constructs the `LandSoilBiogeochemistry` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.

Additional arguments, like parameters and driving atmospheric data, can be passed
in as needed.
"""
function LandSoilBiogeochemistry{FT}(;
    soil_args::NamedTuple = (;),
    soilco2_args::NamedTuple = (;),
) where {FT}

    soil = Soil.EnergyHydrology{FT}(;
        soil_args..., # soil_args must have sources, boundary_conditions, domain, parameters
    )

    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(; soilco2_args...)
    if !(soilco2_args.drivers.soil isa PrognosticSoilforBiogeochemistry)
        throw(
            AssertionError(
                "To run with a prognostic soil energy and hydrology model, the soil driver must be of type PrognosticSoilforBiogeochemistry.",
            ),
        )
    end
    if soil_args.parameters.ν != soilco2_args.drivers.ν
        throw(
            AssertionError(
                "The provided porosity must be the same in both models.",
            ),
        )
    end

    args = (soil, soilco2)
    return LandSoilBiogeochemistry{FT, typeof.(args)...}(args...)
end

struct PrognosticSoilforBiogeochemistry{
    FT,
    F1 <: Function,
    F <: Union{AbstractFloat, ClimaCore.Fields.Field},
} <: Soil.Biogeochemistry.AbstractSoilDriver
    "Carbon content of soil organic matter, of the form f(z::FT, t) where FT <: AbstractFloat"
    soil_organic_carbon::F3
    "Soil porosity (m³ m⁻³)"
    ν::F
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)"
    θ_a100::F
    "Absolute value of the slope of the line relating log(ψ) versus log(θ) (unitless)"
    b::F
end

"""
    soil_temperature(driver::PrognosticSoil, p, _...)
Returns the prognostic soil temperature.
"""
function soil_temperature(driver::PrognosticSoilforBiogeochemistry, p, _...)
    return p.soil.T
end

"""
    soil_moisture(driver::PrognosticSoil, p, _...)

Returns the volumetric liquid fraction, computed by the soil
model from the prognostic liquid and ice fractions.
"""
function soil_moisture(driver::PrognosticSoilforBiogeochemistry, p, _...)
    return p.soil.θ_l
end

"""
    soil_som_C(driver::PrognosticSoilforBiogeochemistry, p, Y, t, z)

Returns the carbon soil organic matter (SOM) at location (z) and time (t) for the prognostic soil
case, which still has a prescribed soil carbon.
"""
function soil_SOM_C(driver::PrognosticSoilforBiogeochemistry, p, Y, t, z)
    return driver.soil_organic_carbon.(z, t)
end

function ClimaLand.get_drivers(model::LandSoilBiogeochemistry)
    bc = model.soil.boundary_conditions.top
    if typeof(bc) <: AtmosDrivenFluxBC{
        <:PrescribedAtmosphere,
        <:AbstractRadiativeDrivers,
        <:Soil.AbstractRunoffModel,
    }
        return (bc.atmos, bc.radiation)
    else
        return (model.soilco2.driver.atmos, nothing)
    end
end
