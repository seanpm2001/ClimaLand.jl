using ClimaLSM.Domains: AbstractDomain
import ClimaLSM.Domains: make_function_space

### Example of component specific domain
"""
    AbstractVegetationDomain{FT} <: AbstractDomain{FT}

An abstract type for vegetation specific domains.
"""
abstract type AbstractVegetationDomain{FT} <: AbstractDomain{FT} end

"""
   RootDomain{FT} <: AbstractVegetationDomain{FT}

Domain for a single bulk plant with roots of varying depths. The user needs
to specify the depths of the root tips as wel as the heights of the
compartments to be modeled within the plant. The compartment heights
are expected to be sorted in ascending order.
"""
struct RootDomain{FT} <: AbstractVegetationDomain{FT}
    "The depth of the root tips, in meters"
    root_depths::Vector{FT}
    "The height of the stem, leaf compartments, in meters"
    compartment_heights::Vector{FT}
    "The surface height"
    z_sfc
end

function make_function_space(domain::RootDomain{FT}) where {FT}
    coord = ClimaCore.Geometry.ZPoint(domain.z_sfc)
    space = ClimaCore.Spaces.PointSpace(coord)
    return space, nothing
end
