module Soil
#=
    Soil

This module contains everything needed to run a soil model
in standalone mode.

The soil model is assumed to have a set of prognostic `Y` and
auxiliary `p` variables, which describe the state of the
soil system. The system is evolved in time by solving
equations of the form

```

\frac{\partial Y}{\partial t} = D(Y, t, p(Y, t; \ldots);\ldots),

```

i.e. partial (or ordinary) differential equations depending
on state `Y`,
auxiliary functions of the state `p`, and other parameters
represented by the ellipses. The operator `D` indicates a generic
nonlinear differential operator.  Not every model
requires auxilary variables, but these are useful for storing
quantities that are needed multiple times per right hand side
evaluation. An example would be the temperature, which is computed
from the prognostic variables of water content and internal
energy.

Currently, both the Richardson Richards Equation (RRE; hydrology alone)
and an integrated soil energy and hydrology model are supported.

Addition of additional versions of soil
models requires defining a model type (of super type
`AbstractSoilModel`), and extending the methods
imported by Models.jl, as needed, for computing the
right hand side functions of the ordinary differential equations
and the functions which update auxiliary variables whenever
the right hand side is evaluated.

This code base assumes that DifferentialEquations.jl
will be used for evolving the system in time,
and that the array-like objected being stepped forward
is a `ClimaCore.Fields.FieldVector`.

The `FieldVector` type is used in order to make use of
`ClimaCore` functionality when solving PDEs (`ClimaCore` handles
all of the spatial discretization and operator functionality)
 and for ease of handling multi-column models.

To simulate land surfaces with multiple components (vegetation,
soil, rivers, etc), the ClimaLand.jl package should be used.
That package will use the methods of this function for advancing
the system forward in time, extending methods as needed to account
for interactions between components.
=#

using ClimaLand
using DocStringExtensions
using LinearAlgebra
using ClimaCore
using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name, ⋅
import UnrolledUtilities
import ..Parameters as LP
import ClimaCore: Fields, Operators, Geometry, Spaces
using Thermodynamics

import ClimaLand.Domains: Column, HybridBox, SphericalShell
import ClimaLand:
    AbstractImExModel,
    make_update_aux,
    make_compute_exp_tendency,
    make_compute_imp_tendency,
    make_update_boundary_fluxes,
    make_update_jacobian,
    prognostic_vars,
    auxiliary_vars,
    prognostic_domain_names,
    auxiliary_domain_names,
    name,
    prognostic_types,
    auxiliary_types,
    AbstractSource,
    source!,
    heaviside,
    surface_temperature,
    surface_specific_humidity,
    surface_albedo,
    surface_emissivity,
    surface_height,
    surface_resistance,
    get_drivers
export RichardsModel,
    RichardsParameters,
    ImplicitEquationJacobian,
    EnergyHydrology,
    EnergyHydrologyParameters,
    AbstractSoilModel,
    AbstractSoilSource,
    PhaseChange



"""
    AbstractSoilSource{FT} <:  ClimaLand.AbstractSource{FT}

An abstract type for types of source terms for the soil equations.

In standalone mode, the only supported source type is freezing and
thawing. ClimaLand.jl creates additional sources to include as
necessary e.g. root extraction (not available in stand alone mode).
"""
abstract type AbstractSoilSource{FT} <: ClimaLand.AbstractSource{FT} end

"""
    AbstractSoilModel{FT} <: ClimaLand.AbstractImExModel{FT}

The abstract type for all soil models.

Currently, we only have plans to support a RichardsModel, simulating
the flow of liquid water through soil via the Richardson-Richards equation,
and a fully integrated soil heat and water model, with phase change.
"""
abstract type AbstractSoilModel{FT} <: ClimaLand.AbstractImExModel{FT} end

ClimaLand.name(::AbstractSoilModel) = :soil
"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Column, _...)
Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a column domain, there are no horizontal
contributions to the right hand side.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Column,
    _...,
) end

"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Union{HybridBox, SphericalShell},
                          lateral_flow::Val{false},
                          _...)
Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a 3D domain, for either the `RichardsModel` or the
`EnergyHydrology` model, if the `lateral_flow` flag is set to false,
there are no horizontal contributions to the right hand side.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
    lateral_flow::Val{false},
    _...,
) end

"""
    append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple
Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple = (srcs..., src)

"""
    append_source(src::Nothing , srcs::Tuple)::Tuple
Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::Nothing, srcs::Tuple)::Tuple = srcs

"""
    ImplicitEquationJacobian{M, S}

A struct containing the necessary information for constructing a tridiagonal
Jacobian matrix (`W`) for solving Richards equation, treating only the vertical
diffusion term implicitly.

`matrix` is a block matrix containing the tri-diagonal matrix `∂ϑres∂ϑ`
in the RichardsModel case.
`solver` is a diagonal solver for the RichardsModel case because our matrix is
block diagonal.

Note that the diagonal, upper diagonal, and lower diagonal entry values
are stored in this struct and updated in place.
"""
struct ImplicitEquationJacobian{M, S}
    "Jacobian matrix stored as a MatrixFields.FieldMatrix"
    matrix::M
    "Solver to use for solving the tridiagonal system"
    solver::S
end

"""
    ImplicitEquationJacobian(
        Y::ClimaCore.Fields.FieldVector;
)

Outer constructor for the ImplicitEquationJacobian Jacobian
matrix struct.

Initializes all variables to zeros.
"""
function ImplicitEquationJacobian(Y::ClimaCore.Fields.FieldVector)
    FT = eltype(Y)
    center_space = axes(Y.soil.ϑ_l)

    # Cosntruct a tridiagonal matrix that will be used as the Jacobian
    tridiag_type = MatrixFields.TridiagonalMatrixRow{FT}
    # Create a field containing a `TridiagonalMatrixRow` at each point
    tridiag_field = Fields.Field(tridiag_type, center_space)
    fill!(parent(tridiag_field), NaN)

    # Get all prognostic vars in soil, and create a tridiagonal matrix for each
    soil_varnames = MatrixFields.top_level_names(Y.soil)
    varnames = map(
        name -> MatrixFields.append_internal_name(@name(soil), name),
        soil_varnames,
    )
    matrix = MatrixFields.FieldMatrix(
        UnrolledUtilities.unrolled_map(
            x -> (x, x) => copy(tridiag_field),
            varnames,
        )...,
    )

    # Set up block diagonal solver for block Jacobian
    alg = MatrixFields.BlockDiagonalSolve()
    solver = MatrixFields.FieldMatrixSolver(alg, matrix, Y)

    return ImplicitEquationJacobian(matrix, solver)
end

Base.similar(w::ImplicitEquationJacobian) = w

function LinearAlgebra.ldiv!(
    x::Fields.FieldVector,
    A::ImplicitEquationJacobian,
    b::Fields.FieldVector,
)
    MatrixFields.field_matrix_solve!(A.solver, x, A.matrix, b)
end

include("./retention_models.jl")
include("./rre.jl")
include("./energy_hydrology.jl")
include("Runoff/Runoff.jl")
using .Runoff
include("./boundary_conditions.jl")
include("./soil_hydrology_parameterizations.jl")
include("./soil_heat_parameterizations.jl")
include("Biogeochemistry/Biogeochemistry.jl")
using .Biogeochemistry
end
