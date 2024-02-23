abstract type AbstractBiomassModel{FT <: AbstractFloat} end

"""
   PrescribedBiomassModel{FT}{F <: AbstractTimeVaryingInput,
                          FS <: Union{AbstractFloat, ClimaCore.Fields.Field},
                          }

A struct containing the prescribed area indices of the plants: LAI varies in time, 
while SAI and RAI are fixed.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedBiomassModel{FT <: AbstractFloat,
    F <: AbstractTimeVaryingInput,
    FS <: Union{AbstractFloat, ClimaCore.Fields.Field}
} <: AbstractBiomassModel{FT}
    "A function of simulation time `t` giving the leaf area index (LAI; m2/m2)"
    LAIfunction::F
    "The stem area index (SAI; m2/m2)"
    SAI::FS
    "The root area index (RAI; m2/m2)"
    RAI::FS
end

function PrescribedBiomassModel(
    LAIfunction::AbstractTimeVaryingInput,
    SAI::FS,
    RAI::FS,
    ) where {FS}
    FT = eltype(SAI)
    PrescribedBiomassModel{FT, typeof(LAIfunction), FS}(LAIfunction, SAI, RAI)
end

ClimaLand.name(model::AbstractBiomassModel) = :biomass
ClimaLand.auxiliary_vars(model::PrescribedBiomassModel) = (:LAI,)
ClimaLand.auxiliary_types(model::PrescribedBiomassModel{FT}) where {FT} =
    (FT,)
ClimaLand.auxiliary_domain_names(::PrescribedBiomassModel) =
    (:surface,)

function update_area_indices!(model::PrescribedBiomassModel, canopy, p, t)
    LAI = p.canopy.biomass.LAI
    evaluate!(LAI, model.LAIfunction, t)
    lai_consistency_check.(canopy.hydraulics.n_stem, canopy.hydraulics.n_leaf, LAI,
                           model.SAI,
                           model.RAI)
end

function area_index(name, model::PrescribedBioMassModel, p)
    if name == :leaf
        p.canopy.biomass.LAI
    elseif name == :stem
        model.SAI
    else
        model.RAI
    end
end

"""
    lai_consistency_check(
        n_stem::Int64,
        n_leaf::Int64,
        LAI::FT,
        SAI::FT,
        RAI::FT
    ) where {FT}

Carries out consistency checks using the area indices supplied and the number of
stem elements being modeled.

Note that it is possible to have a plant with no stem compartments
but with leaf compartments, and that a plant must have leaf compartments
(even if LAI = 0).

Specifically, this checks that:
1. n_leaf > 0
2. if LAI is nonzero or SAI is nonzero, RAI must be nonzero.
3. if SAI > 0, n_stem must be > 0 for consistency. If SAI == 0, n_stem must
be zero.
"""
function lai_consistency_check(
    n_stem::Int64,
    n_leaf::Int64,
    LAI::FT,
    SAI::FT,
    RAI::FT
) where {FT}
    @assert n_leaf > 0
    if LAI > eps(FT) || SAI > eps(FT)
        @assert RAI > eps(FT)
    end
    # If the SAI is zero, there should be no stem compartment
    if SAI < eps(FT)
        @assert n_stem == FT(0)
    else
        # if SAI is > 0, n_stem should be > 0 for consistency
        @assert n_stem > 0
    end

end    
