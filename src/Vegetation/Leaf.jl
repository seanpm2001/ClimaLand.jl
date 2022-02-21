module Leaf

using Photosynthesis, PlantHydraulics, StomataModels, WaterPhysics

using ClimaLSM
using ClimaCore
using UnPack
using DocStringExtensions
import ClimaCore: Fields

using ClimaLSM.Domains: AbstractVegetationDomain, RootDomain
import ClimaLSM:
    AbstractModel,
    initialize_prognostic,
    make_update_aux,
    make_rhs,
    make_ode_function,
    prognostic_vars,
    auxiliary_vars,
    initialize,
    initialize_auxiliary

export  LandLeafModel, LeafModel
abstract type AbstractLeafConductance{FT <: AbstractFloat} end

# Define our own Land model here:
Base.@kwdef struct LandLeafModel{FT}
    "Photosynthesis model"
    psm::Photosynthesis.AbstractPhotoModelParaSet{FT}
    "Stomatal Model"
    sm::StomataModels.AbstractStomatalModel{FT}
    "Leaf Hydraulic System"
    lhs::PlantHydraulics.AbstractHydraulicOrgan{FT}
end

struct LeafModel{FT,PS,CL,AL} <: AbstractModel{FT}
    param_set::PS
    "Canopy Layer"
    clayer::CL
    "Air Layer"
    environ::AL
    "The model name, which must be :leaf"
    model_name::Symbol
end

function LeafModel{FT}(;
    param_set::LandLeafModel{FT},
    clayer::CanopyLayer{FT},
    environ::AirLayer{FT},
) where {FT}
    args = (param_set, clayer, environ)
    return LeafModel{FT, typeof.(args)...}(args..., :leaf)
end

"""
    prognostic_vars(model::LeafModel)

A function which returns the names of the prognostic 
variables of the `LeafModel`.
"""
prognostic_vars(model::LeafModel) = (:g_sw,)
auxiliary_vars(model::LeafModel)  = (:An, :NPQ, :Ci, :APAR)

"""
    make_rhs(model::LeafModel)

A function which creates the rhs! function for the LeafModel.

The rhs! function must comply with a rhs function of OrdinaryDiffEq.jl.
"""
function make_rhs(model::LeafModel)
    function rhs!(dY, Y, p, t)
        @unpack psm, sm, lhs = model.param_set
        @unpack clayer, environ = model
        # Now set the stomatal conductance to the prognostic variable:
        #@show Y.leaf.g_sw[1]
        gas_exchange!(psm, clayer, environ, GswDrive());
 
        # set derivative of prognostic variable (brute force here):
        ∂gs∂t = dgsdt(psm,environ, clayer,environ,sm,1.0);
        #@show t, p.leaf.APAR, ∂gs∂t
        @. dY.leaf.g_sw .= ∂gs∂t[1] 
    end
    return rhs!
end

"""
    make_update_aux(model::RichardsModel)

An extension of the function `make_update_aux`, for the Richardson-
Richards equation. 

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function make_update_aux(model::LeafModel)
    function update_aux!(p, Y, t)
        @unpack clayer = model
        # clayer.g_sw  .= Y.leaf.g_sw[1]
        #@show clayer.g_sw
        p.leaf.An   = clayer.ps.An
        p.leaf.NPQ  = clayer.ps.NPQ
        p.leaf.Ci   = clayer.ps.p_i
        p.leaf.APAR = get_APAR(t)
        clayer.APAR  .= p.leaf.APAR
        clayer.g_sw .= Y.leaf.g_sw[1]
        #@show p.leaf.APAR
        #clayer.ps.APAR = p.leaf.APAR[1];  
    end
    return update_aux!
end

# Dummy timeseries for APAR:
function get_APAR(t::FT) where FT
    # Jump every 30min
    index = FT(floor(t/1000/3))
    return index * FT(500) + FT(20)
end

function dgsdt(
    psm,
    environ,
    clayer::CanopyLayer{FT},
    envir::AirLayer{FT},
    sm::EmpiricalStomatalModel{FT},
    β::FT
) where {FT<:AbstractFloat}
# unpack values
@unpack g_sw, n_leaf = clayer;

    # update g_sw
    dgsdt = similar(g_sw);
    #@show "Haeh?", n_leaf
    for iLF in 1:n_leaf
        gsw_ss_1 = max(sm.g0[iLF],stomatal_conductance(sm, clayer, envir, β, iLF));
        g_sw[iLF] = gsw_ss_1
        gas_exchange!(psm, clayer, environ, GswDrive());
        g_sw[iLF] = gsw_ss_1
        gsw_ss_1 = max(sm.g0[iLF],stomatal_conductance(sm, clayer, envir, β, iLF));
        #@show gsw_ss_1, iLF,sm.g0[iLF]
        dgsdt[iLF] = (gsw_ss_1 - g_sw[iLF]) / clayer.τ_esm;
        #@show dgsdt[iLF],g_sw[iLF]
    end

    return dgsdt
end

end