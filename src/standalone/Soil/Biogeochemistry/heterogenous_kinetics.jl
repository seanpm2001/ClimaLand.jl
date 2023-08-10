struct HeterogenousDecayModel{FT, PS, D, BC, S, DT} <:
    AbstractSoilBiogeochemistryModel{FT}
 "the parameter set"
 parameters::PS
 "the biogeochemistry domain, using ClimaCore.Domains"
 domain::D
 "A tuple of sources, each of type AbstractSource"
 sources::S
 " Drivers"
 drivers::DT
end

ClimaLSM.name(model::HeterogenousDecayModel) = :biogeoC 
ClimaLSM.domain_name(model::HeterogenousDecayModel) = :surface

ClimaLSM.prognostic_vars(::HeterogenousDecayModel) = (:C_top30cm,)
ClimaLSM.prognostic_types(::HeterogenousDecayModel{FT}) where {FT} = (FT,)
# ma = moving average; can rename 
# k is decay rate distribution, defined by mean and variance of a log normal distribution
ClimaLSM.auxiliary_vars(::HeterogenousDecayModel) = (:T_soil_ma, :θ_soil_ma, :k)
ClimaLSM.auxiliary_types(::HeterogenousDecayModel{FT}) where {FT} = (FT, FT, FT)

function ClimaLSM.make_compute_exp_tendency(model::HeterogenousDecayModel)
    function compute_exp_tendency!(dY, Y, p, t)
        @. dY.biogeoC.C_top30cm = - p.biogeoC.k * Y.biogeoC.C_top30cm #+ p.biogeoC.source;
    end
    return compute_exp_tendency!
end

function ClimaLSM.make_update_aux(model::HeterogenousDecayModel)
    function update_aux!(p, Y, t)
        # update moving average quantities using model.drivers, model.parameters and t
        # also access the soil T and moisture from the model.drivers
        # write the update rule
        # Depth averaged; at the current time t
        T_soil = soil_temperature(model.driver.met, p, Y, t, nothing)
        θ_soil = soil_moisture(model.driver.met, p, Y, t, nothing)
        T = p.biogeoC.T_soil_ma
        θ = p.biogeoC.θ_soil_ma 
        # Compute depth average of T_soil

        # Update moving average
        T .= moving_average(T, T_soil)
        # Now use these to compute the distribution parameters
        p.biogeoC.mean_k .= distribution_mean.(T, θ, model.parameters)
        p.biogeoC.std_k .= distribution_std.(T, θ, model.parameters)

        # sample to update k
        # k ∼ P(k)
        N = size(k)
        p.biogeoC.k .= exp.(randn(N) .* p.biogeoC.std_k .+ p.biogeoC.mean_k)

        # update p.biogeoC.source with some parameterized function stored somewhere else
    end
    return compute_exp_tendency!
end




# When prescribed soil T and moisture, how are they supplied? depth resolved? 