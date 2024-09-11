module NeuralSnow

using Flux, Dates
using ClimaCore
using ClimaLand
using ClimaLand.Snow: AbstractDensityModel, SnowModel, SnowParameters
import ClimaLand.Snow: dzdt, snow_density
import ClimaLand.Parameters as LP
using Thermodynamics

using CSV, HTTP, Flux, cuDNN, StatsBase, BSON
ModelTools = Base.get_extension(ClimaLand, :NeuralSnowExt).ModelTools;

export snow_density, dzdt, NeuralDepthModel

"""
    NeuralDepthModel{FT <: AbstractFloat} <: AbstractDensityModel{FT}

Establishes the density parameterization where snow density
is calculated from a neural network determining the rate of change of snow depth, dzdt.
"""
struct NeuralDepthModel{FT} <: AbstractDensityModel{FT}
    z_model::Flux.Chain
    Δt::FT #make constant? Or leave adaptable for different algos?
    type::Symbol
    sol_avg::Vector{ClimaCore.Fields.Field} #if EMA is the best type, then decision on whether to leave this or delete and use dY. var instead
    wind_avg::Vector{ClimaCore.Fields.Field}
    temp_avg::Vector{ClimaCore.Fields.Field}
    hum_avg::Vector{ClimaCore.Fields.Field}
    precip_avg::Vector{ClimaCore.Fields.Field}
    α::FT
end

"""
    get_znetwork(; download_link = "https://caltech.box.com/shared/static/ay7cv0rhuiytrqbongpeq2y7m3cimhm4.bson")

A default function for returning the snow-depth neural network from the paper.
"""
#decision on whether to include this function determines dependency of BSON for the neural extension
function get_znetwork(; download_link = "https://caltech.box.com/shared/static/ay7cv0rhuiytrqbongpeq2y7m3cimhm4.bson")
    z_idx = 1
    p_idx = 7
    nfeatures = 7
    zmodel = ModelTools.make_model(nfeatures, 4, z_idx, p_idx)
    zmodel_state = BSON.load(IOBuffer(HTTP.get(download_link).body))[:zstate]
    Flux.loadmodel!(zmodel, zmodel_state)
    ModelTools.settimescale!(zmodel, 86400.0)
    return zmodel

    #give ability to also load from local file too?
    #BSON.@load "../SnowResearch/cleancode/testhrmodel.bson" zstate
    #Flux.loadmodel!(zmodel, zstate)
end

"""
    NeuralDepthModel{FT}(type::Symbol, Δt::FT; model::Flux.Chain = get_znetwork(), α = nothing)

An outer constructor for the `NeuralDepthModel` density parameterization for usage in a snow model.
Can choose to use custom models via the `model` argument, or a custom weight with `α` if using `type` `:DEMA` or `:EMA`.
the `type` argument dictates the formulation of inputs to the network (moving averages or instantaneous values of inputs).
"""
function NeuralDepthModel(Δt::FT; type::Symbol = :EMA, model::Flux.Chain = get_znetwork(), α = nothing) where {FT}
    usemodel = model
    if FT == Float32
        usemodel = Flux.f32(model)
    elseif FT == Float64
        usemodel = Flux.f64(model)
    end
    if !in(type, [:instantaneous, :daily, :EMA, :DEMA])
        error("Requested type of NeuralDepthModel is not implemented. Existing models are ':instantaneous', ':daily', ':EMA', :'DEMA'")
    end
    num_idx = 1
    weight = FT(0)
    if type == :daily
        #set this based off of Δt, timestepping algo, etc.
        num_idx = 24*4
    end
    if type == :DEMA
        num_idx = 2
    end
    if type == :instantaneous
        weight = FT(1)
    end
    if type in [:EMA, :DEMA]
        #auto-calculate value of weight for rough daily average or according to Δt and model timestep, etc:
        weight = FT(2/97.0)
    end
    if !isnothing(α)
        weight = FT(α)
    end
    sol_avg = Vector{ClimaCore.Fields.Field}(undef, num_idx)
    wind_avg = Vector{ClimaCore.Fields.Field}(undef, num_idx)
    temp_avg = Vector{ClimaCore.Fields.Field}(undef, num_idx)
    hum_avg = Vector{ClimaCore.Fields.Field}(undef, num_idx)
    precip_avg = Vector{ClimaCore.Fields.Field}(undef, num_idx)
    return NeuralDepthModel{FT}(usemodel, Δt, type, sol_avg, wind_avg, temp_avg, hum_avg, precip_avg, weight)
end

"""
    fieldavg(vec::Vector)

Helper function for evaluating the neural network in a pointwise manner over a `ClimaCore.Field`
and returning the output in a broadcastable way.

"""
function eval_nn(model, z::FT, swe::FT, P::FT, T::FT, R::FT, qrel::FT, u::FT)::FT where {FT}
    input = FT.([z, swe, qrel, R, u, T, P])
    return model(input)[1]
end

"""
    fieldavg(vec::Vector)

Helper function for finding the first unassigned index in a vector of `ClimaCore.Field`s.

"""
function first_unassigned(vec::Vector)
    idxs = collect(1:length(vec))
    idx = findfirst((!).(isassigned.(Ref(vec), idxs)))
    return isnothing(idx) ? 0 : idx
end


"""
    fieldavg(vec::Vector)

Helper function for finding the average value of a vector of `ClimaCore.Field`s.

"""
function fieldavg(vec::Vector)
    avg = deepcopy(vec[1]) #need this to avoid editing the actual fields
    n = first_unassigned(vec)
    n = isnothing(n) ? length(vec) : n-1
    if n == 1
        return avg
    else
        for elem in vec[2:n]
            avg .+= elem
        end
    end
    return avg ./ n
end


"""
    dzdt(density::NeuralDepthModel, model::SnowModel{FT}, Y, p, t) where {FT}

Returns the change in snow depth (rate) given the current model state and the `NeuralDepthModel`
density paramterization. The calculation will differ slightly depending on the `type` of `NeuralDepthModel`
passed.

"""
function dzdt(density::NeuralDepthModel, model::SnowModel{FT}, Y, p, t) where {FT}
    #get inputs
    z = Y.snow.Z
    swe = Y.snow.S
    dprecipdt_snow = abs.(p.drivers.P_snow) #need to change downward direction to scalar
    air_temp = p.drivers.T .- model.parameters.earth_param_set.T_freeze
    sol_rad = p.drivers.SW_d
    rel_hum = Thermodynamics.relative_humidity.(Ref(model.atmos.thermo_params), p.drivers.thermal_state)
    wind_speed = p.drivers.u

    #update stored input fields, depending on type:
    local idx    
    for (field, val) in [(density.sol_avg, sol_rad),
        (density.hum_avg, rel_hum),
        (density.temp_avg, air_temp),
        (density.wind_avg, wind_speed),
        (density.precip_avg, dprecipdt_snow)]
        
        idx = first_unassigned(field)
        if idx == 0
            if density.type in [:instantaneous, :EMA, :DEMA]
                field[1] = (density.α .* val) .+ ((FT(1) - density.α) .* field[1])
                if density.type == :DEMA
                    field[2] = (density.α .* (FT(2) .* val .- field[1])) .+ ((FT(1) - density.α) .* field[2])
                end
            elseif density.type == :daily
                popfirst!(field)
                push!(field, val)
            end
            idx += length(field)
        elseif (idx == 2) & (density.type == :DEMA)
            temp = field[1]
            field[1] = (density.α .* val) .+ ((FT(1) - density.α) .* field[1])
            field[2] = (density.α .* (FT(2) .* val .- field[1])) .+ ((FT(1) - density.α) .* temp)
        else
            field[idx] = val
        end
    end

    sol_rad_avg = (density.type == :daily) ? fieldavg(density.sol_avg) : density.sol_avg[idx]
    rel_hum_avg = (density.type == :daily) ? fieldavg(density.hum_avg) : density.hum_avg[idx]
    air_temp_avg = (density.type == :daily) ? fieldavg(density.temp_avg) : density.temp_avg[idx]
    wind_speed_avg = (density.type == :daily) ? fieldavg(density.wind_avg) : density.wind_avg[idx]
    precip_snow = (density.type == :daily) ? dprecipdt_snow : density.precip_avg[idx] #makes better results

    #eval neural network
    dzdt = eval_nn.(Ref(density.z_model), z, swe, precip_snow, air_temp_avg, sol_rad_avg, rel_hum_avg, wind_speed_avg)
    return dzdt
end

"""
    snow_density(density::NeuralDepthModel, SWE::FT, z::FT, parameters::SnowParameters{FT}) where {FT}

Returns the snow density given the current model state and the `NeuralDepthModel`
density paramterization.

"""
function snow_density(density::NeuralDepthModel, SWE::FT, z::FT, parameters::SnowParameters{FT})::FT where {FT}
    ρ_l = FT(LP.ρ_cloud_liq(parameters.earth_param_set))
    if SWE == FT(0) #if there is no snowpack, aka avoid NaN
        return FT(ρ_l)
    end
    ρ_new = SWE / z * ρ_l
    return ρ_new
end


end