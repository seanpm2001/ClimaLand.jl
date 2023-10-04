using Flux, LinearAlgebra
using DataFrames, Dates

"""
    make_model(nfeatures, n, z_idx, p_idx; in_scale, dtype)

Create the neural network to be trained, with initial scaling weights.

# Arguments
- `nfeatures::Int`: indicates number of features.
- `n::Int`: the value of the hyperparameter n.
- `z_idx::Int`: The index of the data vectors pertaining to the depth (z) values.
- `p_idx::Int`: The index of the data vectors pertaining to the precipitation values.
-  `in_scale::Vector{<:Real}`: Optional scaling constants for each input feature.
- `dtype::Type`: Sets type of output model. Default is Float32.
"""
function make_model(nfeatures::Int, n::Int, z_idx::Int, p_idx::Int;
    in_scale::Union{Vector{<:Real}, Nothing} = nothing,
    dtype::Type = Float32)
    in_scales = (isnothing(in_scale)) ? Matrix{dtype}(diagm(ones(nfeatures))) : Matrix{dtype}(diagm((1.0 ./ in_scale)))
    get_relus = Matrix{dtype}([1 0 0; 0 1 0; 1 0 -1])
    get_min = Matrix{dtype}([1 1 -1; 0 1 0])
    get_max = Matrix{dtype}([1 -1])
    model = Chain(
        pred=SkipConnection(
            Chain(
                scale = Dense(in_scales, false, identity),
                l1 = Dense(nfeatures, n * nfeatures, relu),
                l2 = Dense(n * nfeatures, nfeatures, elu),
                l3 = Dense(nfeatures, 1)
            ),
            vcat #returns [predicted value, input...]
        ),
        get_boundaries=Parallel(vcat,
            up_bound=x -> relu.(x[1, :])' .* (x[p_idx+1, :] .> 0)', # = upper = relu(upper)
            low_bound=x -> x[z_idx+1, :]', # = z = relu(z)
            output_pos=x -> x[1, :]'
        ),
        final_scale=Dense(Matrix{dtype}(diagm([1.0, 1.0, 1.0])), false, identity),
        #output_no_thresh = x -> x[3, :]'
        apply_relus=Dense(get_relus, false, relu),  #returns relu(upper) = upper, relu(z) = z, relu(upper - pred_z)
        apply_upper=Dense(get_min, false, relu), #returns relu(min(pred_z, upper)+z), relu(z) = z
        apply_lower=Dense(get_max, false, identity), #returns relu(min(pred_z, upper)+z) - relu(z) = max(min(pred_Z, upper), lower)
    )
    return model
end


"""
    get_model_ps(model)

Return the trainable weights for the developed neural model.

# Arguments
- `model::Chain`: the neural model to be used.
"""
function get_model_ps(model::Chain)
    return Flux.params(Flux.params(model[:pred])[2:7])
end

"""
    settimescale!(model, dt; dtype)

Set the timescale parameter for model usage.

# Arguments
- `model::Chain`: the neural model to be used.
- `dt::Real`: the number of seconds per timestep for usage.
- `dtype::Type`: Sets type, consistent with neural model. Default is Float32.
"""
function settimescale!(model, dt::Real; dtype::Type = Float32)
    model[:final_scale].weight[2, 2] = dtype(1.0 / dt)
end

"""
    setoutscale!(model, scale; dtype)

Set the physical scaling parameter for model usage (i.e. rectifying scaling done on model input).

# Arguments
- `model::Chain`: the neural model to be used.
- `scale::Real`: the scaling parameter to return data to applicable units.
- `dtype::Type`: Sets type, consistent with neural model. Default is Float32.
"""
function setoutscale!(model, scale::Real; dtype::Type = Float32)
    model[:final_scale].weight[3, 3] = dtype(scale)
end

"""
    LRmodel(data, vars, target; dtype, scale_const)

Create a linear regression model on a data frame for comparison to neural model.
Returns the coefficients for the model.

# Arguments
- `data::DataFrame`: The data set to be utilized.
- `vars::Vector{Symbol}`: The input variables to be used in the model creation.
-  `target::Symbol`: The target variable to be used in the model ceation.
- `dtype::Type`: Sets type, consistent with neural model. Default is Float32.
- `scale_const`: Optional scaling constant for model output. Default is 1.0.
"""
function LRmodel(data::DataFrame, vars::Vector{Symbol}, target::Symbol; dtype::Type=Float32, scale_const = 1.0)
    X = Matrix{dtype}(select(data, vars))
    y = Vector{dtype}(data[!, target]) ./ dtype(scale_const)
    constants = [X ones(nrow(data))] \ y
    return constants .* scale_const
end

"""
    LRmodel(x_train, y_train; dtype, scale_const)

Create a linear regression model on a training matrix for comparison to neural model.
Returns the coefficients for the model.
**Note: using the same matrices input to the neural model will require a transpose of x_train, y_train

# Arguments
- `x_train::Matrix`: The input to be utilized.
- `y_train::Vector`: The output data to be utilized.
- `dtype::Type`: Sets type, consistent with neural model. Default is Float32.
- `scale_const`: Optional scaling constant for model output. Default is 1.0.
"""
function LRmodel(x_train::Matrix, y_train::Vector; dtype::Type = Float32, scale_const = 1.0)
    #using x_train from neural input will require a transpose of x_train, y_train
    return Vector{dtype}(([x_train ones(size(xtrain)[1])] \ y_train) .* scale_const)
end

"""
    eval(model, input)

Evaluate a created model on a given input vector.

# Arguments
- `model::Chain`: A neural model to be used for prediction.
- `input`: The input data used to generate a prediction.
"""
function eval(model::Chain, input)
    return model(input)
end

"""
    eval(model, input)

Evaluate a created model on a given input vector.

# Arguments
- `model::Vector{<:Real}`: Linear regression coefficients used for prediction.
- `input`: The input data used to generate a prediction.
"""
function eval(model::Vector{<:Real}, input)
    #requires input matrix to be the same orientation as that for the neural model
    return model[1:end-1]' * input .+ model[end]
end


"""
    make_timeseries(model, timeseries, dt; predictvar, timevar, inputvars, dtype, hole_thresh)

Generate a predicted timeseries given forcing data and the timestep present in that data (holes acceptable).

# Arguments
- `model`: The model used for forecasting (can be any model with a defined "eval" call).
- `timeseries::DataFrame`: The input data frame used to generate predictions, including a time variable.
- `dt::Period`: The unit timestep present in the dataframe (i.e. daily dataframe, dt = Day(1) or Second(86400)).
- `predictvar::Symbol`: The variable to predict from the timeseries. Default is :z.
- `timevar::Symbol`: The variable giving the time of each forcing. Default is :date.
- `inputvars::Vector{Symbol}`: The variables (in order), to extract from the data to use for predictions.
Default is [:z, :SWE, :rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg, :dprecipdt_snow] like the paper.
- `dtype::Type`: The data type required for input to the model. Default is Float32.
- `hole_thresh::Int`: The acceptable number of "holes" in the timeseries for the model to skip over. Default is 30.
"""
function make_timeseries(model, timeseries::DataFrame, dt::Period;
    predictvar::Symbol = :z,
    timevar::Symbol = :date,
    inputvars::Vector{Symbol} = [:z, :SWE, :rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg, :dprecipdt_snow],
    dtype::Type=Float32,
    hole_thresh::Int=30)
    forcings = Matrix{dtype}(select(timeseries, inputvars))
    pred_idx = findfirst(inputvars .== predictvar)
    check_dates = (timeseries[2:end, timevar] - timeseries[1:end-1, timevar]) ./ dt
    pred_vals = zeros(nrow(timeseries) - 1)
    pred_series = zeros(nrow(timeseries))
    pred_series[1] = timeseries[1, predictvar]
    countresets = 0
    for j in 2:length(pred_series)
        input = forcings[j-1, :]
        input[pred_idx] = pred_series[j-1]
        pred = eval(model, input)[1]
        pred_vals[j-1] = pred
        nperiods = check_dates[j-1]
        new_val = pred_series[j-1] + nperiods * Dates.value(Second(dt)) * pred
        pred_series[j] = (nperiods <= hole_thresh) ? max(0.0, new_val) : timeseries[j, predictvar]  #the "max" is only in the case of holes for a non-negative system and does not generalize
        #pred_series[j] = (nperiods <= hole_thresh) ? new_val : true_series[j]  #for showing no thresholds
        if nperiods > hole_thresh
            countresets += 1
        end
    end
    return pred_series, pred_vals, countresets
end