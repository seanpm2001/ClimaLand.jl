using Flux
using DataFrames, Dates, StatsBase
using Flux.Data: DataLoader
using Flux: loadparams!

"""
    show_mem_usage(ret)

Print the current RAM usage (percentage).

# Arguments
- `ret::Bool`: indicates whether to return the percentage as a number.
"""
function show_mem_usage(ret::Bool = false)
    used_mem = ((Sys.total_memory() / 2^20) - (Sys.free_memory() / 2^20)) / 1000.0
    perc = 100.0 * (1.0 - Sys.free_memory()/Sys.total_memory())
    print("RAM Usage: ", round(used_mem, digits = 3), " GB (", round(perc, digits = 2), "%)\n")
    if ret
        return used_mem
    end
end

"""
    make_data(data, input_vars, target, target_scale; testidx, dtype)

Create training and testing matrices for the models from prepped input data.

# Arguments
- `data::DataFrame`: the data to be used for training and validation data.
- `input_vars::Vector{Symbol}`: The set of features to be extracted as input data
- `target::Symbol`: The feature to be extracted as the target variable
- `target_scale::Real`: The scaling to apply to the output data
- `testidx': the set of indices to use as test/validation data, or "nothing". Default is "nothing"
- `dtype::Type`: The data type consistent with the model. Default is Float32.
"""
function make_data(data::DataFrame, input_vars::Vector{Symbol}, target::Symbol, target_scale::Real; testidx=nothing, dtype::Type=Float32)
    vardata = select(data, cat(dims=1, input_vars, target))
    traindata = (isnothing(testidx)) ? vardata : vardata[(!).(testidx), :]
    x_train = Matrix{dtype}(traindata[!, Not(target)])'
    y_train = Vector{dtype}(traindata[!, target])' ./ dtype(target_scale)
    if isnothing(testidx)
        return x_train, y_train
    end
    testdata = vardata[testidx, :]
    x_test = Matrix{dtype}(testdata[!, Not(target)])'
    y_test = Vector{dtype}(testdata[!, target])' ./ dtype(target_scale)
    return x_train, y_train, x_test, y_test
end

"""
    nash_sutcliffe(pred, truth)

Calculate the Nash-Sutcliffe efficiency (NSE) between a predicted and true timeseries.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function nashsutcliffe(pred::Vector{<:Real}, truth::Vector{<:Real})
    return 1.0 - (sum((truth .- pred) .^ 2) / sum((truth .- mean(truth)) .^ 2))
end

"""
    rmse(pred, truth)

Calculate the Root-Mean-Square Error (RMSE) between a predicted and true output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function rmse(pred::Vector{<:Real}, truth::Vector{<:Real})
    return sqrt(Flux.Losses.mse(pred, truth))
end

"""
    bias(pred, truth)

Calculate the bias between a predicted and true output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function bias(pred::Vector{<:Real}, truth::Vector{<:Real})
    return mean(pred .- truth)
end

"""
    std_resid(pred, truth)

Calculate the standard deviation of the residuals between a predicted and true output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function std_resid(pred::Vector{<:Real}, truth::Vector{<:Real})
    return std(pred .- truth)
end

"""
    l1(pred, truth)

Calculate the Mean Absolute Error (L1 Norm) between a predicted and true output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function l1(pred::Vector{<:Real}, truth::Vector{<:Real})
    return Flux.Losses.mae(pred, truth)
end

"""
    trendslope(pred, truth)

Calculate the slope of the trendline between a predicted and true output, with or without an offset.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
-  `offset:Bool`: indicates whether to include an offset instead of a direct proportionality.
"""
function trendslope(pred::Vector{<:Real}, truth::Vector{<:Real}; offset::Bool=false)
    if offset
        return [truth ones(length(truth))] \ pred
    end
    return [truth zeros(length(truth))] \ pred
end

"""
    med_percent_err(pred, truth)

Calculate the median of the percent error between all corresponding elements of a true and predicted output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function med_percent_err(pred::Vector{<:Real}, truth::Vector{<:Real})
    temp = abs.((pred .- truth)) ./ abs.(truth)
    temp = temp[temp.<Inf]
    return median(filter(!isnan, temp))
end

"""
    pack_percent_err(pred, truth)

Returns the L1 error of a predicted and true output, divided by the mean of all true nonzero values, as a measure of relative performance.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function pack_percent_err(pred::Vector{<:Real}, truth::Vector{<:Real})
    return l1(pred, truth) ./ mean(truth[truth.>0.0])
end

"""
    pearson_cor(pred, truth)

Calculate the correlation coefficient between a true and predicted output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
"""
function pearson_cor(pred::Vector{<:Real}, truth::Vector{<:Real})
    return cor(pred, truth)
end

"""
    get_scores(pred, truth; timeseries)

Return a vector of scoring metrics between a true and predicted output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
- `timeseries::Bool`: indicates whether to include metrics for timeseries data
"""
function get_scores(pred::Vector{<:Real}, truth::Vector{<:Real}; timeseries=false)
    scores = zeros(0)
    push!(scores, l1(pred, truth))
    push!(scores, rmse(pred, truth))
    push!(scores, bias(pred, truth))
    push!(scores, std_resid(pred, truth))
    push!(scores, trendslope(pred, truth)[1])
    push!(scores, med_percent_err(pred, truth))
    if timeseries
        push!(scores, nashsutcliffe(pred, truth))
        push!(scores, pack_percent_err(pred, truth))
    end
    return scores
end


"""
    display_scores(pred, truth; timeseries)

Display formatted scoring metrics between a true and predicted output.

# Arguments
- `pred::Vector{<:Real}`: the predicted series.
- `truth::Vector{<:Real}`: the true series.
- `timeseries::Bool`: indicates whether to include metrics for timeseries data
"""
function display_scores(pred::Vector{<:Real}, truth::Vector{<:Real}; timeseries=false)
    print("\n******** SCORES: *******\n")
    print("MAE:          ", l1(pred, truth), "\n")
    print("RMSE:         ", rmse(pred, truth), "\n")
    print("bias:         ", bias(pred, truth), "\n")
    print("std_resid:    ", std_resid(pred, truth), "\n")
    print("trend [m, b]: ", trendslope(pred, truth), "\n")
    print("median % err: ", med_percent_err(pred, truth), "\n")
    if timeseries
        print("NSE:          ", nashsutcliffe(pred, truth), "\n")
        print("pack % err:   ", pack_percent_err(pred, truth), "\n")
    end
end

"""
    feature_importance(data, input_vars, target, loss; dtype, ret, verbose)

Calculate the feature importance metrics via a random shuffling of input features (direct model outuput).

# Arguments
- `data::DataFrame`: the dataframe including input and target data.
- `input_vars::Vector{Symbol}`: the variables used for prediction.
- `target::Symbol`: the variable used as a target.
- `loss`: the loss function used to evaluate model performance.
- `dtype::Type`: the data type consistent with model output. Default is Float32.
- `ret::Bool`: indicates whether to return the scores as a vector. Default is true.
- `verbose::Bool`: indicates whether to print the scores as they are generated. Default is false.
"""
function feature_importance(data::DataFrame, input_vars::Vector{Symbol}, target::Symbol, loss;
    dtype::Type = Float32,
    ret::Bool = true,
    verbose::Bool = false)
    #temp = model[:final_scale].weight[3,3]
    input_base = Matrix{dtype}(select(data, input_vars))'
    output = Vector{dtype}(data[!, target])'
    #setoutscale!(model, 1.0);
    baseline = loss(input_base, output)
    scores = Vector{dtype}()
    for feature in input_vars
        newtrain = deepcopy(data)
        shuffle_feature = Flux.Random.shuffle(newtrain[!, feature])
        newtrain[!, feature] .= shuffle_feature
        newinput = Matrix{dtype}(select(newtrain, input_vars))'
        newloss = loss(newinput, output)
        if verbose print("Importance Score of ", string(feature), ": ", newloss / baseline, "\n") end
        push!(scores, newloss/baseline)
    end
    #setoutscale!(model, temp);
    if ret
        return scores
    end
end

"""
    feature_importance_series(data, input_vars, target, loss; dtype, ret, verbose)

Calculate the feature importance metrics via a random shuffling of input features (model outuput as timeseries).

# Arguments
- `data::DataFrame`: the dataframe including input and target data.
- `input_vars::Vector{Symbol}`: the variables used for prediction.
- `target::Symbol`: the variable used as a target.
- `loss`: the loss function used to evaluate model performance.
- `dtype::Type`: the data type consistent with model output. Default is Float32.
- `ret::Bool`: indicates whether to return the scores as a vector. Default is true.
- `verbose::Bool`: indicates whether to print the scores as they are generated. Default is false.
"""
function feature_importance_series(data::DataFrame, input_vars::Vector{Symbol}, target::Symbol, model, loss; dtype = Float32, ret = true, verbose = false)
    scores = zeros(length(input_vars))
    for site in unique(data[!, :id])
        sitedata = data[data[!, :id] .== site, :]
        true_out = sitedata[!, target]
        pred_zs_base, _, _ = make_timeseries(model, sitedata, target, :date, input_vars, Second(86400))
        baseline = loss(pred_zs_base, true_out)
        for i in 1:length(input_vars)
            feature = input_vars[i]
            newdata = deepcopy(sitedata)
            shuffle_feature = Flux.Random.shuffle(newdata[!, feature])
            newdata[!, feature] .= shuffle_feature
            pred_zs_new, _, _ = make_timeseries(model, newdata, target, :date, input_vars, Second(86400))
            newloss = loss(pred_zs_new, true_out)
            scores[i] += newloss/baseline
        end
    end
    scores /= length(unique(data[!, :id]))
    for i in 1:length(scores)
        if verbose print("Importance Score of ", string(input_vars[i]), ": ", scores[i], "\n") end
    end
    if ret
        return scores
    end
end

"""
    custom_loss(x, y, model, n1, n2)

Creates a loss function to be used during training specified by two hyperparameters n1, n2, as outlined in the paper.

# Arguments
- `x`: the input values.
- `y`: output values to compare to.
- `model`: the model over which to evaluate the loss function.
- `n1::Int`: the hyperparameter dictating the scaling of mismatch error.
- `n2::Int`: the hyperparameter dictating the weighting of a given mismatch error by the target magnitude.
"""
function custom_loss(x, y, model, n1, n2)
    sum(abs.(model(x) - y) .^ n1 .* (1.0 .+ abs.(y) .^ n2)) ./ length(y)
end

"""
    trainmodel!(model, ps, x_train, y_train, n1, n2; nepochs, opt, verbose, cb)

A training function for a neural model, permitting usage of a callback function.

# Arguments
- `model`: the model used for training.
- `ps`: the model parameters that will be trained.
- `x_train`: the input training data to be used in training.
- `y_train`: the target data to be used in training.
- `n1`: the scaling hypermarameter used to generate custom loss functions.
- `n2`: the weighting hyperparameter used to generate custom loss functions.
- `nepochs::Int`: the number of epochs. Default is 100.
- `opt`: the Flux optimizer to be used. Default is RMSProp()
- `verbose::Bool`: indicates whether to print the training loss every 10 epochs
- `cb`: Allows utlization of a callback function (must take no input arguments)
"""
function trainmodel!(model, ps, x_train, y_train, n1, n2;
    nepochs::Int = 100,
    opt = Flux.Optimise.RMSProp(),
    verbose = false,
    cb = Nothing)
    train_loader = DataLoader((x_train, y_train), batchsize=nbatch, partial=false, shuffle=true)
    loss(x, y) = custom_loss(x, y, model, n1, n2)
    for epoch in 1:nepochs
        for (x, y) in train_loader
            Flux.train!(loss, ps, [(x, y)], opt)
        end
        if verbose & (epoch%10 == 0)
            print("Epoch: ", epoch, " | training loss: ", loss(x_train, y_train), "\n")
        end
        cb()
    end    
end