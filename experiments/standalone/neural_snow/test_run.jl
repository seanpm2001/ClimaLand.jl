print("\nCode Starting, importing libraries!\n")
include("data_tools.jl")
include("model_tools.jl")
include("train_test_tools.jl")
include("display_tools.jl")
using DataFrames, CSV, Dates

print("\nRunning Evaluation!\n")
const nepochs = 100
n = 5
n1 = 1
n2 = 1

pred_vars = [
    :z,
    :SWE,
    :rel_hum_avg,
    :sol_rad_avg,
    :wind_speed_avg,
    :air_temp_avg,
    :dprecipdt_snow,
]
target = :dzdt
z_idx = 1
p_idx = 7
Δt = Second(86400)
const nfeatures = length(pred_vars)

data_filename = "cleandata.csv"
data = CSV.read(data_filename, DataFrame)
usedata = prep_data(data)
out_scale = maximum(abs.(usedata[!, target]))
in_scales = std.(eachcol(select(usedata, pred_vars)))
x_train, y_train = make_data(usedata, pred_vars, target, out_scale)

model = make_model(nfeatures, n, z_idx, p_idx, in_scale = in_scales)
ps = get_model_ps(model)

#training loop:
settimescale!(model, Dates.value(Δt) * out_scale)
setoutscale!(model, 1.0)
print("\nTraining model!\n")
trainmodel!(model, ps, x_train, y_train, n1, n2, verbose = true)

#validation/usage
print("\nGenerating Timeseries!\n")
setoutscale!(model, out_scale)
settimescale!(model, Dates.value(Δt))
for site in unique(data[!, :id])
    sitedata = usedata[usedata[!, :id] .== site, :]
    true_series = sitedata[!, :z]
    pred_series, _, _ = make_timeseries(model, sitedata, Δt)
    print("\nSITE: ", site)
    display_scores(pred_series, true_series, timeseries = true)
    siteplot(
        sitedata[!, :date],
        pred_series,
        true_series,
        "SITE: " * string(site),
    )
end
