print("Code Starting!\n")
print("Importing Libraries!\n")
include("data_tools.jl")
using Plots
############## Helper functions for code ##################

function desc(input)
    return show(describe(input), allrows = true, allcols = true)
end

function siteplot(data, id)
    vars = ["SWE", "z", "rel_hum_avg", "sol_rad_avg", "wind_speed_avg", "air_temp_avg", "dprecipdt"]
    for var in vars
        t = string(id)*": "*var
        display(plot(data[!, :date], data[!, Symbol(var)], title = t))
    end
end

const inch2meter = 0.0254
const kmphr2mps = 5.0/18.0

filter_val = Dict{Symbol,Tuple{Real,Real}}(
    :SWE => (0.0, 250.0),
    :z => (0.0, 420.0),
    :precip => (0.0, 250.0),
    :rel_hum_avg => (10.0, 100.0),
    :sol_rad_avg => (0.0, 1500.0),
    :wind_speed_avg => (0.0, 216.0),
    :air_temp_avg => (-55.0, 60.0)
)

scales = Dict{Symbol, Real}(
    :SWE => inch2meter,
    :z => inch2meter,
    :precip => inch2meter,
    :rel_hum_avg => 0.01,
    :wind_speed_avg => kmphr2mps,
)

good_stations = Dict{Int, Tuple{String, String}}(
    1030 => ("start", "2016-01-01"),
    1053 => ("2010-01-01", "end"),
    1083 => ("2013-01-01", "end"),
    1105 => ("start", "2012-06-01"),
    1122 => ("start", "end"),
    1123 => ("start", "end"),
    1159 => ("start", "end"),
    1168 => ("start", "end"),
    1170 => ("start", "end"),
    1254 => ("2018-01-01", "end"),
    1286 => ("start", "end"),
    306 => ("2013-01-01", "2020-01-01"),
    316 => ("start", "end"),
    344 => ("start", "end"),
    367 => ("2022-09-01", "2023-04-01"),
    395 => ("start", "end"),
    457 => ("2008-01-01", "end"),
    482 => ("start", "2010-01-01"),
    491 => ("2013-01-01", "end"),
    532 => ("2013-01-01", "2022-01-01"),
    551 => ("start", "end"),
    571 => ("start", "end"),
    599 => ("start", "end"),
    608 => ("start", "2018-01-01"),
    613 => ("2007-01-01", "2015-01-01"),
    665 => ("start", "2022-01-01"),
    734 => ("start", "end"),
    737 => ("start", "end"),
    744 => ("2014-01-01", "2016-01-01"),
    832 => ("start", "end"),
    845 => ("2014-01-01", "end"),
    854 => ("2019-01-01", "end"),
    857 => ("start", "end"),
    921 => ("2019-01-01", "end"),
    922 => ("2015-01-01", "2018-01-01"),
    942 => ("2009-01-01", "2018-01-01"),
    969 => ("2011-01-01", "end"),
    974 => ("start", "end"),
    978 => ("2005-01-01", "2018-01-01"),
)

station_metadata_site = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customMultipleStationReport,metric/daily/start_of_period/network=%22SNTL%22%7Cname/0,0/stationId,state.code,elevation,latitude,longitude"
metadata = CSV.read(HTTP.get(station_metadata_site).body, DataFrame, comment="#", delim=",")
metacols = ["id", "state", "elev", "lat", "lon"]
DataFrames.rename!(metadata, Symbol.(metacols))

allsites = Any[]
for site in sort(collect(keys(good_stations)))
    state = metadata[metadata[!, :id].==site, :state][1]
    start_date = good_stations[site][1]
    end_date = good_stations[site][2]
    print("\nSITE: ", site, " (", state, ") | FROM ", start_date, " TO ", end_date, "\n")
    daily = apply_bounds(sitedata_daily(site, state, start = start_date, finish = end_date), filter_val)
    hourly = apply_bounds(sitedata_hourly(site, state, start = start_date, finish = end_date), filter_val)
    hourly_d = hourly2daily(hourly)
    daily = scale_cols(rectify_daily_hourly(daily, hourly_d), scales)
    #daily = scale_cols(hourly, scales)
    daily_clean = daily[completecases(daily), :]
    daily_clean = makediffs(daily_clean, Day(1))
    #daily_clean = makediffs(daily_clean, Hour(1))
    good_vals = daily_clean[!, :dprecipdt] .>= 0.0
    daily_clean = daily_clean[good_vals, Not(:precip)]
    #daily_r = rolldata(daily_clean, Day(1), 7)
    siteplot(daily_clean, site)
    desc(daily_clean)
    print("\nSIZE: ", nrow(daily_clean), "\n")

    daily_clean[!, :id] .= site
    daily_clean[!, :elev] .= metadata[metadata[!, :id].==site, :elev][1]
    daily_clean[!, :lat] .= metadata[metadata[!, :id].==site, :lat][1]
    daily_clean[!, :lon] .= metadata[metadata[!, :id].==site, :lon][1]

    push!(allsites, daily_clean)
end

totaldata = deepcopy(allsites[1])
for site in allsites[2:end]
    append!(totaldata, site)
end
#CSV.write("newcleanfiletestdata.csv", totaldata)