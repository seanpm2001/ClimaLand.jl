import ClimaAnalysis

"""
    get_data_sources(diagnostics_folder_path)

Get the necessary data and do the preprocessing for computing the leaderboard. Return
`sim_var_dict` which contains simulation data, `obs_var_dict` which contains observational
data, `mask_dict` which contains masks, `compare_vars_biases_plot_extrema` which contains
the ranges for the colormaps.

To add a variable for the leaderboard, add a key-value pair to the dictionary
`sim_var_dict` whose key is the short name of the variable and the value is a function
that returns a `OutputVar` Any preprocessing is done in the function which includes unit
conversion and shifting the dates.

Then, add a key-value pair to the dictionary `obs_var_dict` whose key is the same short
name as before and the value is a function that takes in a start date and returns a
`OutputVar`. Any preprocessing is done in the function.

Next, add a key-value pair to the dictionary `mask_dict` whose key is the same short name
as before and the value is a function that takes in a `OutputVar` representing simulation
data and a `OutputVar` representing observational data and returns a masking function or
`nothing` if no masking function is needed. The masking function is used to correctly
normalize the global bias and global RMSE.

Finally, add a key-value pair to the dictionary `compare_vars_biases_plot_extrema` whose
key is the same short name as before and the value is a tuple of floats which determine
the range of the bias plots.
"""
function get_data_sources(diagnostics_folder_path)
    # For plotting bias plots
    compare_vars_biases_plot_extrema = Dict(
        "et" => (-0.00001, 0.00001),
        "gpp" => (-0.000005, 0.000005),
        "lwu" => (-40.0, 40.0),
    )

    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()

    for short_name in ["et", "lwu"]
        sim_var_dict[short_name] =
            () -> begin
                sim_var = get(
                    ClimaAnalysis.SimDir(diagnostics_folder_path),
                    short_name = short_name,
                )
                # Remove the line below later when start_date is added to the diagnostics
                sim_var =
                    ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
                return sim_var
            end
    end

    sim_var_dict["gpp"] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "gpp",
            )
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end

    # Dict for loading in observational data
    obs_var_dict = Dict{String, Any}()
    obs_var_dict["et"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                ClimaLand.Artifacts.ilamb_dataset_path(;
                    context = "evspsbl_MODIS_et_0.5x0.5.nc",
                ),
                "et",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            ClimaAnalysis.units(obs_var) == "kg/m2/s" &&
                (obs_var = ClimaAnalysis.set_units(obs_var, "kg m^-2 s^-1"))
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)
            return obs_var
        end

    obs_var_dict["gpp"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                ClimaLand.Artifacts.ilamb_dataset_path(;
                    context = "gpp_FLUXCOM_gpp.nc",
                ),
                "gpp",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            ClimaAnalysis.dim_units(obs_var, "lon") == "degree" &&
                (obs_var.dim_attributes["lon"]["units"] = "degrees_east")
            ClimaAnalysis.dim_units(obs_var, "lat") == "degree" &&
                (obs_var.dim_attributes["lat"]["units"] = "degrees_north")
            # converting from `g C m-2 day-1` in obs to `mol CO2 m^-2 s^-1` in sim
            obs_var = ClimaAnalysis.convert_units(
                obs_var,
                "mol CO2 m^-2 s^-1",
                conversion_function = units -> (units / 86400.0) / 12.011,
            )
            obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)
            return obs_var
        end

    obs_var_dict["lwu"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                ClimaLand.Artifacts.ilamb_dataset_path(;
                    context = "rlus_CERESed4.2_rlus.nc",
                ),
                "rlus",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            ClimaAnalysis.units(obs_var) == "W m-2" &&
                (obs_var = ClimaAnalysis.set_units(obs_var, "W m^-2"))
            return obs_var
        end

    # Dict for loading in masks
    mask_dict = Dict{String, Any}()

    mask_dict["et"] =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.make_lonlat_mask(
                ClimaAnalysis.slice(
                    obs_var,
                    time = ClimaAnalysis.times(obs_var) |> first,
                );
                set_to_val = isnan,
            )
        end

    mask_dict["gpp"] =
        (sim_var, obs_var) -> begin
            return ClimaAnalysis.make_lonlat_mask(
                ClimaAnalysis.slice(
                    obs_var,
                    time = ClimaAnalysis.times(obs_var) |> first,
                );
                set_to_val = isnan,
            )
        end

    mask_dict["lwu"] = (sim_var, obs_var) -> begin
        return nothing
    end

    return sim_var_dict,
    obs_var_dict,
    mask_dict,
    compare_vars_biases_plot_extrema
end
