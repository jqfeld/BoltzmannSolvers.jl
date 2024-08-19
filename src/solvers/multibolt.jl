struct MultiBolt <: Solver end


function default_swarm_names(::MultiBolt)
    default_names =
        [
        "E_N" => :reduced_field,
        "muN_FLUX" => :reduced_mobility,
        "DLN_FLUX" => :longitudinal_diff_coef,
        "DTN_FLUX" => :transversal_diff_coef,
        "alpha_eff_N" => :reduced_townsend_alpha_coef,
        "avg_en" => :mean_energy,
    ]
end


function load_raw_dataframe(::MultiBolt, source; kwargs...)
    swarm_param_files = [
        "muN_FLUX.txt",
        "DTN_FLUX.txt",
        "DLN_FLUX.txt",
        "avg_en.txt",
        "alpha_eff_N.txt",
    ] .|> x -> joinpath(source, x)

    df = innerjoin(
        CSV.read.(swarm_param_files, DataFrame,
            comment="#",
            delim="\t",
        )...,
        on="E_N"
    )

    rates_directory = joinpath(source, "PerGas")
    gas_list = readdir(rates_directory)
    for gas in gas_list
        rate_filenames = readdir(joinpath(rates_directory, gas))
        df = map(rate_filenames) do x
            file = joinpath(rates_directory, gas, x)
            line = filter(startswith("# Process:"), readlines(file))[1]
            lhs, dir, rhs, type = match(r"(?:Process:\s*)(.*?)(->|<->)(.*?)(?:,\s*(Attachment|Excitation|Ionization|Elastic|Effective))$", line)
            if dir != "->"
                error("Only '->' implemented.")
            end
            reaction_name = replace("$(lhs)-->$(rhs)", "E" => "e", " " => "")
            if startswith(x, "alpha") 
                reaction_name = "alpha($reaction_name)"
            end
            if startswith(x, "eta") 
                reaction_name = "eta($reaction_name)"
            end
            CSV.read(file, DataFrame,
                comment="#",
                delim="\t",
                skipto=7,
                header=["E_N", reaction_name]
            )
        end |> x -> innerjoin(df, x..., on="E_N")
    end

    return df
end

"""
    parse_reaction_names(::Multibolt, source)

Do nothing. For MultiBolt it is more convenient to parse the reaction process already in `load_raw_dataframe`.
"""
parse_reaction_names(::MultiBolt, _) = Pair[]
