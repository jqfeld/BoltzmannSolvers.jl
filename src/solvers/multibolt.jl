struct MultiBolt <: Solver end

MULTIBOLT_NAMES = Dict{String,Symbol}(
    [
    "E_N" => :reduced_field,
    # "RedDiff((ms)^-1)" => :reduced_diffusion_coef,
    "muN_FLUX" => :reduced_mobility,
    # "DriftVelocity(ms^-1)" => :drift_velocity,
    # "RedTow(m^2)" => :reduced_townsend_ion_coef,
    # "RedAtt(m^2)" => :reduced_attachment_coef,
    # "RedDiffE(eV(ms)^-1)" => :reduced_energy_diffusion_coef,
    # "RedMobE(eV(msV)^-1)" => :reduced_energy_mobility,
    "avg_en" => :mean_energy,
    # "CharE(eV)" => :characteristic_energy,
    # "EleTemp(eV)" => :electron_temperature,
])

function load_dataframe(::MultiBolt, source, names=MULTIBOLT_NAMES)
    swarm_param_files = [
        "muN_FLUX.txt",
        "avg_en.txt"] .|> x -> joinpath(source, x)

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
            CSV.read(joinpath(rates_directory, gas, x), DataFrame,
                comment="#",
                delim="\t",
                skipto=7,
                header=[
                    "E_N",
                    (
                        length(gas_list) > 1 ?
                        replace(gas, "_" => "") * "_" :
                        ""
                    ) *
                    replace(x,
                        ".txt" => "",
                        "alpha_N_0" => "alpha_N",
                        "_" => "",
                        # count=2
                    )]
            )
        end |> x -> innerjoin(df, x..., on="E_N")
    end

    return rename!(df, names)
end

