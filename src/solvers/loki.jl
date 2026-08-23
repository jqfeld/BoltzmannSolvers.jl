struct LoKI <: Solver end



function default_swarm_names(l::LoKI)
    default_names = [
        "RedField(Td)" => :reduced_field,
        "RedDiff((ms)^-1)" => :reduced_diffusion_coef,
        "RedMob((msV)^-1)" => :reduced_mobility,
        "DriftVelocity(ms^-1)" => :drift_velocity,
        "RedTow(m^2)" => :reduced_townsend_alpha_coef,
        "RedAtt(m^2)" => :reduced_attachment_coef,
        "RedDiffE(eV(ms)^-1)" => :reduced_energy_diffusion_coef,
        "RedMobE(eV(msV)^-1)" => :reduced_energy_mobility,
        "MeanE(eV)" => :mean_energy,
        "CharE(eV)" => :characteristic_energy,
        "EleTemp(eV)" => :electron_temperature,
    ]
end

# Given a parsed collision (lhs, dir, rhs, type), build the rate-coefficient
# column-name pairs ("R<id>_ine(m^3s^-1)" [, "R<id>_sup(m^3s^-1)"] => reaction
# name). `type` only matters for the "Effective" special case — every other
# type string (Vibrational, Rotational, Excitation, Ionization, Elastic, ...)
# falls through the same `dir`-based branch. Shared between the `.txt`
# reading path (parsed from a `#`-comment line) and the `.h5` reading path
# (parsed from a per-reaction `description` field, see
# `ext/BoltzmannSolversHDF5Ext.jl`) so both produce identical reaction-name
# conventions regardless of source format.
function _loki_reaction_columns(id, lhs, dir, rhs, type)
    type == "Effective" && return ["R$(id)_ine(m^3s^-1)" => "Effective($lhs)"]
    dir == "->" && return ["R$(id)_ine(m^3s^-1)" => "$(lhs)-->$(rhs)"]
    if dir == "<->"
        return [
            "R$(id)_ine(m^3s^-1)" => "$(lhs)-->$(rhs)",
            "R$(id)_sup(m^3s^-1)" => "$(lhs)<--$(rhs)",
        ]
    end
    return Pair{String,String}[]
end

function _parse_loki_txt_reactions(source)
    reaction_names = Pair{String,String}[]
    rate_table_file = joinpath(source, "lookUpTableRateCoeff.txt")
    lines = filter(startswith('#'), readlines(rate_table_file))

    for l in lines
        # non-greedy lhs up to the first arrow: species labels may contain
        # '-' themselves (e.g. "O(-,gnd)", "O2(A3Su+_C3Du_c1Su-)")
        m = match(r"(\d+)\s+(.+?)(<->|->|<-)(.*),(.+?)\s", l)
        isnothing(m) && continue
        id, lhs, dir, rhs, type = m
        append!(reaction_names, _loki_reaction_columns(id, lhs, dir, rhs, type))
    end
    return reaction_names
end

"""
    _load_loki_hdf5(source)

Overridden by `ext/BoltzmannSolversHDF5Ext.jl` once `HDF5.jl` is loaded.
"""
_load_loki_hdf5(source) = error(
    "Reading LoKI-B HDF5 output requires HDF5.jl — run `using HDF5` " *
    "before calling load_dataframe(LoKI(), \"...h5\")."
)

"""
    _parse_loki_hdf5_reactions(source)

Overridden by `ext/BoltzmannSolversHDF5Ext.jl` once `HDF5.jl` is loaded.
"""
_parse_loki_hdf5_reactions(source) = error(
    "Reading LoKI-B HDF5 output requires HDF5.jl — run `using HDF5` " *
    "before calling load_dataframe(LoKI(), \"...h5\")."
)

"""
    parse_reaction_names(::LoKI, source)

`source` is either a directory containing `lookUpTableRateCoeff.txt` (LoKI-B's
plain-text lookup-table output) or a single `.h5` file (LoKI-B's HDF5 export —
requires `using HDF5`, see [`_load_loki_hdf5`](@ref)).
"""
function parse_reaction_names(l::LoKI, source, kwargs...)
    isdir(source) && return _parse_loki_txt_reactions(source)
    isfile(source) && return _parse_loki_hdf5_reactions(source)
    error("LoKI: `source` must be a directory (lookUpTable*.txt) or a single HDF5 file.")
end


function _load_loki_txt(source)
    swarm_table_file = joinpath(source, "lookUpTableSwarm.txt")
    df_swarm = CSV.read(swarm_table_file, DataFrame,
        comment="#",
        delim=" ",
        ignorerepeated=true
    )

    rate_table_file = joinpath(source, "lookUpTableRateCoeff.txt")
    df_rate_coef = CSV.read(rate_table_file, DataFrame,
        comment="#",
        delim=" ",
        ignorerepeated=true
    )

    df = innerjoin(df_swarm, df_rate_coef, on="RedField(Td)")

    # lookUpTablePower.txt (per-channel electron energy gain/loss balance) is
    # a newer addition to LoKI-B's output and isn't guaranteed to exist in
    # every simulation directory (e.g. older runs, or power-balance output
    # disabled), so it's joined in only when present rather than required.
    power_table_file = joinpath(source, "lookUpTablePower.txt")
    if isfile(power_table_file)
        df_power = CSV.read(power_table_file, DataFrame,
            comment="#",
            delim=" ",
            ignorerepeated=true
        )
        # "RelPwrBalance" is written with a literal trailing "%" glued onto
        # the number (e.g. "6.95053583371272e-13%"), which defeats CSV.jl's
        # own Float64 auto-detection and leaves the column as a String.
        if eltype(df_power[!, "RelPwrBalance"]) <: AbstractString
            df_power[!, "RelPwrBalance"] = parse.(Float64, rstrip.(df_power[!, "RelPwrBalance"], '%'))
        end
        df = innerjoin(df, df_power, on="RedField(Td)")
    end

    return df
end

"""
    load_raw_dataframe(::LoKI, source; kwargs...)

`source` is either a directory containing LoKI-B's plain-text lookup tables
(`lookUpTableSwarm.txt`, `lookUpTableRateCoeff.txt`, and optionally
`lookUpTablePower.txt`) or a single `.h5` file — LoKI-B's HDF5 export of the
same simulation (requires `using HDF5`, see [`_load_loki_hdf5`](@ref)).
"""
function load_raw_dataframe(l::LoKI, source;
    kwargs...
)
    isdir(source) && return _load_loki_txt(source)
    isfile(source) && return _load_loki_hdf5(source)
    error("LoKI: `source` must be a directory (lookUpTable*.txt) or a single HDF5 file.")
end
