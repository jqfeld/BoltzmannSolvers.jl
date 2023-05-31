struct LoKI <: Solver end



function default_names(::LoKI, names=Dict{String,Symbol}();
    replace=false)
    default_names = Dict{String,Symbol}(
        [
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
    ])
    if replace == true
        return names
    else
        return merge(default_names, names)

    end
end

function parse_reaction_names(::LoKI, source, replace_strings=[])
    reaction_dict = Dict{String,Symbol}()
    rate_table_file = joinpath(source, "lookUpTableRateCoeff.txt")
    lines = filter(startswith('#'), readlines(rate_table_file))
    for l in lines
        m = match(r"(\d+)\s+([^\-\<]+)(<->|->|<-)(.*),(.+?)\s", l)
        if isnothing(m)
            continue
        end
        id, lhs, dir, rhs, type = m
        if type == "Effective"
            key = "R$(id)_ine(m^3s^-1)"
            value = "Effective($lhs)"
            reaction_dict[key] = Symbol(replace(value, replace_strings...)) |> normalize_reaction_name!
            continue
        end
        if dir == "->"
            key = "R$(id)_ine(m^3s^-1)"
            value = "$(lhs)-->$(rhs)"
            reaction_dict[key] = Symbol(replace(value, replace_strings...)) |> normalize_reaction_name!
        elseif dir == "<->"
            key = "R$(id)_ine(m^3s^-1)"
            value = "$(lhs)-->$(rhs)"
            reaction_dict[key] = Symbol(replace(value, replace_strings...)) |> normalize_reaction_name!
            key = "R$(id)_sup(m^3s^-1)"
            value = "$(lhs)<--$(rhs)"
            reaction_dict[key] = Symbol(replace(value, replace_strings...)) |> normalize_reaction_name!
        end
    end
    return reaction_dict
end


function load_dataframe(::LoKI, source;
    replace_strings=[],
    names=default_names(LoKI())
)
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
    rename!(df_rate_coef, parse_reaction_names(LoKI(), source, replace_strings))
    return rename!(
        innerjoin(df_swarm, df_rate_coef, on="RedField(Td)"),
        names
    )
end

