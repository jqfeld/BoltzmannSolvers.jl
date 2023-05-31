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

function parse_reaction_names(l::LoKI, source, kwargs...)

    reaction_names = Pair[]
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
            push!(reaction_names, key => value)
            continue
        end
        if dir == "->"
            key = "R$(id)_ine(m^3s^-1)"
            value = "$(lhs)-->$(rhs)"
            push!(reaction_names, key => value)
        elseif dir == "<->"
            key = "R$(id)_ine(m^3s^-1)"
            value = "$(lhs)-->$(rhs)"
            push!(reaction_names, key => value)
            key = "R$(id)_sup(m^3s^-1)"
            value = "$(lhs)<--$(rhs)"
            push!(reaction_names, key => value)
        end
    end
    return reaction_names

end


function load_raw_dataframe(l::LoKI, source;
    kwargs...
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

    return innerjoin(df_swarm, df_rate_coef, on="RedField(Td)")
end

