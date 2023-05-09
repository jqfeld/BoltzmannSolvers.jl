struct LoKI <: Solver end


LOKI_NAMES = Dict{String,Symbol}(
    [
    "RedField(Td)" => :reduced_field,
    "RedDiff((ms)^-1)" => :reduced_diffusion_coef,
    "RedMob((msV)^-1)" => :reduced_mobility,
    "DriftVelocity(ms^-1)" => :drift_velocity,
    "RedTow(m^2)" => :reduced_townsend_ion_coef,
    "RedAtt(m^2)" => :reduced_attachment_coef,
    "RedDiffE(eV(ms)^-1)" => :reduced_energy_diffusion_coef,
    "RedMobE(eV(msV)^-1)" => :reduced_energy_mobility,
    "MeanE(eV)" => :mean_energy,
    "CharE(eV)" => :characteristic_energy,
    "EleTemp(eV)" => :electron_temperature,
])

function load_dataframe(::LoKI, source, names=LOKI_NAMES)
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
    return rename!(
        innerjoin(df_swarm, df_rate_coef, on="RedField(Td)"),
        names
    )
end

