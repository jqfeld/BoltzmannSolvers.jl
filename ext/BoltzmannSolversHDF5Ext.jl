module BoltzmannSolversHDF5Ext
using BoltzmannSolvers
using HDF5
using DataFrames

# Extracts (lhs, dir, rhs, type) from a LoKI-B collision `description` string
# read directly from the HDF5 `rateCoefficients` dataset, e.g.
# "e+N2(X,v=0)<->e+N2(X,v=1),Vibrational". Same shape as the regex the `.txt`
# reading path matches against a `#`-comment line (src/solvers/loki.jl), just
# anchored at both ends since `description` is already a clean, isolated
# string (no leading id, no trailing padding/`#`).
function _parse_loki_hdf5_description(desc::AbstractString)
    m = match(r"^([^\-\<]+)(<->|->|<-)(.*),(.+?)\s*$", desc)
    isnothing(m) && return nothing
    return m[1], m[2], m[3], m[4]
end

function BoltzmannSolvers._load_loki_hdf5(source::AbstractString)
    h5open(source, "r") do f
        g = f["electronKinetics"]

        reduced_field = vec(read(g["reducedField"]))
        swarm = vec(read(g["swarmParameters"]))

        df = DataFrame(
            "RedField(Td)" => reduced_field,
            "MeanE(eV)" => getproperty.(swarm, :meanEnergy),
            "CharE(eV)" => getproperty.(swarm, :characEnergy),
            "EleTemp(eV)" => getproperty.(swarm, :Te),
            "DriftVelocity(ms^-1)" => getproperty.(swarm, :driftVelocity),
            "RedMob((msV)^-1)" => getproperty.(swarm, :redMobility),
            "RedDiff((ms)^-1)" => getproperty.(swarm, :redDiffCoeff),
            "RedMobE(eV(msV)^-1)" => getproperty.(swarm, :redMobilityEnergy),
            "RedDiffE(eV(ms)^-1)" => getproperty.(swarm, :redDiffCoeffEnergy),
            "RedTow(m^2)" => getproperty.(swarm, :redTownsendCoeff),
            "RedAtt(m^2)" => getproperty.(swarm, :redAttCoeff),
        )

        # rateCoefficients: (n_reactions, 1, n_conditions), one compound
        # record per (reaction, condition) with the coefficient(s) and a
        # self-describing `description` — decide ine-only vs ine+sup from the
        # parsed `dir`, not from whether `sup_coeff` happens to be zero (a
        # genuinely reversible reaction can have a momentarily negligible
        # backward rate at some condition).
        rate = read(g["rateCoefficients"])
        n_reactions, _, n_conditions = size(rate)
        for r in 1:n_reactions
            id = rate[r, 1, 1].rate_id
            parsed = _parse_loki_hdf5_description(rate[r, 1, 1].description)
            isnothing(parsed) && continue
            lhs, dir, rhs, type = parsed
            df[!, "R$(id)_ine(m^3s^-1)"] = [rate[r, 1, c].ine_coeff for c in 1:n_conditions]
            if dir == "<->"
                df[!, "R$(id)_sup(m^3s^-1)"] = [rate[r, 1, c].sup_coeff for c in 1:n_conditions]
            end
        end

        # powerBalanceSummary: (n_conditions, 1, 3), 3rd axis = (Net, Gain,
        # Loss) — confirmed by cross-matching against the already-tested
        # lookUpTablePower.txt values (see plan/PR notes). Optional, same as
        # lookUpTablePower.txt being an optional join on the `.txt` path.
        if haskey(g, "powerBalanceSummary")
            pbs = read(g["powerBalanceSummary"])
            net, gain, loss = pbs[:, 1, 1], pbs[:, 1, 2], pbs[:, 1, 3]

            df[!, "PowerField(eVm^3s^-1)"] = getproperty.(net, :Field)
            for (field, prefix) in (
                (:Elastic, "Ela"), (:CAR, "CAR"), (:Electronic, "Ele"),
                (:Vibrational, "Vib"), (:Rotational, "Rot"),
            )
                df[!, "Pwr$(prefix)Gain(eVm^3s^-1)"] = getproperty.(gain, field)
                df[!, "Pwr$(prefix)Loss(eVm^3s^-1)"] = getproperty.(loss, field)
                df[!, "Pwr$(prefix)Net(eVm^3s^-1)"] = getproperty.(net, field)
            end
            # Not gain/loss-split in the `.txt` convention either — single
            # (Net-slot) columns only.
            df[!, "PwrIon(eVm^3s^-1)"] = getproperty.(net, :Ionization)
            df[!, "PwrAtt(eVm^3s^-1)"] = getproperty.(net, :Attachment)
            df[!, "PwrGroth(eVm^3s^-1)"] = getproperty.(net, :eDensGrowth)
            df[!, "PwrBalance(eVm^3s^-1)"] = getproperty.(net, :Balance)
            # The Gain slot's "Balance" field holds the *relative* balance as
            # a fraction; the `.txt` file's "RelPwrBalance" is the same value
            # expressed as a percentage.
            df[!, "RelPwrBalance"] = getproperty.(gain, :Balance) .* 100
        end

        return df
    end
end

function BoltzmannSolvers._parse_loki_hdf5_reactions(source::AbstractString)
    reaction_names = Pair{String,String}[]
    h5open(source, "r") do f
        rate = read(f["electronKinetics"]["rateCoefficients"])
        n_reactions = size(rate, 1)
        for r in 1:n_reactions
            id = rate[r, 1, 1].rate_id
            parsed = _parse_loki_hdf5_description(rate[r, 1, 1].description)
            isnothing(parsed) && continue
            lhs, dir, rhs, type = parsed
            append!(reaction_names, BoltzmannSolvers._loki_reaction_columns(id, lhs, dir, rhs, type))
        end
    end
    return reaction_names
end

end
