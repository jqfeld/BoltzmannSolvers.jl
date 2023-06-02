module BoltzmannSolversPlasmaSpeciesExt
using BoltzmannSolvers
using PlasmaSpecies


"""
    normalize_reaction_names(x)

Normalize the reaction name using `PlasmaSpecies.jl`. If `PlasmaSpecies.jl` is not loaded, this is a NOP.
PlasmaSpecies is loaded!
"""
function BoltzmannSolvers.normalize_reaction_names(str::String)
    out = str
    try
        out = string(parse_reaction(str)) 
    catch e
        @warn "Error while parsing" str out
    end
    return out
end

end
