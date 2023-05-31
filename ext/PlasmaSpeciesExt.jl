module PlasmaSpeciesExt
using BoltzmannSolvers
using PlasmaSpecies


"""
    normalize_reaction_name!(x)

Normalize the reaction name using `PlasmaSpecies.jl`. If `PlasmaSpecies.jl` is not loaded, this is a NOP.
PlasmaSpecies is loaded!
"""
function BoltzmannSolvers.normalize_reaction_name!(str::String)
    string(parse_reaction(str)) 
end

end
