module BoltzmannSolvers

using CSV, DataFrames

include("solvers.jl")
include("interpolation.jl")

load_dataframe(::Solver, source) = error("Not implemented for this solver.")
export load_dataframe, default_names
export create_interpolation

"""
    normalize_reaction_name!(x)

Normalize the reaction name using `PlasmaSpecies.jl`. 
If `PlasmaSpecies.jl` is not loaded, this is a NOP.
PlasmaSpecies is NOT loaded!
"""
normalize_reaction_name!(x) = x

end
