module BoltzmannSolvers

using CSV, DataFrames

include("solvers.jl")
include("interpolation.jl")

load_raw_dataframe(::Solver, source) = error("Not implemented for this solver.")
default_swarm_names(::Solver) = error("Not implemented for this solver.")


function load_dataframe(s::S, source; kwargs...) where S <: Solver
    df = load_raw_dataframe(s, source; kwargs...)
    rename!(df, default_swarm_names(s)...)
    reaction_names = parse_reaction_names(s, source)
    rename!(df, reaction_names...)
    return df
end 

export load_dataframe, default_names
export create_interpolation

end
