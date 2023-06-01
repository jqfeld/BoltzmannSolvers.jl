module BoltzmannSolvers

using CSV, DataFrames

include("solvers.jl")
include("interpolation.jl")

load_raw_dataframe(::Solver, source) = error("Not implemented for this solver.")
default_swarm_names(::Solver) = error("Not implemented for this solver.")

normalize(x) = x

function load_dataframe(s::S, source; 
    replacements=nothing, 
    normalize=false,
    kwargs...) where S <: Solver
    df = load_raw_dataframe(s, source; kwargs...)
    rename!(df, default_swarm_names(s)...)
    reaction_names = parse_reaction_names(s, source)
    rename!(df, reaction_names...)
    if !isnothing(replacements)
        rename!(df, names(df) .=> replace.(names(df), replacements))
    end 

    if normalize
        rename!(df, names(df) .=> normalize.(names(df)))
    end

    return df
end 

export load_dataframe
export create_interpolation

end
