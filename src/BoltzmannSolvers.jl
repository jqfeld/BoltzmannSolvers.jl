module BoltzmannSolvers

using CSV, DataFrames

include("solvers.jl")
include("interpolation.jl")

load_dataframe(::Solver, source) = error("Not implemented for this solver.")
export load_dataframe, default_names
export create_interpolation

end
