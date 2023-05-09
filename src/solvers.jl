abstract type Solver end

include("solvers/loki.jl")
export LoKI

include("solvers/multibolt.jl")
export MultiBolt

#TODO: implement bolsig and multibolt
include("solvers/bolsig.jl")


