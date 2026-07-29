abstract type Solver end

include("solvers/loki.jl")
export LoKI

include("solvers/multibolt.jl")
export MultiBolt

include("solvers/bolsig.jl")
export BOLSIG

include("solvers/bolsig_input.jl")
export BOLSIGInput, BOLSIGReadCollisions, BOLSIGConditions, BOLSIGSaveResults
export BOLSIGRunSpec, BOLSIGFixedRun, BOLSIGExplicitRun, BOLSIGSeriesSegment, BOLSIGSeriesRun, BOLSIGRun2D
export BOLSIGVariable, ReducedFieldVar, MeanEnergyVar, MaxwellianEnergyVar
export BOLSIGSeriesType, LinearSeries, QuadraticSeries, ExponentialSeries
export read_bolsig_input, write_bolsig_input

include("solvers/bolsig_run.jl")
export BOLSIGRunResult, run_bolsig


