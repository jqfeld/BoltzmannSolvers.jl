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
export BOLSIGRunResult, run_solver

include("solvers/multibolt_input.jl")
export MultiBoltInput, MultiBoltSpecies, MultiBoltXsecScale
export MultiBoltSweep, MultiBoltSweepStyle, MultiBoltLinearSweep, MultiBoltLogSweep, MultiBoltRegularSweep, MultiBoltDefinedSweep
export MultiBoltModelType, HDModel, HDGEModel, HDGE01Model, SSTModel
export MultiBoltSweepVariable, ENTdSweep, TKSweep, pTorrSweep, BinFracSweep, NuSweep, NTermsSweep
export MultiBoltScatteringModel, IsotropicScattering, IdealForwardScattering, ScreenedCoulombScattering
export MultiBoltInterpMethod, LinearInterp, LogarithmicInterp

include("solvers/multibolt_run.jl")
export MultiBoltRunResult


