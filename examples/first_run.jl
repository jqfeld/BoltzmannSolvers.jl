using BoltzmannSolvers
using DataFrames

# BOLSIG+'s own first example (_research/input_data/bolsig/input-examples.dat):
# 10%/90% Ar/He mixture, E/N swept 0.1-1000 Td via two chained RUNSERIES
# segments (31 exponentially-spaced points from 0.1-100 Td, then 9
# linearly-spaced from 200-1000 Td — see BOLSIGSeriesRun's docstring for why
# BOLSIG+ chains segments like this instead of one range).
input = BOLSIGInput(;
    collisions = [
        BOLSIGReadCollisions("SigloDataBase-LXCat-04Jun2013.txt", ["Ar"], true),
        BOLSIGReadCollisions("SigloDataBase-LXCat-04Jun2013.txt", ["He"], true),
    ],
    conditions = BOLSIGConditions(;
        reduced_field = 10.0,   # placeholder — BOLSIGSeriesRun selects its
                                # variable by code, not a VAR marker, so this
                                # value is overridden by the scan below
        n_grid_points = 400,
        gas_fractions = [0.1, 0.9],
    ),
    runs = [
        BOLSIGSeriesRun([
            BOLSIGSeriesSegment(ReducedFieldVar, 0.1, 100.0, 31, ExponentialSeries),
            BOLSIGSeriesSegment(ReducedFieldVar, 200.0, 1000.0, 9, LinearSeries),
        ]),
    ],
    save = BOLSIGSaveResults(file="example1.dat", format=3),
)

# bolsig_path resolves from the BOLSIGMINUS_PATH environment variable if not
# given explicitly here — see run_bolsig's docstring (src/solvers/bolsig_run.jl).
result = run_bolsig(input; collision_dir=joinpath(@__DIR__, "data"))
result.success || error("BOLSIG+ run failed:\n$(result.log)")

df = load_dataframe(BOLSIG(), result.output_files[1])
@show size(df)
@show df.reduced_field[[1, end]]
@show df.mean_energy[[1, end]]
