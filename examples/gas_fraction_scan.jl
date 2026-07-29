using BoltzmannSolvers
using DataFrames

# Same run as first_run.jl (10%/90% Ar/He mixture, E/N swept 0.1-1000 Td via
# two chained RUNSERIES segments), repeated for four different Ar/He mixture
# ratios — each run gets its own fresh temp directory from run_bolsig, so
# reusing the same output filename across runs is safe (no collisions).
function build_input(gas_fractions)
    return BOLSIGInput(;
        collisions = [
            BOLSIGReadCollisions("SigloDataBase-LXCat-04Jun2013.txt", ["Ar"], true),
            BOLSIGReadCollisions("SigloDataBase-LXCat-04Jun2013.txt", ["He"], true),
        ],
        conditions = BOLSIGConditions(;
            reduced_field = 10.0,   # placeholder — BOLSIGSeriesRun selects its
                                    # variable by code, not a VAR marker, so this
                                    # value is overridden by the scan below
            n_grid_points = 400,
            gas_fractions = gas_fractions,
        ),
        runs = [
            BOLSIGSeriesRun([
                BOLSIGSeriesSegment(ReducedFieldVar, 0.1, 100.0, 31, ExponentialSeries),
                BOLSIGSeriesSegment(ReducedFieldVar, 200.0, 1000.0, 9, LinearSeries),
            ]),
        ],
        save = BOLSIGSaveResults(file="scan.dat", format=3),
    )
end

gas_fraction_sets = [[0.1, 0.9], [0.15, 0.85], [0.2, 0.8], [0.25, 0.75]]
collision_dir = joinpath(@__DIR__, "data")

dataframes = map(gas_fraction_sets) do gas_fractions
    input = build_input(gas_fractions)
    result = run_bolsig(input; collision_dir)
    result.success || error("BOLSIG+ run failed for gas_fractions=$gas_fractions:\n$(result.log)")
    return load_dataframe(BOLSIG(), result.output_files[1])
end

for (gas_fractions, df) in zip(gas_fraction_sets, dataframes)
    println("Ar/He = $(gas_fractions): size(df) = $(size(df)), mean_energy at E/N=10Td ≈ $(df.mean_energy[argmin(abs.(df.reduced_field .- 10.0))])")
end
