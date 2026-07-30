using BoltzmannSolvers
using DataFrames
using AbstractGPs
using KernelFunctions: RowVecs, SEKernel, with_lengthscale
using Statistics: mean, std

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

# ── Gaussian process regression: (reduced_field, Ar fraction) -> everything else ──
#
# Pool the four runs into one training set. `gas_fractions` is per-run, not
# per-row, so it's broadcast into a new column before stacking. `reduced_field`
# is used in log10 space as a GP input (it spans 0.1-1000 Td — a squared-
# exponential kernel needs a roughly-uniform notion of "distance" between
# points, which log-spacing gives here since the RUNSERIES scan itself is
# log-spaced over most of that range).
for (gas_fractions, df) in zip(gas_fraction_sets, dataframes)
    df.ar_fraction = fill(gas_fractions[1], nrow(df))
end
combined = vcat(dataframes...)

Xmatrix = hcat(log10.(combined.reduced_field), combined.ar_fraction)
Xdata = RowVecs(Xmatrix)
output_columns = setdiff(names(combined), ["reduced_field", "ar_fraction"])

# One independent GP per output column (the standard approach for
# multi-output regression absent a specific reason to model cross-column
# correlations) — squared-exponential kernel with fixed, hand-picked
# lengthscales (one per input dimension), no hyperparameter fitting.
# BOLSIG+ is deterministic (the same input always gives the same output), so
# there's no real observation noise to model either — `JITTER` below is a
# tiny nugget added purely so the covariance matrix's Cholesky factorization
# is numerically well-behaved, not a statistical noise term. Lengthscales
# are chosen to roughly match the spacing of the training inputs: the E/N
# scan is log-spaced from 0.1-1000 Td (log10 spans [-1, 3]), and the four
# Ar-fraction training points are spaced 0.05 apart.
const LENGTHSCALES = [0.5, 0.05]   # [log10(reduced_field), ar_fraction]
const SIGNAL_STD = 1.0             # kernel output scale, in standardized (z-score) units
const JITTER = 1e-6                # numerical-stability nugget, not modeled noise

function fit_gp(X, y::AbstractVector{<:Real})
    σy, μy = std(y), mean(y)
    # Columns that are exactly constant (e.g. an inactive reaction channel,
    # always 0 across every run) carry no signal to condition a GP on.
    (σy == 0 || !isfinite(σy)) && return (; constant=μy)
    yz = (y .- μy) ./ σy

    kernel = SIGNAL_STD^2 * with_lengthscale(SEKernel(), LENGTHSCALES)
    fx = GP(kernel)(X, JITTER)
    return (; post=posterior(fx, yz), μy, σy)
end

# Returns (predicted mean, predicted std), both back in the column's
# original (unstandardized) units. `mean_and_var` gives the GP's predictive
# variance at Xnew — how uncertain the model is, not just its point estimate.
function predict_gp(fit, Xnew)
    haskey(fit, :constant) && return fit.constant, 0.0
    μz, σz2 = mean_and_var(fit.post(Xnew))
    return fit.μy + fit.σy * only(μz), fit.σy * sqrt(only(σz2))
end

println("\nConditioning a GP per output column ($(length(output_columns)) columns)...")
gp_models = Dict(col => fit_gp(Xdata, Float64.(combined[!, col])) for col in output_columns)

# ── Two nearby points: actual vs. predicted vs. GP uncertainty ─────────────
#
# Two Ar/He mixtures close to one another (neither in the training set) at
# the same E/N — comparing both the point prediction *and* the GP's own
# uncertainty at two nearby locations shows the uncertainty behaving as a GP
# should, not just whether a single prediction happens to land close.
demo_columns = ("mean_energy", "reduced_mobility", "reduced_townsend_alpha_coef")

function run_and_extract(gas_fractions; e_n_target=10.0)
    result = run_bolsig(build_input(gas_fractions); collision_dir)
    result.success || error("BOLSIG+ run failed for gas_fractions=$gas_fractions:\n$(result.log)")
    df = load_dataframe(BOLSIG(), result.output_files[1])
    row = argmin(abs.(df.reduced_field .- e_n_target))
    point = RowVecs(reshape([log10(df.reduced_field[row]), gas_fractions[1]], 1, 2))
    return df, row, point
end

fractions_A = [0.18, 0.82]
fractions_B = [0.19, 0.81]
df_A, row_A, point_A = run_and_extract(fractions_A)
df_B, row_B, point_B = run_and_extract(fractions_B)

function compare_points(label)
    println("\n$label")
    for col in demo_columns
        pred_A, std_A = predict_gp(gp_models[col], point_A)
        pred_B, std_B = predict_gp(gp_models[col], point_B)
        println("  $col:")
        println("    Ar=$(fractions_A[1]): actual=$(df_A[row_A, col]), predicted=$pred_A ± $std_A")
        println("    Ar=$(fractions_B[1]): actual=$(df_B[row_B, col]), predicted=$pred_B ± $std_B")
    end
end

compare_points("Before conditioning on either point:")

# Condition on point A's actual data — one real observation added to the
# training set, keeping the *original* standardization (μy, σy) fixed rather
# than recomputed, since it doesn't change just because one more point was
# observed. Point A should now predict ~exactly (query == an observed
# point, up to the numerical JITTER), and point B — close to A but still
# unobserved — should get both a better point estimate and a smaller
# uncertainty purely from proximity to the new observation.
for col in demo_columns
    fit = gp_models[col]
    haskey(fit, :constant) && continue   # nothing to condition on a known constant
    X_ext = RowVecs(vcat(Xmatrix, reshape([log10(df_A.reduced_field[row_A]), fractions_A[1]], 1, 2)))
    y_ext = vcat((combined[!, col] .- fit.μy) ./ fit.σy, [(df_A[row_A, col] - fit.μy) / fit.σy])
    kernel = SIGNAL_STD^2 * with_lengthscale(SEKernel(), LENGTHSCALES)
    fx = GP(kernel)(X_ext, JITTER)
    gp_models[col] = (; post=posterior(fx, y_ext), fit.μy, fit.σy)
end

compare_points("After conditioning on Ar=$(fractions_A[1]):")
