using BoltzmannSolvers
using DataFrames
using Printf

# Pure Ar plasma, E/N swept over the same 7 points with both BOLSIG+ and
# MultiBolt, using the *same* cross-section file for both. Biagi_Ar.txt
# (examples/data/, copied from _research/MultiBolt/cross-sections/) turns
# out to be plain LXCat format — the same block structure BOLSIG+'s own
# SigloDataBase-LXCat-04Jun2013.txt uses (confirmed by running bolsigminus
# against it directly) — so sharing one file here isolates the comparison to
# the two solvers' numerical methods (BOLSIG+'s two-term expansion vs.
# MultiBolt's default HD+GE multi-term expansion), rather than also mixing
# in a difference in the underlying collision data.
en_td_values = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]
collision_dir = joinpath(@__DIR__, "data")
xsec_file = joinpath(collision_dir, "Biagi_Ar.txt")

# ── BOLSIG+ ──────────────────────────────────────────────────────────────
# BOLSIGExplicitRun is BOLSIG+'s arbitrary-list scan (as opposed to
# BOLSIGSeriesRun's formula-based ranges) — the only way to hit exactly the
# same E/N points MultiBolt's `def` sweep below uses.
bolsig_input = BOLSIGInput(;
    collisions = [BOLSIGReadCollisions("Biagi_Ar.txt", ["Ar"], true)],
    conditions = BOLSIGConditions(; reduced_field=nothing, gas_fractions=[1.0]),
    runs       = [BOLSIGExplicitRun(reshape(en_td_values, :, 1))],
    save       = BOLSIGSaveResults(file="ar_scan.dat", format=3),
)
bolsig_result = run_solver(bolsig_input; collision_dir)
bolsig_result.success || error("BOLSIG+ run failed:\n$(bolsig_result.log)")
bolsig_df = load_dataframe(BOLSIG(), bolsig_result.output_files[1])

# ── MultiBolt ────────────────────────────────────────────────────────────
# MultiBoltDefinedSweep is the direct analog of BOLSIGExplicitRun above.
multibolt_config = MultiBoltInput(;
    cross_section_files = [xsec_file],
    species             = [MultiBoltSpecies("Ar", 1.0)],
    export_name         = "ar_scan",
    sweep               = MultiBoltSweep(ENTdSweep, MultiBoltDefinedSweep(en_td_values)),
)
multibolt_result = run_solver(multibolt_config; multibolt_path=get(ENV, "MULTIBOLT_PATH", nothing))
multibolt_result.success || error("MultiBolt run failed:\n$(multibolt_result.log)")
multibolt_df = load_dataframe(MultiBolt(), multibolt_result.output_dir)

# ── Compare ──────────────────────────────────────────────────────────────
# Both solvers write rows in the order the E/N values were given, and both
# normalize their native column names via `default_swarm_names` — but only
# :reduced_field, :mean_energy, :reduced_mobility, and
# :reduced_townsend_alpha_coef end up named identically on both sides.
# Diffusion is not among them: BOLSIG+'s two-term model produces a single
# isotropic `:reduced_diffusion_coef`, while MultiBolt's multi-term model
# splits it into separate `:longitudinal_diff_coef`/`:transversal_diff_coef`
# — different physical quantities, not just different names, so they're
# left out of this comparison rather than forced together.
compare_columns = [:mean_energy, :reduced_mobility, :reduced_townsend_alpha_coef]

@assert bolsig_df.reduced_field == multibolt_df.reduced_field == en_td_values

# Below the ~11.5 eV excitation / 15.76 eV ionization thresholds, alpha is
# essentially zero on both sides — a relative difference there is dividing
# noise by noise, not a meaningful number, so it's shown as "-" instead.
println("E/N scan comparison (pure Ar, shared Biagi_Ar.txt cross sections):\n")
for col in compare_columns
    println(col, ":")
    @printf("  %-10s %-14s %-14s %s\n", "E/N (Td)", "BOLSIG+", "MultiBolt", "rel. diff")
    for (i, en) in enumerate(en_td_values)
        b, m = bolsig_df[i, col], multibolt_df[i, col]
        rel_diff_str = isnan(m) || max(abs(b), abs(m)) < 1e-25 ? "-" : @sprintf("%.1f%%", 100 * (m - b) / b)
        @printf("  %-10g %-14.4g %-14.4g %s\n", en, b, m, rel_diff_str)
    end
    println()
end
