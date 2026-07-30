using BoltzmannSolvers
using DataFrames
using Test

# The real multibolt_linux binary and its bundled cross-section files aren't
# committed to the repo — only present when working from this project's own
# `_research/` scratch directory. Skip (with a warning) rather than fail
# when they're not available, same pattern as test_bolsig_run.jl.
_default_multibolt_path() = joinpath(@__DIR__, "..", "_research", "MultiBolt", "bin", "multibolt_linux")
_default_xsec_dir() = joinpath(@__DIR__, "..", "_research", "MultiBolt", "cross-sections")

const MULTIBOLT_PATH = get(ENV, "MULTIBOLT_PATH", _default_multibolt_path())
const XSEC_DIR = _default_xsec_dir()

# Deliberately small/fast settings (few grid points, loose convergence) —
# this test only needs to confirm the plumbing (config -> CLI args -> real
# run -> readable output), not physically accurate results.
_fast_settings(; kwargs...) = MultiBoltInput(;
    N_terms=2, Nu=50, initial_eV_max=50.0, conv_err=1e-3, iter_max=20, iter_min=4,
    kwargs...,
)

if isfile(MULTIBOLT_PATH) && isdir(XSEC_DIR)
    @testset "run_multibolt: single-species run is readable" begin
        config = _fast_settings(;
            cross_section_files=[joinpath(XSEC_DIR, "Biagi_Ar.txt")],
            species=[MultiBoltSpecies("Ar", 1.0)],
            export_name="test_single",
            sweep=MultiBoltSweep(ENTdSweep, MultiBoltDefinedSweep([100.0])),
        )
        result = run_multibolt(config; multibolt_path=MULTIBOLT_PATH)

        @test result.success
        @test result.exit_code == 0
        @test isdir(result.output_dir)

        df = load_dataframe(MultiBolt(), result.output_dir)
        @test size(df, 1) == 1
        @test df.reduced_field[1] ≈ 100.0
        @test "mean_energy" in names(df)
    end

    @testset "run_multibolt: two-species defined sweep is readable" begin
        config = _fast_settings(;
            cross_section_files=[joinpath(XSEC_DIR, "Biagi_N2.txt"), joinpath(XSEC_DIR, "Biagi_Ar.txt")],
            species=[MultiBoltSpecies("N2", 0.5), MultiBoltSpecies("Ar", 0.5)],
            export_name="test_sweep",
            sweep=MultiBoltSweep(ENTdSweep, MultiBoltDefinedSweep([50.0, 100.0, 200.0])),
        )
        result = run_multibolt(config; multibolt_path=MULTIBOLT_PATH)

        @test result.success
        df = load_dataframe(MultiBolt(), result.output_dir)
        @test size(df, 1) == 3
        @test df.reduced_field == [50.0, 100.0, 200.0]
        # reaction columns from both gases should be present
        @test any(n -> occursin("N2", n), names(df))
        @test any(n -> occursin("Ar", n), names(df))
    end

    @testset "run_multibolt error paths" begin
        # bin_frac sweep requires exactly 2 species
        bad_config = _fast_settings(;
            cross_section_files=[joinpath(XSEC_DIR, "Biagi_N2.txt")],
            species=[MultiBoltSpecies("N2", 1.0)],
            export_name="bad",
            sweep=MultiBoltSweep(BinFracSweep, MultiBoltRegularSweep(0.0, 0.2, 1.0)),
        )
        @test_throws ErrorException run_multibolt(bad_config; multibolt_path=MULTIBOLT_PATH)

        good_config = _fast_settings(;
            cross_section_files=[joinpath(XSEC_DIR, "Biagi_Ar.txt")],
            species=[MultiBoltSpecies("Ar", 1.0)],
            export_name="test_no_path",
        )
        withenv("MULTIBOLT_PATH" => nothing) do
            @test_throws ErrorException run_multibolt(good_config)
        end

        withenv("MULTIBOLT_PATH" => MULTIBOLT_PATH) do
            result = run_multibolt(good_config)
            @test result.success
        end
    end
else
    @warn "Skipping run_multibolt tests — multibolt_linux binary/cross-sections not found (expected outside this project's own _research/ checkout)" MULTIBOLT_PATH XSEC_DIR
end
