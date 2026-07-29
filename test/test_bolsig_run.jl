using BoltzmannSolvers
using Test

const INPUT_DATA = joinpath(@__DIR__, "data", "bolsig_input")
const REFERENCE_DATA = joinpath(@__DIR__, "data", "bolsig")

# The real bolsigminus binary (27 MB, platform-specific) and its
# cross-section database aren't committed to the repo — only present when
# working from this project's own `_research/` scratch directory. Skip
# (with a warning) rather than fail when they're not available, same
# pattern as the HDF5-not-loaded guard in test_loki.jl.
_default_bolsig_path() = joinpath(@__DIR__, "..", "_research", "bin", "bolsigminus")
_default_collision_dir() = joinpath(@__DIR__, "..", "_research", "bolsig")

const BOLSIG_PATH = get(ENV, "BOLSIGMINUS_PATH", _default_bolsig_path())
const COLLISION_DIR = _default_collision_dir()

if isfile(BOLSIG_PATH) && isdir(COLLISION_DIR)
    @testset "run_bolsig reproduces the reference outputs exactly" begin
        for (input_name, reference_name) in (
            ("example1_input.dat", "example1.dat"),
            ("example4_input.dat", "example4.dat"),
        )
            input = read_bolsig_input(joinpath(INPUT_DATA, input_name))
            result = run_bolsig(input; bolsig_path=BOLSIG_PATH, collision_dir=COLLISION_DIR)

            @test result.success
            @test length(result.output_files) == 1
            @test isfile(result.output_files[1])
            @test read(result.output_files[1], String) == read(joinpath(REFERENCE_DATA, reference_name), String)

            # exit code 2 is BOLSIG+'s *normal* termination (see module notes
            # in src/solvers/bolsig_run.jl) — assert it's specifically 2,
            # not just "nonzero", so a real regression (e.g. a different
            # failure exit code) still gets caught.
            @test result.exit_code == 2
        end
    end

    @testset "run_bolsig error paths" begin
        input = read_bolsig_input(joinpath(INPUT_DATA, "example1_input.dat"))

        @test_throws ErrorException run_bolsig(input; bolsig_path=BOLSIG_PATH, collision_dir=tempdir())

        withenv("BOLSIGMINUS_PATH" => nothing) do
            @test_throws ErrorException run_bolsig(input; collision_dir=COLLISION_DIR)
        end

        withenv("BOLSIGMINUS_PATH" => BOLSIG_PATH) do
            result = run_bolsig(input; collision_dir=COLLISION_DIR)
            @test result.success
        end
    end
else
    @warn "Skipping run_bolsig tests — bolsigminus binary/database not found (expected outside this project's own _research/ checkout)" BOLSIG_PATH COLLISION_DIR
end
