# BoltzmannSolvers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jqfeld.github.io/BoltzmannSolvers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jqfeld.github.io/BoltzmannSolvers.jl/dev/)
[![Build Status](https://github.com/jqfeld/BoltzmannSolvers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jqfeld/BoltzmannSolvers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jqfeld/BoltzmannSolvers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jqfeld/BoltzmannSolvers.jl)

WIP package for working with electron Boltzmann-equation solver codes in Julia: loading their swarm-parameter/rate-coefficient output into a normalized `DataFrame`, and — for BOLSIG+ and MultiBolt — configuring and running the solvers themselves.

## Installation

```julia
]add https://github.com/jqfeld/BoltzmannSolvers.jl.git
```

## Solvers

Three solvers are supported. Reading output always goes through the same entry point, `load_dataframe(solver, source; replacements, normalize, kwargs...)`, which normalizes each solver's native column names to a shared set of symbols (e.g. `:reduced_field`, `:mean_energy`, `:reduced_mobility`) so results from different solvers can be compared directly.

### `LoKI()`

`source` is either a directory of LoKI-B's plain-text lookup tables (`lookUpTableSwarm.txt`/`lookUpTableRateCoeff.txt`, plus `lookUpTablePower.txt` — per-channel electron energy gain/loss balance — when present, all joined on `RedField(Td)`), or a single `.h5` file, LoKI-B's HDF5 export of the same data.

Reading `.h5` files requires `using HDF5` first — `HDF5.jl` is a weak dependency (`ext/BoltzmannSolversHDF5Ext.jl`), since it pulls in a heavy transitive dependency tree most plain-text users don't need. `JLD2.jl` is *not* a substitute here: it's Julia's own serialization format, not a general HDF5 reader, and errors on this file's `reducedField` dataset.

### `MultiBolt()`

`source` is a run's export directory (`muN_FLUX.txt`/`DTN_FLUX.txt`/`DLN_FLUX.txt`/`avg_en.txt`/`alpha_eff_N.txt` joined on `E_N`, plus a `PerGas/<gas>/*` directory of per-reaction-rate files).

MultiBolt has no input *file* format — it's entirely command-line-argument driven, so there's no `read_multibolt_input`/`write_multibolt_input` pair the way BOLSIG+ has:

- `MultiBoltInput(; kwargs...)` models one run's configuration (cross-section files, species/fractions, sweep, scattering models, convergence settings, ...).
- `run_multibolt(config::MultiBoltInput; multibolt_path=nothing)` translates it into CLI args and runs the binary in a fresh temp directory. `multibolt_path` resolves from the keyword argument, then the `MULTIBOLT_PATH` environment variable. Returns a `MultiBoltRunResult` (`workdir`, `output_dir`, `success`, `exit_code`, `log`), with `success = exit_code == 0 && isdir(output_dir)` — MultiBolt exits `0` on a normal run.

**Caveat**: a bare invocation with no `--sweep_option`/`--sweep_style` at all reliably crashes `multibolt_linux` with `std::bad_alloc`, regardless of grid/convergence settings — a real bug in the binary itself. `MultiBoltInput`'s `sweep = nothing` default ("just run once at the given conditions") is handled by transparently routing it through an equivalent 1-point `def` sweep at the current `EN_Td`, so this doesn't surface as a footgun.

### `BOLSIG()`

`source` is a single BOLSIG+ output file (not a directory). Auto-detects which of BOLSIG+'s three result layouts the file uses: the condensed `R#`/`A#`/`C#` indexed-table format, the verbose block-per-quantity format (values wrapped across multiple lines), or the single-run `R#`-header name/value report. 2D parametric scans (varying two conditions at once) are not yet supported.

BOLSIG+ *input* scripts (the `.dat` files that configure a run) can also be read/written, independent of the output-reading interface:

- `read_bolsig_input(path) -> BOLSIGInput` / `write_bolsig_input(path, input::BOLSIGInput)`. `BOLSIGInput` models one standalone run — collision data sources, the `CONDITIONS` block (any field may be `VAR`, BOLSIG+'s placeholder for a swept variable), one or more run specs (`BOLSIGFixedRun`/`BOLSIGExplicitRun`/`BOLSIGSeriesRun`/`BOLSIGRun2D`), and a `BOLSIGSaveResults`. `BOLSIGSaveResults.format` predicts which output layout you'll get back: `1`→single-run report, `2`→condensed, `3`→verbose (`4`/`5`/`6` — Energy/SIGLO/PLASIMO — aren't read by `BOLSIG()`). Doesn't model BOLSIG+'s full scripting language (e.g. concatenating multiple scripts via `CLEARCOLLISIONS`/`CLEARRUNS` resets in one file) — one `BOLSIGInput` is one script.
- `run_bolsig(input::BOLSIGInput; bolsig_path=nothing, collision_dir=pwd())` runs the real `bolsigminus` binary: writes `input` to a fresh temp directory, places a copy/symlink of every referenced collision file there (BOLSIG+'s input format can't reference them by absolute or subdirectory-relative path — its own `/`-comment convention breaks both), and runs the binary there. `bolsig_path` resolves from the keyword argument, then the `BOLSIGMINUS_PATH` environment variable. Returns a `BOLSIGRunResult` (`workdir`, `output_files`, `success`, `exit_code`, `log`).

**Caveat**: `success` is based on whether the expected output file(s) exist, not on `exit_code` — BOLSIG+ always tries to read further commands from stdin after finishing its input file and exits with code `2` when there's nothing more to read, on every run, success or failure alike. That's normal behavior, not something to work around.

## Examples

```julia
julia> using BoltzmannSolvers
julia> df = load_dataframe(MultiBolt(), "/some/path/to/solver/output", replacements=["VIBV1" => "v=1"])
julia> meanE = create_interpolation(df, :mean_energy, :reduced_field)
julia> meanE(100) # value of :mean_energy at 100 Td
```

Reading a LoKI-B HDF5 export:

```julia
julia> using BoltzmannSolvers, HDF5
julia> df = load_dataframe(LoKI(), "/some/path/to/simulation.h5")
```

Generating and running a BOLSIG+ input script (needs the real `bolsigminus` binary and cross-section database on disk — set `BOLSIGMINUS_PATH` or pass `bolsig_path`):

```julia
julia> using BoltzmannSolvers
julia> input = BOLSIGInput(;
           collisions = [BOLSIGReadCollisions("SigloDataBase-LXCat-04Jun2013.txt", ["Ar"], true)],
           # BOLSIGSeriesRun identifies its variable by code, not a VAR marker,
           # so reduced_field here is just a placeholder (overridden by the scan).
           conditions = BOLSIGConditions(; reduced_field=1.0, gas_fractions=[1.0]),
           runs       = [BOLSIGSeriesRun([BOLSIGSeriesSegment(ReducedFieldVar, 1.0, 1000.0, 50, ExponentialSeries)])],
           save       = BOLSIGSaveResults(file="results.dat", format=3),
       )
julia> result = run_bolsig(input; bolsig_path="/opt/bolsig/bolsigminus", collision_dir="/opt/bolsig/xsdata")
julia> result.success || error("BOLSIG+ run failed:\n$(result.log)")
julia> df = load_dataframe(BOLSIG(), result.output_files[1])
```

Configuring and running MultiBolt (needs the real `multibolt_linux`/`multibolt_win64.exe` binary on disk — set `MULTIBOLT_PATH` or pass `multibolt_path`):

```julia
julia> using BoltzmannSolvers
julia> config = MultiBoltInput(;
           cross_section_files = ["/opt/multibolt/cross-sections/Biagi_N2.txt", "/opt/multibolt/cross-sections/Biagi_Ar.txt"],
           species             = [MultiBoltSpecies("N2", 0.5), MultiBoltSpecies("Ar", 0.5)],
           export_name         = "my_run",
           sweep               = MultiBoltSweep(ENTdSweep, MultiBoltDefinedSweep([50.0, 100.0, 200.0])),
       )
julia> result = run_multibolt(config; multibolt_path="/opt/multibolt/multibolt_linux")
julia> result.success || error("MultiBolt run failed:\n$(result.log)")
julia> df = load_dataframe(MultiBolt(), result.output_dir)
```

Runnable scripts (`julia --project=examples examples/<script>.jl`) are in `examples/`:

- `first_run.jl` — builds and runs a single BOLSIG+ E/N scan for an Ar/He mixture.
- `gas_fraction_scan.jl` — repeats that scan across four Ar/He ratios, then fits a fixed-hyperparameter Gaussian process (over reduced field × gas fraction) to interpolate between them, with an uncertainty-before/after-conditioning demo.
- `argon_solver_comparison.jl` — runs the same pure-Ar E/N scan through both BOLSIG+ and MultiBolt (sharing one cross-section file between them) and compares the results.

## TODO

- `MultiBolt`/`BOLSIG` real-run tests (`test/test_multibolt_run.jl`/`test/test_bolsig_run.jl`) are skip-if-absent — they only run when the real (uncommitted, platform-specific) solver binaries are available locally under `_research/`, so CI never actually exercises them.
