# BoltzmannSolvers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jqfeld.github.io/BoltzmannSolvers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jqfeld.github.io/BoltzmannSolvers.jl/dev/)
[![Build Status](https://github.com/jqfeld/BoltzmannSolvers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jqfeld/BoltzmannSolvers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jqfeld/BoltzmannSolvers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jqfeld/BoltzmannSolvers.jl)

WIP package to use the output of various Boltzmann solver codes in Julia. 

Currently, reading the output from three solvers are implemented:

- `LoKI()` — `source` is either a directory of LoKI-B's plain-text lookup
  tables (`lookUpTableSwarm.txt`/`lookUpTableRateCoeff.txt`, plus
  `lookUpTablePower.txt` — per-channel electron energy gain/loss balance —
  when present, all joined on `RedField(Td)`), or a single `.h5` file —
  LoKI-B's HDF5 export of the same data. Reading `.h5` files requires
  `using HDF5` first (`HDF5.jl` is a weak dependency, loaded via
  `ext/BoltzmannSolversHDF5Ext.jl`, since it pulls in a heavy transitive
  dependency tree that most users reading plain-text output don't need).
  Note: `JLD2.jl` does *not* work for this — it's Julia's own serialization
  format, not a general third-party-HDF5 reader, and errors on this file's
  `reducedField` dataset specifically.
- `MultiBolt()` — `source` is a run's export directory (`muN_FLUX.txt`/
  `DTN_FLUX.txt`/`DLN_FLUX.txt`/`avg_en.txt`/`alpha_eff_N.txt` joined on
  `E_N`, plus a `PerGas/<gas>/*` directory of per-reaction-rate files).

  MultiBolt has no input *file* format at all (unlike BOLSIG+) — it's
  entirely command-line-argument driven (confirmed by running the real
  `multibolt_linux --help` and its own bundled example scripts), so there's
  no `read_multibolt_input` counterpart to `read_bolsig_input`:
  `MultiBoltInput(; kwargs...)` models one run's configuration, and
  `run_multibolt(config::MultiBoltInput; multibolt_path=nothing)` runs it
  directly (translating the struct into CLI args, exporting to a fresh temp
  directory) — no intermediate file, no `write_multibolt_input`.
  `multibolt_path` resolves from the keyword argument, then the
  `MULTIBOLT_PATH` environment variable. Returns a `MultiBoltRunResult`
  (`workdir`, `output_dir`, `success`, `exit_code`, `log`) —
  `success = exit_code == 0 && isdir(output_dir)`; unlike BOLSIG+,
  MultiBolt exits `0` on a normal run (confirmed directly), no quirky
  "always nonzero" behavior to work around.

  One real bug in the binary itself, found by testing (not documented
  anywhere): a truly bare run — no `--sweep_option`/`--sweep_style` at all —
  reliably crashes `multibolt_linux` with `std::bad_alloc`, regardless of
  grid/convergence settings. `MultiBoltInput`'s `sweep = nothing` (its
  default, meaning "just run once at the given conditions") is handled by
  transparently routing it through an equivalent 1-point `def` sweep at the
  current `EN_Td` instead, so this doesn't surface as a footgun.
- `BOLSIG()` — `source` is a single BOLSIG+ output file (not a directory).
  Auto-detects which of BOLSIG+'s three result layouts the file uses: the
  condensed `R#`/`A#`/`C#` indexed-table format, the verbose
  block-per-quantity format (values wrapped across multiple lines), or the
  single-run `R#`-header name/value report. 2D parametric scans (varying two
  conditions at once) are not yet supported.

  BOLSIG+ *input* scripts (the `.dat` files that configure a run) can also
  be read/written, independent of the `Solver` output-reading interface:
  `read_bolsig_input(path) -> BOLSIGInput` and
  `write_bolsig_input(path, input::BOLSIGInput)`. `BOLSIGInput` models one
  standalone run — collision data sources, the `CONDITIONS` block (any
  field may be `VAR`, BOLSIG+'s placeholder for a swept variable), one or
  more run specs (`BOLSIGFixedRun`/`BOLSIGExplicitRun`/`BOLSIGSeriesRun`/
  `BOLSIGRun2D`), and a `BOLSIGSaveResults`. `BOLSIGSaveResults.format`
  predicts which output layout above you'll get back: `1`→single-run report,
  `2`→condensed, `3`→verbose (`4`/`5`/`6` — Energy/SIGLO/PLASIMO — aren't
  read by `BOLSIG()`). Doesn't model BOLSIG+'s full scripting language (e.g.
  concatenating multiple scripts with `CLEARCOLLISIONS`/`CLEARRUNS` resets
  in one file) — one `BOLSIGInput` is one script.

  `run_bolsig(input::BOLSIGInput; bolsig_path=nothing, collision_dir=pwd())`
  runs the real `bolsigminus` binary against a `BOLSIGInput`: writes it to a
  fresh temp directory, places a copy/symlink of every referenced collision
  file there (BOLSIG+'s input format can't reference them by absolute or
  subdirectory-relative path — its own `/`-comment convention breaks both),
  and runs the binary there. `bolsig_path` resolves from the keyword
  argument, then the `BOLSIGMINUS_PATH` environment variable. Returns a
  `BOLSIGRunResult` with `workdir`, `output_files`, `success`, `exit_code`,
  and `log`. **`success` is based on whether the expected output file(s)
  exist, not on `exit_code`** — BOLSIG+ always tries to read further
  commands from stdin after finishing its input file and errors out with
  exit code `2` when there's nothing more to read, on every run, success or
  failure alike; that isn't a bug to work around, just how the binary
  behaves when run non-interactively.

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

Generating a BOLSIG+ input script:

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
julia> write_bolsig_input("run.dat", input)
```

Running it (needs the real `bolsigminus` binary and cross-section database
on disk — set `BOLSIGMINUS_PATH` or pass `bolsig_path`):

```julia
julia> result = run_bolsig(input; bolsig_path="/opt/bolsig/bolsigminus", collision_dir="/opt/bolsig/xsdata")
julia> result.success || error("BOLSIG+ run failed:\n$(result.log)")
julia> df = load_dataframe(BOLSIG(), result.output_files[1])
```

Configuring and running MultiBolt (needs the real `multibolt_linux`/
`multibolt_win64.exe` binary on disk — set `MULTIBOLT_PATH` or pass
`multibolt_path`):

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

## TODO

- `MultiBolt`/`BOLSIG` real-run tests (`test/test_multibolt_run.jl`/
  `test/test_bolsig_run.jl`) are skip-if-absent — they only run when the
  real (uncommitted, platform-specific) solver binaries are available
  locally under `_research/`, so CI never actually exercises them.
