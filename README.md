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
- `MultiBolt()`
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

## TODO

- Add test data for `MultiBolt` solver output (e.g. under `test/data/`, following the same pattern as `test/data/loki/` and `test/data/bolsig/`) and real tests in `test/runtests.jl` — `LoKI` and `BOLSIG` now have real tests against real solver output, but `MultiBolt`'s `load_raw_dataframe`/`default_swarm_names` remain untested.
