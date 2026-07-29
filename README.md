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

## TODO

- Add test data for `MultiBolt` solver output (e.g. under `test/data/`, following the same pattern as `test/data/loki/` and `test/data/bolsig/`) and real tests in `test/runtests.jl` — `LoKI` and `BOLSIG` now have real tests against real solver output, but `MultiBolt`'s `load_raw_dataframe`/`default_swarm_names` remain untested.
