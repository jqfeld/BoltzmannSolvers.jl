# BoltzmannSolvers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jqfeld.github.io/BoltzmannSolvers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jqfeld.github.io/BoltzmannSolvers.jl/dev/)
[![Build Status](https://github.com/jqfeld/BoltzmannSolvers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jqfeld/BoltzmannSolvers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jqfeld/BoltzmannSolvers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jqfeld/BoltzmannSolvers.jl)

WIP package to use the output of various Boltzmann solver codes in Julia. 

Currently, reading the output from three solvers are implemented:

- `LoKI()`
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

## TODO

- Add synthetic test data for `LoKI`/`MultiBolt` solver output (e.g. under `test/`, similar to `LXCat.jl`'s `test/test_data.txt`) and real tests in `test/runtests.jl` — currently that file is an empty placeholder, so `load_raw_dataframe`/`default_swarm_names`/`parse_reaction_names` for both solvers are untested.
