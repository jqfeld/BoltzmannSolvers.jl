# Running the real MultiBolt binary from a MultiBoltInput. Behavior below
# was determined by running the actual `multibolt_linux` binary directly
# (its own `--help` and `bin/script_examples/sh_script/*.sh`), not from
# documentation alone:
#
#   - Configuration is entirely CLI arguments (see multibolt_input.jl) — no
#     input file, so unlike BOLSIG+ there's no file-format quirk to work
#     around: both `--LXCat_Xsec_fid` and `--export_location` accept
#     absolute paths directly (tested), so no collision-file
#     copying/symlinking is needed the way `run_solver(::BOLSIGInput; ...)`
#     needs it.
#   - MultiBolt exits with code 0 on a normal completed run (tested) — no
#     analog of BOLSIG+'s "always exits 2, that's normal" quirk.

"""
    MultiBoltRunResult

Result of [`run_solver`](@ref)`(::MultiBoltInput; ...)`.

- `workdir` — the temporary directory used as `--export_location`.
- `output_dir` — `joinpath(workdir, config.export_name)`, where MultiBolt
  actually wrote its output — pass this directly to
  `load_dataframe(MultiBolt(), ...)`.
- `success` — `exit_code == 0 && isdir(output_dir)`.
- `exit_code`, `log` — the process's exit code and combined stdout/stderr,
  for diagnostics when `success` is `false`.
"""
struct MultiBoltRunResult
    workdir::String
    output_dir::String
    success::Bool
    exit_code::Int
    log::String
end

_resolve_multibolt_path(path::AbstractString) = path
function _resolve_multibolt_path(::Nothing)
    haskey(ENV, "MULTIBOLT_PATH") && return ENV["MULTIBOLT_PATH"]
    error(
        "run_solver: no MultiBolt binary path given — pass `multibolt_path=\"/path/to/multibolt_linux\"` " *
        "or set the MULTIBOLT_PATH environment variable."
    )
end

"""
    run_solver(config::MultiBoltInput; multibolt_path=nothing) -> MultiBoltRunResult

Run MultiBolt with `config`, exporting to a fresh temporary directory.

`multibolt_path` — path to the `multibolt_linux`/`multibolt_win64.exe`
executable. Resolved from this keyword argument first, then the
`MULTIBOLT_PATH` environment variable; errors clearly if neither is given.

# Example

```julia
result = run_solver(config; multibolt_path="/opt/multibolt/multibolt_linux")
result.success || error("MultiBolt run failed:\\n\$(result.log)")
df = load_dataframe(MultiBolt(), result.output_dir)
```
"""
function run_solver(config::MultiBoltInput; multibolt_path::Union{AbstractString,Nothing}=nothing)
    path = abspath(_resolve_multibolt_path(multibolt_path))
    workdir = mktempdir()
    args = _multibolt_args(config, workdir)

    io = IOBuffer()
    exit_code = try
        proc = run(pipeline(Cmd([path; args]); stdin=devnull, stdout=io, stderr=io); wait=true)
        proc.exitcode
    catch e
        e isa Base.ProcessFailedException ? only(e.procs).exitcode : rethrow()
    end

    output_dir = joinpath(workdir, config.export_name)
    success = exit_code == 0 && isdir(output_dir)
    return MultiBoltRunResult(workdir, output_dir, success, exit_code, String(take!(io)))
end
