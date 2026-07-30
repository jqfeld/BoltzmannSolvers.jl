# Running the real BOLSIG+ binary from a BOLSIGInput. Behavior below was
# determined by running the actual `bolsigminus` binary directly, not from
# documentation (BOLSIG+'s input/output format docs don't cover process
# invocation at all):
#
#   - The input file path is `argv[1]`.
#   - After finishing the given input file, BOLSIG+ always tries to read
#     further commands from stdin, and — since there's nothing more to
#     read — errors with a Fortran "End of file" runtime error and exits
#     with code 2. This happens on *every* successful run, so exit code 2
#     is normal, not a failure signal; `stdin` must be explicitly connected
#     to `devnull` (not inherited) so this happens immediately rather than
#     the process hanging on a real terminal's stdin.
#   - Because BOLSIG+'s own input format uses "/" as an inline-comment
#     delimiter, a collision-database filename can't be an absolute path
#     (misread as an empty field) or even a subdirectory-relative path
#     (breaks its fixed-format READCOLLISIONS parsing) — only a bare
#     filename resolved against the working directory works. So every
#     referenced collision file must be placed, under its exact bare
#     filename, directly alongside the generated input file.
#   - Success must be judged by whether the expected output file(s) exist,
#     not by exit code.

"""
    BOLSIGRunResult

Result of [`run_solver`](@ref)`(::BOLSIGInput; ...)`.

- `workdir` — the temporary directory the run happened in (`input.dat`,
  symlinked/copied collision files, `bolsiglog.txt`, and any output files
  all live here).
- `output_files` — absolute paths to every file BOLSIG+ was configured to
  write (`input.save.file` if present, plus any `BOLSIGRun2D.output_file`).
- `success` — whether every entry in `output_files` exists after running.
  *Not* derived from `exit_code` (see module notes — exit code 2 is BOLSIG+'s
  normal termination, not a failure indicator).
- `exit_code`, `log` — the process's exit code and combined stdout/stderr,
  for diagnostics when `success` is `false`.
"""
struct BOLSIGRunResult
    workdir::String
    output_files::Vector{String}
    success::Bool
    exit_code::Int
    log::String
end

_resolve_bolsig_path(path::AbstractString) = path
function _resolve_bolsig_path(::Nothing)
    haskey(ENV, "BOLSIGMINUS_PATH") && return ENV["BOLSIGMINUS_PATH"]
    error(
        "run_solver: no BOLSIG+ binary path given — pass `bolsig_path=\"/path/to/bolsigminus\"` " *
        "or set the BOLSIGMINUS_PATH environment variable."
    )
end

"""
    run_solver(input::BOLSIGInput; bolsig_path=nothing, collision_dir=pwd()) -> BOLSIGRunResult

Write `input` to a temporary directory (via [`write_bolsig_input`](@ref),
which also validates it), place every collision file it references
alongside it, and run BOLSIG+ there.

- `bolsig_path` — path to the `bolsigminus` executable. Resolved from this
  keyword argument first, then the `BOLSIGMINUS_PATH` environment variable;
  errors clearly if neither is given.
- `collision_dir` — directory to look up each `BOLSIGReadCollisions.file`
  bare filename in (BOLSIG+'s input format can't reference collision files
  by absolute or subdirectory-relative path — see module notes). Defaults
  to the current working directory.

# Example

```julia
result = run_solver(input; bolsig_path="/opt/bolsig/bolsigminus")
result.success || error("BOLSIG+ run failed:\\n\$(result.log)")
df = load_dataframe(BOLSIG(), result.output_files[1])
```
"""
function run_solver(input::BOLSIGInput; bolsig_path::Union{AbstractString,Nothing}=nothing, collision_dir::AbstractString=pwd())
    # abspath: `Cmd(...; dir=workdir)` below spawns with `workdir` as the
    # working directory, so a relative `bolsig_path` (or a relative
    # `collision_dir`, resolved via `joinpath` further down) must be
    # resolved against *this* process's cwd first, not `workdir`'s.
    path = abspath(_resolve_bolsig_path(bolsig_path))
    collision_dir = abspath(collision_dir)
    workdir = mktempdir()

    # Multiple READCOLLISIONS entries commonly share one database file (e.g.
    # two species read from the same LXCat file) — place each unique bare
    # filename only once.
    for file in unique(c.file for c in input.collisions)
        src = joinpath(collision_dir, file)
        isfile(src) || error(
            "run_solver: collision file '$file' not found in collision_dir=\"$collision_dir\"."
        )
        dst = joinpath(workdir, file)
        try
            symlink(abspath(src), dst)
        catch
            cp(src, dst)   # fall back if symlinks aren't permitted (e.g. some restricted filesystems)
        end
    end

    write_bolsig_input(joinpath(workdir, "input.dat"), input)

    output_files = String[]
    isnothing(input.save) || push!(output_files, joinpath(workdir, input.save.file))
    for run in input.runs
        run isa BOLSIGRun2D && push!(output_files, joinpath(workdir, run.output_file))
    end

    io = IOBuffer()
    exit_code = try
        proc = run(pipeline(Cmd(`$path input.dat`; dir=workdir); stdin=devnull, stdout=io, stderr=io); wait=true)
        proc.exitcode
    catch e
        e isa Base.ProcessFailedException ? only(e.procs).exitcode : rethrow()
    end

    success = all(isfile, output_files)
    return BOLSIGRunResult(workdir, output_files, success, exit_code, String(take!(io)))
end
