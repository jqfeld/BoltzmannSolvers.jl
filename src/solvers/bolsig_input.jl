# Reading/writing BOLSIG+ *input* scripts (as opposed to bolsig.jl, which
# reads BOLSIG+ *output*). Models one coherent, standalone run configuration
# — not the full generality of BOLSIG+'s scripting language (which supports
# concatenating multiple independent scripts in one file via
# CLEARCOLLISIONS/CLEARRUNS resets). See `_research/bolsig/input-examples.dat`
# for the reference format this was built from.

struct BOLSIGReadCollisions
    file::String
    species::Vector{String}
    extrapolate::Bool
end

# Every field is `Union{T,Nothing}` — `nothing` means the line held `VAR`
# (BOLSIG+'s placeholder for "this is a swept variable", resolved by
# whichever `BOLSIGRunSpec` is used — see module notes on `write_bolsig_input`
# for which run kinds actually rely on `VAR` markers). Only `reduced_field`,
# `angular_frequency`, and `gas_fractions` are ever seen as `VAR` in the
# reference examples, but nothing in the format restricts it to those, so
# every field allows it uniformly rather than hardcoding which are "allowed".
Base.@kwdef struct BOLSIGConditions
    reduced_field::Union{Float64,Nothing}   # required — no sensible default for the primary physical condition
    angular_frequency::Union{Float64,Nothing} = 0.0
    eb_angle_cosine::Union{Float64,Nothing} = 0.0
    gas_temperature::Union{Float64,Nothing} = 300.0
    excitation_temperature::Union{Float64,Nothing} = 300.0
    transition_energy::Union{Float64,Nothing} = 0.0
    ionization_degree::Union{Float64,Nothing} = 0.0
    plasma_density::Union{Float64,Nothing} = 1e18
    ion_charge_parameter::Union{Float64,Nothing} = 1.0
    ion_neutral_mass_ratio::Union{Float64,Nothing} = 1.0
    coulomb_collision_model::Union{Int,Nothing} = 3
    energy_sharing::Union{Int,Nothing} = 1
    growth_model::Union{Int,Nothing} = 1
    maxwellian_mean_energy::Union{Float64,Nothing} = 0.0
    n_grid_points::Union{Int,Nothing} = 400
    manual_grid::Union{Int,Nothing} = 0
    manual_max_energy::Union{Float64,Nothing} = 200.0
    precision::Union{Float64,Nothing} = 1e-10
    convergence::Union{Float64,Nothing} = 1e-4
    max_iterations::Union{Int,Nothing} = 1000
    gas_fractions::Vector{Union{Float64,Nothing}}   # length = total species count across all `collisions`
    normalize_fractions::Bool = true
end

# The ordered list of scalar `BOLSIGConditions` fields, matching the fixed
# textual order of BOLSIG+'s CONDITIONS block. Shared between the parser,
# writer, and VAR-position validation so all three stay in sync.
const _BOLSIG_CONDITIONS_SCALAR_FIELDS = (
    :reduced_field, :angular_frequency, :eb_angle_cosine, :gas_temperature,
    :excitation_temperature, :transition_energy, :ionization_degree,
    :plasma_density, :ion_charge_parameter, :ion_neutral_mass_ratio,
    :coulomb_collision_model, :energy_sharing, :growth_model,
    :maxwellian_mean_energy, :n_grid_points, :manual_grid,
    :manual_max_energy, :precision, :convergence, :max_iterations,
)

abstract type BOLSIGRunSpec end

"""
    BOLSIGFixedRun()

A bare `RUN`: a single run at the current `BOLSIGConditions` (no `VAR`
fields need to be resolved).
"""
struct BOLSIGFixedRun <: BOLSIGRunSpec end

"""
    BOLSIGExplicitRun(values)

`RUN` followed by explicit data rows — BOLSIG+'s manual/arbitrary-list scan.
`values` is `n_runs × n_var_fields`; column order matches the order `VAR`
appears in `BOLSIGConditions` (see `_BOLSIG_CONDITIONS_SCALAR_FIELDS`, then
`gas_fractions`). Multiple `BOLSIGExplicitRun`s may appear in one
`BOLSIGInput`'s `runs` (BOLSIG+ allows repeated `RUN` blocks that all
enqueue rows in sequence).
"""
struct BOLSIGExplicitRun <: BOLSIGRunSpec
    values::Matrix{Float64}
end

@enum BOLSIGVariable ReducedFieldVar=1 MeanEnergyVar=2 MaxwellianEnergyVar=3
@enum BOLSIGSeriesType LinearSeries=1 QuadraticSeries=2 ExponentialSeries=3

"""
    BOLSIGSeriesSegment(variable, min, max, count, type)

One `RUNSERIES` block: a formula-based range scan over `variable` (identified
by code, not by a `VAR` marker in `BOLSIGConditions`).
"""
struct BOLSIGSeriesSegment
    variable::BOLSIGVariable
    min::Float64
    max::Float64
    count::Int
    type::BOLSIGSeriesType
end

"""
    BOLSIGSeriesRun(segments)

One or more chained `RUNSERIES` blocks — BOLSIG+ concatenates consecutive
`RUNSERIES` blocks targeting the same or different ranges into one combined
scan (e.g. the reference `example1.dat`'s 40 points come from two segments,
31 exponentially-spaced + 9 linearly-spaced).
"""
struct BOLSIGSeriesRun <: BOLSIGRunSpec
    segments::Vector{BOLSIGSeriesSegment}
end

"""
    BOLSIGRun2D(variable1, min1, max1, num1, type1, min2, max2, num2, type2, output_file)

A `RUN2D` block: a 2D parametric scan. The first variable is identified by
code (like `BOLSIGSeriesSegment`); the second is whichever
`BOLSIGConditions` field is marked `VAR` (there's no numeric code for
arbitrary second dimensions, e.g. a gas-mixture fraction). Writes its own
`output_file` directly, bypassing `BOLSIGSaveResults`.

Note: `BOLSIG()` (`src/solvers/bolsig.jl`) does not read RUN2D's 2D output
format — this type only supports *generating* the input, not reading the
result back.
"""
struct BOLSIGRun2D <: BOLSIGRunSpec
    variable1::BOLSIGVariable
    min1::Float64
    max1::Float64
    num1::Int
    type1::BOLSIGSeriesType
    min2::Float64
    max2::Float64
    num2::Int
    type2::BOLSIGSeriesType
    output_file::String
end

"""
    BOLSIGSaveResults(; file, format, kwargs...)

`format` is BOLSIG+'s output-layout code (1-6). The three layouts
[`load_dataframe`](@ref)`(`[`BOLSIG`](@ref)`(), ...)` actually reads:
`1` (Run by run) → the single-run report format, `2` (Combined) → the
condensed `R#`/`A#`/`C#` table format, `3` (E/N) → the verbose
block-per-quantity format. `4` (Energy), `5` (SIGLO), and `6` (PLASIMO) are
not read by `BOLSIG()`.
"""
Base.@kwdef struct BOLSIGSaveResults
    file::String
    format::Int
    conditions::Bool = true
    transport_coefficients::Bool = true
    rate_coefficients::Bool = true
    reverse_rate_coefficients::Bool = false
    energy_loss_coefficients::Bool = false
    distribution_function::Bool = true
    skip_failed_runs::Bool = false
    include_cross_sections::Bool = true
end

"""
    BOLSIGInput(; noscreen=true, nologfile=true, collisions, conditions, runs, save=nothing)

Configuration for one standalone BOLSIG+ run script. `save` is `nothing`
only when every entry in `runs` is a [`BOLSIGRun2D`](@ref) (which writes its
own output file and doesn't need a `SAVERESULTS` block).

See [`read_bolsig_input`](@ref)/[`write_bolsig_input`](@ref) to parse/emit
the actual `.dat` script format.
"""
Base.@kwdef struct BOLSIGInput
    noscreen::Bool = true
    nologfile::Bool = true
    collisions::Vector{BOLSIGReadCollisions}
    conditions::BOLSIGConditions
    runs::Vector{BOLSIGRunSpec}
    save::Union{BOLSIGSaveResults,Nothing} = nothing
end

# ── Parsing ──────────────────────────────────────────────────────────────────

const _BOLSIG_INPUT_KEYWORDS = Set([
    "READCOLLISIONS", "CLEARCOLLISIONS", "CONDITIONS", "CLEARRUNS", "RUN",
    "RUNSERIES", "RUN2D", "SAVERESULTS", "END", "/NOSCREEN", "/NOLOGFILE",
])

# The "value part" of a line is everything before its first "/" (which may
# itself be followed by more "/" characters as part of a descriptive
# comment, e.g. "10.  / Electric field / N (Td)" — hence `limit=2`).
_bolsig_input_value_part(line::AbstractString) = strip(first(split(line, "/"; limit=2)))

# Flattens a raw input file into a sequence of "meaningful" tokens: keyword
# lines (kept verbatim) and the value-part of data lines. Blank lines, "!"
# comments, and "/" comments are dropped — except the two bare "/NOSCREEN"/
# "/NOLOGFILE" directives, which are themselves keywords despite the
# leading "/".
function _bolsig_input_tokens(path::AbstractString)
    result = String[]
    for line in readlines(path)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "!") && continue
        if startswith(s, "/")
            (s == "/NOSCREEN" || s == "/NOLOGFILE") && push!(result, s)
            continue
        end
        push!(result, _bolsig_input_value_part(s))
    end
    return result
end

_bolsig_or_var(s::AbstractString) = s == "VAR" ? nothing : s
_bolsig_parse_f(s) = s === nothing ? nothing : parse(Float64, s)
_bolsig_parse_i(s) = s === nothing ? nothing : parse(Int, s)

# Parses the 22-field CONDITIONS block starting at `toks[i]`. `n_species` is
# the total species count across every `READCOLLISIONS` seen so far, needed
# to know how many gas-fraction values to consume (they may be split across
# several lines via a trailing "&" continuation marker).
function _parse_bolsig_conditions(toks, i::Int, n_species::Int)
    start = i
    values = Dict{Symbol,Any}()
    for f in _BOLSIG_CONDITIONS_SCALAR_FIELDS
        raw = _bolsig_or_var(toks[i])
        is_int_field = f in (:coulomb_collision_model, :energy_sharing, :growth_model,
                             :n_grid_points, :manual_grid, :max_iterations)
        values[f] = is_int_field ? _bolsig_parse_i(raw) : _bolsig_parse_f(raw)
        i += 1
    end

    fractions = Union{Float64,Nothing}[]
    while length(fractions) < n_species
        line = toks[i]
        i += 1
        continued = endswith(line, "&")
        line = continued ? strip(line[1:end-1]) : line
        for tok in split(line)
            push!(fractions, tok == "VAR" ? nothing : parse(Float64, tok))
        end
    end

    normalize_fractions = parse(Int, toks[i]) == 1
    i += 1

    conditions = BOLSIGConditions(; values..., gas_fractions=fractions, normalize_fractions)
    return conditions, i - start
end

function _parse_bolsig_saveresults(toks, i::Int)
    start = i
    file = toks[i]; i += 1
    format = parse(Int, toks[i]); i += 1
    flags = ntuple(_ -> (v = parse(Int, toks[i]) == 1; i += 1; v), 8)
    save = BOLSIGSaveResults(;
        file, format,
        conditions=flags[1], transport_coefficients=flags[2],
        rate_coefficients=flags[3], reverse_rate_coefficients=flags[4],
        energy_loss_coefficients=flags[5], distribution_function=flags[6],
        skip_failed_runs=flags[7], include_cross_sections=flags[8],
    )
    return save, i - start
end

"""
    read_bolsig_input(path) -> BOLSIGInput

Parse a BOLSIG+ input script into a [`BOLSIGInput`](@ref). Expects one
standalone script (see `BOLSIGInput`'s docstring) — a file that
concatenates multiple scripts via `CLEARCOLLISIONS`/`CLEARRUNS` resets
(like `_research/bolsig/input-examples.dat`) should be split into separate
files first.
"""
function read_bolsig_input(path::AbstractString)
    toks = _bolsig_input_tokens(path)
    n = length(toks)
    i = 1

    noscreen = false
    nologfile = false
    collisions = BOLSIGReadCollisions[]
    conditions = nothing
    runs = BOLSIGRunSpec[]
    save = nothing

    while i <= n
        kw = toks[i]
        if kw == "/NOSCREEN"
            noscreen = true
            i += 1
        elseif kw == "/NOLOGFILE"
            nologfile = true
            i += 1
        elseif kw == "CLEARCOLLISIONS"
            empty!(collisions)
            i += 1
        elseif kw == "CLEARRUNS"
            empty!(runs)
            i += 1
        elseif kw == "READCOLLISIONS"
            file = toks[i + 1]
            species = String.(split(toks[i + 2]))
            extrapolate = parse(Int, toks[i + 3]) == 1
            push!(collisions, BOLSIGReadCollisions(file, species, extrapolate))
            i += 4
        elseif kw == "CONDITIONS"
            n_species = sum(c -> length(c.species), collisions; init=0)
            conditions, consumed = _parse_bolsig_conditions(toks, i + 1, n_species)
            i += 1 + consumed
        elseif kw == "RUN"
            j = i + 1
            rows = Vector{Float64}[]
            while j <= n && !(toks[j] in _BOLSIG_INPUT_KEYWORDS)
                push!(rows, parse.(Float64, split(toks[j])))
                j += 1
            end
            push!(runs, isempty(rows) ? BOLSIGFixedRun() : BOLSIGExplicitRun(permutedims(reduce(hcat, rows))))
            i = j
        elseif kw == "RUNSERIES"
            variable = BOLSIGVariable(parse(Int, toks[i + 1]))
            mm = parse.(Float64, split(toks[i + 2]))
            count = parse(Int, toks[i + 3])
            type = BOLSIGSeriesType(parse(Int, toks[i + 4]))
            seg = BOLSIGSeriesSegment(variable, mm[1], mm[2], count, type)
            if !isempty(runs) && runs[end] isa BOLSIGSeriesRun
                push!(runs[end].segments, seg)
            else
                push!(runs, BOLSIGSeriesRun([seg]))
            end
            i += 5
        elseif kw == "RUN2D"
            variable1 = BOLSIGVariable(parse(Int, toks[i + 1]))
            r = parse.(Float64, split(toks[i + 2]))
            nn = parse.(Int, split(toks[i + 3]))
            tt = parse.(Int, split(toks[i + 4]))
            output_file = toks[i + 5]
            push!(runs, BOLSIGRun2D(
                variable1, r[1], r[2], nn[1], BOLSIGSeriesType(tt[1]),
                r[3], r[4], nn[2], BOLSIGSeriesType(tt[2]), output_file,
            ))
            i += 6
        elseif kw == "SAVERESULTS"
            save, consumed = _parse_bolsig_saveresults(toks, i + 1)
            i += 1 + consumed
        elseif kw == "END"
            break
        else
            error("read_bolsig_input: unrecognized token '$(kw)' at position $i.")
        end
    end

    isnothing(conditions) && error("read_bolsig_input: no CONDITIONS block found in $path.")

    return BOLSIGInput(; noscreen, nologfile, collisions, conditions, runs, save)
end

# ── Writing ──────────────────────────────────────────────────────────────────

# The ordered list of BOLSIGConditions positions currently marked VAR — the
# scalar fields (in their fixed textual order) followed by any VAR gas
# fractions. Shared between validation and (if ever needed) introspection.
function _bolsig_var_positions(c::BOLSIGConditions)
    positions = Symbol[]
    for f in _BOLSIG_CONDITIONS_SCALAR_FIELDS
        getfield(c, f) === nothing && push!(positions, f)
    end
    for (k, v) in enumerate(c.gas_fractions)
        v === nothing && push!(positions, Symbol("gas_fractions[$k]"))
    end
    return positions
end

# Only BOLSIGExplicitRun/BOLSIGRun2D rely on VAR-marker *position* to know
# which physical quantity a column/second-axis refers to (BOLSIGSeriesRun's
# variable is identified by an explicit numeric code instead, and never
# needs a matching VAR marker — confirmed against the reference examples:
# example1_input.dat's RUNSERIES scans E/N while its CONDITIONS entry for
# "Electric field / N (Td)" is a concrete placeholder, not VAR).
function _validate_bolsig_input(input::BOLSIGInput)
    var_positions = _bolsig_var_positions(input.conditions)
    for run in input.runs
        if run isa BOLSIGExplicitRun
            n_var = length(var_positions)
            size(run.values, 2) == n_var || error(
                "write_bolsig_input: a BOLSIGExplicitRun has $(size(run.values, 2)) " *
                "column(s) but `conditions` marks $n_var field(s) as VAR " *
                "($(var_positions)) — these must match so BOLSIG+ knows which " *
                "column corresponds to which field."
            )
        elseif run isa BOLSIGRun2D
            length(var_positions) == 1 || error(
                "write_bolsig_input: a BOLSIGRun2D's second variable is identified " *
                "by whichever `conditions` field is marked VAR, but " *
                "$(length(var_positions)) are marked ($(var_positions)) — exactly " *
                "one is required."
            )
        end
    end
    isnothing(input.save) && any(r -> !(r isa BOLSIGRun2D), input.runs) && error(
        "write_bolsig_input: `save` is `nothing` but `runs` contains a non-RUN2D " *
        "entry — only BOLSIGRun2D writes its own output file; everything else " *
        "needs a BOLSIGSaveResults."
    )
    return nothing
end

_bolsig_value_str(v::Nothing) = "VAR"
_bolsig_value_str(v::Bool) = v ? "1" : "0"
_bolsig_value_str(v::Real) = string(v)
_bolsig_value_str(v::Enum) = string(Int(v))

_write_bolsig_field(io, value, comment) = println(io, _bolsig_value_str(value), "   / ", comment)

function _write_bolsig_conditions(io, c::BOLSIGConditions)
    comments = Dict(
        :reduced_field => "Electric field / N (Td)",
        :angular_frequency => "Angular field frequency / N (m3/s)",
        :eb_angle_cosine => "Cosine of E-B field angle",
        :gas_temperature => "Gas temperature (K)",
        :excitation_temperature => "Excitation temperature (K)",
        :transition_energy => "Transition energy (eV)",
        :ionization_degree => "Ionization degree",
        :plasma_density => "Plasma density (1/m3)",
        :ion_charge_parameter => "Ion charge parameter",
        :ion_neutral_mass_ratio => "Ion/neutral mass ratio",
        :coulomb_collision_model => "e-e momentum effects & modified Coulomb logarithm: 0=No&No; 1=Yes&No; 2=No&Yes; 3=Yes&Yes*",
        :energy_sharing => "Energy sharing: 1=Equal*; 2=One takes all",
        :growth_model => "Growth: 1=Temporal*; 2=Spatial; 3=Not included; 4=Grad-n expansion",
        :maxwellian_mean_energy => "Maxwellian mean energy (eV)",
        :n_grid_points => "# of grid points",
        :manual_grid => "Manual grid: 0=No; 1=Linear; 2=Parabolic",
        :manual_max_energy => "Manual maximum energy (eV)",
        :precision => "Precision",
        :convergence => "Convergence",
        :max_iterations => "Maximum # of iterations",
    )
    for f in _BOLSIG_CONDITIONS_SCALAR_FIELDS
        _write_bolsig_field(io, getfield(c, f), comments[f])
    end
    println(io, join(_bolsig_value_str.(c.gas_fractions), "  "), "   / Gas composition fractions")
    _write_bolsig_field(io, c.normalize_fractions, "Normalize composition to unity: 0=No; 1=Yes")
end

_write_bolsig_run(io, run::BOLSIGFixedRun) = (println(io, "RUN"); println(io))

function _write_bolsig_run(io, run::BOLSIGExplicitRun)
    println(io, "RUN")
    for row in eachrow(run.values)
        println(io, join(row, "  "))
    end
    println(io)
end

function _write_bolsig_run(io, run::BOLSIGSeriesRun)
    for seg in run.segments
        println(io, "RUNSERIES")
        _write_bolsig_field(io, seg.variable, "Variable: 1=E/N; 2=Mean energy; 3=Maxwellian energy")
        println(io, seg.min, "  ", seg.max, "   / Min Max")
        _write_bolsig_field(io, seg.count, "Number")
        _write_bolsig_field(io, seg.type, "Type: 1=Linear; 2=Quadratic; 3=Exponential")
        println(io)
    end
end

function _write_bolsig_run(io, run::BOLSIGRun2D)
    println(io, "RUN2D")
    _write_bolsig_field(io, run.variable1, "First variable: 1=E/N; 2=Mean energy; 3=Maxwellian energy")
    println(io, run.min1, "  ", run.max1, "  ", run.min2, "  ", run.max2, "   / Min1 Max1 Min2 Max2")
    println(io, run.num1, "  ", run.num2, "   / Num1 Num2")
    println(io, _bolsig_value_str(run.type1), "  ", _bolsig_value_str(run.type2),
            "   / Type1 Type2: 1=Linear; 2=Quadratic; 3=Exponential")
    println(io, run.output_file, "   / Output file")
    println(io)
end

function _write_bolsig_save(io, save::BOLSIGSaveResults)
    println(io, "SAVERESULTS")
    println(io, save.file, "   / File")
    _write_bolsig_field(io, save.format, "Format: 1=Run by run; 2=Combined; 3=E/N; 4=Energy; 5=SIGLO; 6=PLASIMO")
    _write_bolsig_field(io, save.conditions, "Conditions: 0=No; 1=Yes")
    _write_bolsig_field(io, save.transport_coefficients, "Transport coefficients: 0=No; 1=Yes")
    _write_bolsig_field(io, save.rate_coefficients, "Rate coefficients: 0=No; 1=Yes")
    _write_bolsig_field(io, save.reverse_rate_coefficients, "Reverse rate coefficients: 0=No; 1=Yes")
    _write_bolsig_field(io, save.energy_loss_coefficients, "Energy loss coefficients: 0=No; 1=Yes")
    _write_bolsig_field(io, save.distribution_function, "Distribution function: 0=No; 1=Yes")
    _write_bolsig_field(io, save.skip_failed_runs, "Skip failed runs: 0=No; 1=Yes")
    _write_bolsig_field(io, save.include_cross_sections, "Include cross sections: 0=No; 1=Yes")
    println(io)
end

"""
    write_bolsig_input(path, input::BOLSIGInput)

Write `input` as a valid BOLSIG+ input script. Validates that `conditions`'
`VAR` markers are consistent with what `runs` actually declares as swept
(a [`BOLSIGExplicitRun`](@ref)'s column count, a [`BOLSIGRun2D`](@ref)'s
second variable) before writing — [`BOLSIGSeriesRun`](@ref) identifies its
variable by code and doesn't need or use a `VAR` marker.
"""
function write_bolsig_input(path::AbstractString, input::BOLSIGInput)
    _validate_bolsig_input(input)

    io = IOBuffer()
    input.noscreen && (println(io, "/NOSCREEN"); println(io))
    input.nologfile && (println(io, "/NOLOGFILE"); println(io))

    for c in input.collisions
        println(io, "READCOLLISIONS")
        println(io, c.file, "   / File")
        println(io, join(c.species, "  "), "   / Species")
        _write_bolsig_field(io, c.extrapolate, "Extrapolate: 0= No 1= Yes")
        println(io)
    end

    println(io, "CONDITIONS")
    _write_bolsig_conditions(io, input.conditions)
    println(io)

    for run in input.runs
        _write_bolsig_run(io, run)
    end

    isnothing(input.save) || _write_bolsig_save(io, input.save)

    println(io, "END")

    write(path, take!(io))
    return path
end
