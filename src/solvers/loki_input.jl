# Reading LoKI-B `.in` setup files (as opposed to loki.jl, which reads
# LoKI *output*). The format is indentation-based with `%` comments,
# `key: value` scalars (MATLAB arithmetic, `linspace`/`logspace` allowed),
# `- item` lists and `- key = value` pair lists. `read_loki_input` parses
# the whole file and maps the parts that describe an electron-kinetics
# (Boltzmann) run into a `LoKIInput`; everything else stays available
# under `raw`. Species and states are kept as plain strings — resolving
# them is a consumer's job (BoltzmannSolvers has no species types).

# ── Value parsing ────────────────────────────────────────────────────────────

# Safe evaluator for MATLAB-style values: numbers, + - * / ^, and the
# `linspace(a,b,n)`/`logspace(a,b,n)` range builders LoKI setup files use
# for condition scans. Anything else falls back to the raw string.
_loki_eval(x::Number) = Float64(x)
function _loki_eval(e::Expr)
    if e.head === :vect
        return Float64[_loki_eval(a) for a in e.args]
    end
    e.head === :call || error("unsupported expression")
    op = e.args[1]
    if op === :linspace || op === :logspace
        a, b = _loki_eval(e.args[2]), _loki_eval(e.args[3])
        n = Int(_loki_eval(e.args[4]))
        r = range(a, b; length=n)
        return op === :linspace ? collect(r) : exp10.(r)
    end
    args = [_loki_eval(a) for a in e.args[2:end]]
    op === :+ && return sum(args)
    op === :- && return length(args) == 1 ? -args[1] : args[1] - args[2]
    op === :* && return prod(args)
    op === :/ && return args[1] / args[2]
    op === :^ && return args[1]^args[2]
    error("unsupported operator")
end
_loki_eval(x) = error("unsupported token")

function _loki_value(str::AbstractString)
    s = strip(str)
    s == "true" && return true
    s == "false" && return false
    v = try
        _loki_eval(Meta.parse(s))
    catch
        nothing
    end
    return v === nothing ? String(s) : v
end

# ── Structure parsing ────────────────────────────────────────────────────────

# Split a `- key = value` list entry at the first `=` outside parentheses
# (state labels contain `=` themselves, e.g. `N2(X,v=*) = 1.0`); `nothing`
# for plain `- item` entries.
function _split_pair_entry(s::AbstractString)
    depth = 0
    for (i, c) in pairs(s)
        c == '(' && (depth += 1)
        c == ')' && (depth -= 1)
        if c == '=' && depth == 0
            return (String(strip(s[1:prevind(s, i)])),
                    String(strip(s[nextind(s, i):end])))
        end
    end
    return nothing
end

# Parse the indented structure into nested Dicts / Vectors (same scheme as
# LoKI-GM setup files: a block is either a map of `key: ...` entries or a
# list of `- ...` entries).
function _parse_loki_setup(lines::Vector{String})
    items = Tuple{Int,String}[]
    for line in lines
        code = rstrip(first(split(line, '%'; limit=2)))
        isempty(strip(code)) && continue
        indent = length(code) - length(lstrip(code))
        push!(items, (indent, strip(code)))
    end
    node, _ = _parse_loki_block(items, 1, -1)
    return node
end

function _parse_loki_block(items, i, parent_indent)
    dict = Dict{String,Any}()
    list = Any[]
    while i <= length(items)
        indent, content = items[i]
        indent <= parent_indent && break
        if startswith(content, "- ")
            entry = strip(content[3:end])
            kv = _split_pair_entry(entry)
            push!(list, kv === nothing ? _loki_value(entry) :
                        (kv[1] => _loki_value(kv[2])))
            i += 1
        else
            m = match(r"^([^:]+):\s*(.*)$", content)
            m === nothing && error("read_loki_input: cannot parse line \"$content\"")
            key = String(strip(m[1]))
            val = strip(m[2])
            if isempty(val)
                child, i = _parse_loki_block(items, i + 1, indent)
                dict[key] = child
            else
                dict[key] = _loki_value(val)
                i += 1
            end
        end
    end
    isempty(dict) || isempty(list) ||
        error("read_loki_input: mixed list/map block in setup file")
    return (isempty(dict) ? list : dict), i
end

function _loki_get(d, path...)
    x = d
    for k in path
        x isa Dict || return nothing
        x = get(x, k, nothing)
    end
    return x
end

# Scalars and vectors are interchangeable in LoKI conditions
# (`reducedField: 10` vs `reducedField: logspace(-3,3,100)`).
_loki_vector(::Nothing, default) = default
_loki_vector(x::Number, _) = [Float64(x)]
_loki_vector(x::AbstractVector, _) = Float64.(x)

# `LXCatFiles: file.txt` (scalar) and the list form are both allowed.
_loki_strings(::Nothing) = String[]
_loki_strings(x::AbstractString) = [String(x)]
_loki_strings(x::AbstractVector) = String.(x)

_loki_scalar(x, default) = x === nothing ? default : x

"""
    LoKISmartGrid(; min_eedf_decay=20.0, max_eedf_decay=25.0, update_factor=0.05)

LoKI-B's adaptive energy-grid settings (`numerics.energyGrid.smartGrid`):
the grid boundary is moved until the EEDF falls between `min_eedf_decay`
and `max_eedf_decay` decades below its maximum.
"""
Base.@kwdef struct LoKISmartGrid
    min_eedf_decay::Float64 = 20.0
    max_eedf_decay::Float64 = 25.0
    update_factor::Float64 = 0.05
end

"""
    LoKIInput(; kwargs...)

Configuration of a LoKI-B electron-kinetics run, as described by the
`workingConditions` and `electronKinetics` sections of a `.in` setup file
(see [`read_loki_input`](@ref)). Species, states, and model choices are
plain strings in LoKI's own notation; `raw` holds the full parsed setup
tree (chemistry, output, GUI, ... configuration).

Conditions given as `linspace`/`logspace`/vector scans are kept as
vectors (`reduced_field`, `electron_temperature`).
"""
Base.@kwdef struct LoKIInput
    # workingConditions
    reduced_field::Vector{Float64} = Float64[]      # Td
    excitation_frequency::Float64 = 0.0             # Hz
    gas_pressure::Float64 = NaN                     # Pa
    gas_temperature::Float64 = 300.0                # K
    electron_temperature::Vector{Float64} = Float64[]  # eV (prescribedEedf scans)
    electron_density::Float64 = NaN                 # m^-3

    # electronKinetics
    eedf_type::String = "boltzmann"                 # boltzmann | prescribedEedf
    shape_parameter::Union{Float64,Nothing} = nothing  # 1=Maxwellian .. 2=Druyvesteyn
    ionization_operator::String = "conservative"    # conservative | oneTakesAll | equalSharing | usingSDCS
    growth_model::String = "temporal"               # temporal | spatial
    include_ee::Bool = false
    lxcat_files::Vector{String} = String[]
    lxcat_extra_files::Vector{String} = String[]
    car_gases::Vector{String} = String[]
    fractions::Vector{Pair{String,Float64}} = Pair{String,Float64}[]
    # population entries verbatim: `"state" => number`, `"state" => "function"`
    # (e.g. "boltzmannPopulation@gasTemperature"), or a bare file-reference
    # string — resolving non-numeric entries is left to the consumer.
    populations::Vector{Any} = Any[]

    # numerics
    max_energy::Float64 = 1.0                       # eV
    cell_number::Int = 1000
    smart_grid::Union{LoKISmartGrid,Nothing} = nothing
    mixing_parameter::Float64 = 0.7
    max_eedf_rel_error::Float64 = 1e-9

    raw::Dict{String,Any} = Dict{String,Any}()
end

"""
    read_loki_input(path) -> LoKIInput

Parse a LoKI-B `.in` setup file. The `workingConditions`,
`electronKinetics` (cross-section files, gas fractions, state populations,
growth/ionization/e-e model choices) and `numerics` sections map onto
[`LoKIInput`](@ref) fields; the complete parsed tree — including sections
this struct does not model (`chemistry`, `output`, state energies and
statistical weights, ...) — is kept in `raw`.
"""
function read_loki_input(path::AbstractString)
    tree = _parse_loki_setup(readlines(path))
    tree isa Dict || error("read_loki_input: $path is not a map at top level.")

    wc(k...) = _loki_get(tree, "workingConditions", k...)
    ek(k...) = _loki_get(tree, "electronKinetics", k...)
    ek() === nothing && error("read_loki_input: no electronKinetics section in $path.")

    smart = ek("numerics", "energyGrid", "smartGrid")
    smart_grid = smart === nothing ? nothing : LoKISmartGrid(;
        min_eedf_decay=Float64(_loki_scalar(_loki_get(smart, "minEedfDecay"), 20.0)),
        max_eedf_decay=Float64(_loki_scalar(_loki_get(smart, "maxEedfDecay"), 25.0)),
        update_factor=Float64(_loki_scalar(_loki_get(smart, "updateFactor"), 0.05)),
    )

    fractions = Pair{String,Float64}[]
    for p in something(ek("gasProperties", "fraction"), Any[])
        p isa Pair && last(p) isa Number || error(
            "read_loki_input: cannot parse gas fraction entry \"$p\".")
        push!(fractions, String(first(p)) => Float64(last(p)))
    end

    return LoKIInput(;
        reduced_field=_loki_vector(wc("reducedField"), Float64[]),
        excitation_frequency=Float64(_loki_scalar(wc("excitationFrequency"), 0.0)),
        gas_pressure=Float64(_loki_scalar(wc("gasPressure"), NaN)),
        gas_temperature=Float64(_loki_scalar(wc("gasTemperature"), 300.0)),
        electron_temperature=_loki_vector(wc("electronTemperature"), Float64[]),
        electron_density=Float64(_loki_scalar(wc("electronDensity"), NaN)),
        eedf_type=String(_loki_scalar(ek("eedfType"), "boltzmann")),
        shape_parameter=(sp = ek("shapeParameter"); sp === nothing ? nothing : Float64(sp)),
        ionization_operator=String(_loki_scalar(ek("ionizationOperatorType"), "conservative")),
        growth_model=String(_loki_scalar(ek("growthModelType"), "temporal")),
        include_ee=_loki_scalar(ek("includeEECollisions"), false),
        lxcat_files=_loki_strings(ek("LXCatFiles")),
        lxcat_extra_files=_loki_strings(ek("LXCatFilesExtra")),
        car_gases=_loki_strings(ek("CARgases")),
        fractions,
        populations=collect(Any, something(ek("stateProperties", "population"), Any[])),
        max_energy=Float64(_loki_scalar(ek("numerics", "energyGrid", "maxEnergy"), 1.0)),
        cell_number=Int(_loki_scalar(ek("numerics", "energyGrid", "cellNumber"), 1000)),
        smart_grid,
        mixing_parameter=Float64(_loki_scalar(
            ek("numerics", "nonLinearRoutines", "mixingParameter"), 0.7)),
        max_eedf_rel_error=Float64(_loki_scalar(
            ek("numerics", "nonLinearRoutines", "maxEedfRelError"), 1e-9)),
        raw=tree,
    )
end
