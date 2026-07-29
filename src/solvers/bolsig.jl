struct BOLSIG <: Solver end

# ── Shared low-level helpers ─────────────────────────────────────────────────

# BOLSIG+ occasionally drops the "E" in an exponent when a fixed-width field
# would otherwise overflow with a 3-digit exponent, e.g. "-0.852184-208"
# instead of "-0.852184E-208". Handle both forms.
function _bolsig_parse_float(tok::AbstractString)
    m = match(r"^([+-]?[0-9.]+)([+-]\d{2,})$", tok)
    return parse(Float64, isnothing(m) ? tok : "$(m[1])E$(m[2])")
end

# Whitespace-split `line` and parse every token as a float; `nothing` if the
# line isn't purely numeric (used to distinguish data rows from block titles).
function _bolsig_numeric_tokens(line::AbstractString)
    toks = split(strip(line))
    isempty(toks) && return nothing
    try
        return _bolsig_parse_float.(toks)
    catch
        return nothing
    end
end

# Splits a "name    value" label line into (name, value::Union{Float64,Nothing}).
# Non-greedy so descriptive names containing internal runs of 2+ spaces (e.g.
# "Ar    Excitation    11.50 eV", which has no trailing bare value) are kept
# whole rather than truncated at the first internal gap.
const _BOLSIG_TRAILING_NUMBER =
    r"^(.*?)(?:\s{2,}([+-]?[0-9][0-9.]*(?:[Ee][+-]?\d+|[+-]\d{2,})?))?$"

function _bolsig_split_label(s::AbstractString)
    m = match(_BOLSIG_TRAILING_NUMBER, rstrip(s))
    name = strip(m[1])
    value = isnothing(m[2]) ? nothing : _bolsig_parse_float(m[2])
    return name, value
end

# Strips a leading "C<digits>" reaction-index prefix, if present, e.g.
# "C1    Ar    Effective (momentum)" -> "Ar    Effective (momentum)".
_bolsig_strip_reaction_index(s::AbstractString) = strip(replace(s, r"^C\d+\s+" => ""))

# BOLSIG+'s own descriptions aren't always unique — e.g. two distinct
# collisions can both display as "N2    Excitation    0.29 eV" when their
# thresholds round to the same two decimals. Disambiguate rather than
# silently dropping one via a dict-key/DataFrame-column collision. `seen`
# must be scoped to whatever set of names needs to be jointly unique (one
# condensed-format section, one full verbose-format file, one single-run
# block) and reused across every name assigned within that scope.
function _bolsig_dedupe_name!(seen::Dict{String,Int}, name::AbstractString)
    count = get(seen, name, 0) + 1
    seen[name] = count
    return count == 1 ? name : "$name (dup$(count))"
end

"""
    _bolsig_normalize_process(desc) -> String

Normalize a raw BOLSIG+ collision description into a compact reaction string,
e.g. `"Ar    Effective (momentum)"` -> `"Effective(Ar)"`,
`"Ar    Excitation    11.50 eV"` -> `"Excitation(Ar,11.50eV)"`. Falls back to
the original string, unmodified, if the description doesn't match the
expected `<species> <type> [threshold eV]` shape rather than guessing.
"""
function _bolsig_normalize_process(desc::AbstractString)
    # A "(dupN)" suffix (from _bolsig_dedupe_name!, when BOLSIG+ itself gives
    # two distinct collisions the identical description) must survive
    # normalization too, or both would collapse back onto the same name.
    m_dup = match(r"^(.*?)( \(dup\d+\))$", desc)
    base, suffix = isnothing(m_dup) ? (desc, "") : (m_dup[1], m_dup[2])

    m = match(r"^(\S+)\s+(\w+)(?:\s*\(momentum\))?(?:\s+([\d.]+)\s*eV)?", base)
    isnothing(m) && return desc
    species, kind, threshold = m[1], m[2], m[3]
    normalized = isnothing(threshold) ? "$(kind)($(species))" : "$(kind)($(species),$(threshold)eV)"
    return normalized * suffix
end

# ── Format detection ─────────────────────────────────────────────────────────
# BOLSIG+ can emit (at least) three structurally different result layouts:
#   :condensed   — "R#"/"A#"/"C#" indexed tables ( Conditions/Transport
#                  coefficients/Rate coefficients sections)
#   :single_run  — one "R#" header per run, followed by plain "name  value"
#                  lines (no table wrapping, since there's only one value)
#   :verbose     — one block per quantity, title line + a value table
#                  wrapped across multiple lines, keyed off a leading
#                  "Electric field / N (Td)" block
function _detect_bolsig_format(lines::AbstractVector{<:AbstractString})
    any(l -> occursin(r"^\s*R#\s", l), lines) && return :condensed
    any(l -> occursin(r"^R\d+\s*$", strip(l)), lines) && return :single_run
    return :verbose
end

# ── Condensed format (R#/A#/C# indexed tables) ───────────────────────────────

# Parses one condensed-format section (" Conditions"/" Transport
# coefficients"/"Rate coefficients (m3/s)") starting at line index `i`
# (the section title itself). Returns (column_names, data, reaction_names)
# where `data` is an (n_runs × n_columns) Float64 matrix aligned with
# `column_names`, and `reaction_names` lists the subset of `column_names`
# that came from a "C#" label (i.e. are reaction rate coefficients).
function _bolsig_condensed_section(lines, i)
    n = lastindex(lines)
    idx = i + 1
    label_map = Dict{String,String}()
    label_order = String[]
    while idx <= n
        m = match(r"^\s*([AC]\d+)\s+(.*)$", lines[idx])
        isnothing(m) && break
        name, _ = _bolsig_split_label(m[2])
        label_map[m[1]] = name
        push!(label_order, m[1])
        idx += 1
    end

    raw_tokens = split(strip(lines[idx]))
    keys = String[]
    for t in raw_tokens
        if startswith(t, "(") && !isempty(keys)
            keys[end] = keys[end] * " " * t
        else
            push!(keys, t)
        end
    end
    idx += 1

    rows = Vector{Float64}[]
    while idx <= n && !isempty(strip(lines[idx]))
        push!(rows, _bolsig_parse_float.(split(strip(lines[idx]))))
        idx += 1
    end

    seen = Dict{String,Int}()
    deduped = Dict(id => _bolsig_dedupe_name!(seen, label_map[id]) for id in label_order)

    column_names = [get(deduped, k, k) for k in keys]
    data = permutedims(reduce(hcat, rows))
    reaction_names = [deduped[k] for k in label_order if startswith(k, "C")]
    return column_names, data, reaction_names
end

function _load_bolsig_condensed(lines)
    i_transport = findfirst(l -> strip(l) == "Transport coefficients", lines)
    i_rate = findfirst(l -> startswith(strip(l), "Rate coefficients"), lines)
    (isnothing(i_transport) || isnothing(i_rate)) &&
        error("BOLSIG+ condensed format: couldn't find 'Transport coefficients'/'Rate coefficients' sections.")

    names_t, data_t, _ = _bolsig_condensed_section(lines, i_transport)
    names_r, data_r, reaction_names = _bolsig_condensed_section(lines, i_rate)

    df_t = DataFrame(data_t, names_t)
    df_r = DataFrame(data_r, names_r)
    "E/N (Td)" in names(df_r) && select!(df_r, Not("E/N (Td)"))

    return innerjoin(df_t, df_r, on="R#"), reaction_names
end

# ── Verbose format (block-per-quantity, wrapped multi-line tables) ──────────

function _load_bolsig_verbose(lines)
    n = lastindex(lines)
    i = findfirst(l -> strip(l) == "Electric field / N (Td)", lines)
    isnothing(i) &&
        error("BOLSIG+ verbose format: couldn't find the 'Electric field / N (Td)' results header.")

    columns = Dict{String,Vector{Float64}}()
    seen = Dict{String,Int}()
    reaction_names = String[]
    idx = i
    while idx <= n
        while idx <= n && isempty(strip(lines[idx]))
            idx += 1
        end
        idx > n && break

        title = strip(lines[idx])
        idx += 1
        is_reaction = occursin(r"^C\d+\s", title)
        is_reaction && (title = _bolsig_strip_reaction_index(title))

        while idx <= n && isempty(strip(lines[idx]))
            idx += 1
        end
        idx > n && break

        toks = _bolsig_numeric_tokens(lines[idx])
        if isnothing(toks)
            # Secondary units/descriptor line (e.g. "Rate coefficient (m3/s)")
            # — folded into the title for non-reaction blocks; dropped for
            # reaction blocks so their name matches the condensed format's.
            is_reaction || (title = title * "  " * strip(lines[idx]))
            idx += 1
            while idx <= n && isempty(strip(lines[idx]))
                idx += 1
            end
            idx > n && break
        end

        values = Float64[]
        while idx <= n
            toks = _bolsig_numeric_tokens(lines[idx])
            isnothing(toks) && break
            append!(values, toks)
            idx += 1
        end
        isempty(values) && continue

        key = _bolsig_dedupe_name!(seen, title)
        columns[key] = values
        is_reaction && push!(reaction_names, key)
    end

    ncond = length(columns["Electric field / N (Td)"])
    for (k, v) in columns
        length(v) == ncond ||
            error("BOLSIG+ verbose format: column '$k' has $(length(v)) values, expected $ncond.")
    end

    return DataFrame(columns), reaction_names
end

# ── Single-run report format ("R#" header + plain name/value lines) ────────

function _load_bolsig_single_run(lines)
    n = lastindex(lines)
    run_starts = findall(l -> occursin(r"^R\d+\s*$", strip(l)), lines)
    isempty(run_starts) &&
        error("BOLSIG+ single-run format: couldn't find any 'R#' run header.")

    rows = Dict{String,Float64}[]
    reaction_names = Set{String}()
    for (k, start) in enumerate(run_starts)
        stop = k < length(run_starts) ? run_starts[k + 1] - 1 : n
        row = Dict{String,Float64}()
        seen = Dict{String,Int}()   # scoped per run block — see _bolsig_dedupe_name!
        for idx in (start + 1):stop
            s = strip(lines[idx])
            (isempty(s) || occursin(r"^-+$", s) || startswith(s, "Rate coefficients")) && continue
            name, value = _bolsig_split_label(s)
            isnothing(value) && continue
            is_reaction = occursin(r"^C\d+\s", name)
            is_reaction && (name = _bolsig_strip_reaction_index(name))
            name = _bolsig_dedupe_name!(seen, name)
            row[name] = value
            is_reaction && push!(reaction_names, name)
        end
        push!(rows, row)
    end

    allnames = sort!(collect(reduce(union, keys.(rows))))
    df = DataFrame([nm => Vector{Union{Missing,Float64}}(missing, length(rows)) for nm in allnames])
    for (i, row) in enumerate(rows), (k, v) in row
        df[i, k] = v
    end

    return df, collect(reaction_names)
end

# ── Top-level dispatch ───────────────────────────────────────────────────────

function _load_bolsig(lines)
    fmt = _detect_bolsig_format(lines)
    fmt === :condensed && return _load_bolsig_condensed(lines)
    fmt === :single_run && return _load_bolsig_single_run(lines)
    return _load_bolsig_verbose(lines)
end

"""
    load_raw_dataframe(::BOLSIG, source; kwargs...)

Read a single BOLSIG+ output file (`source` is a file path, unlike `LoKI`/
`MultiBolt` which take a directory — BOLSIG+ writes one output file per run).
Auto-detects which of BOLSIG+'s three result layouts the file uses (condensed
`R#`/`A#`/`C#` tables, a verbose block-per-quantity wrapped-table format, or a
single-run `R#`-header name/value report) and returns one row per scanned
condition (or a single row, for the single-run format).
"""
function load_raw_dataframe(::BOLSIG, source; kwargs...)
    df, _ = _load_bolsig(readlines(source))
    return df
end

function default_swarm_names(::BOLSIG)
    return [
        "E/N (Td)" => :reduced_field,
        "Electric field / N (Td)" => :reduced_field,
        "Mean energy (eV)" => :mean_energy,
        "Mobility *N (1/m/V/s)" => :reduced_mobility,
        "Diffusion coefficient *N (1/m/s)" => :reduced_diffusion_coef,
        "Energy mobility *N (1/m/V/s)" => :reduced_energy_mobility,
        "Energy diffusion coef. *N (1/m/s)" => :reduced_energy_diffusion_coef,
        "Townsend ioniz. coef. alpha/N (m2)" => :reduced_townsend_alpha_coef,
    ]
end

"""
    parse_reaction_names(::BOLSIG, source)

Map each rate-coefficient column produced by [`load_raw_dataframe`](@ref) to a
normalized reaction string via [`_bolsig_normalize_process`](@ref).
"""
function parse_reaction_names(::BOLSIG, source)
    _, reaction_names = _load_bolsig(readlines(source))
    return [name => _bolsig_normalize_process(name) for name in reaction_names]
end
