# Configuration for a MultiBolt run and its translation into command-line
# arguments (as opposed to bolsig_input.jl, which reads/writes an actual
# input *file* — MultiBolt has no such format; it's entirely CLI-argument
# driven, confirmed by running `multibolt_linux --help` and the tool's own
# `bin/script_examples/sh_script/*.sh`). `_multibolt_args` plays the role
# `write_bolsig_input` plays for BOLSIG+: turning a config struct into
# exactly what the external binary consumes.

@enum MultiBoltModelType HDModel HDGEModel HDGE01Model SSTModel
@enum MultiBoltSweepVariable ENTdSweep TKSweep pTorrSweep BinFracSweep NuSweep NTermsSweep
@enum MultiBoltScatteringModel IsotropicScattering IdealForwardScattering ScreenedCoulombScattering
@enum MultiBoltInterpMethod LinearInterp LogarithmicInterp

_multibolt_model_str(m::MultiBoltModelType) = m == HDModel ? "HD" :
    m == HDGEModel ? "HD+GE" : m == HDGE01Model ? "HD+GE_01" : "SST"

_multibolt_sweep_var_str(v::MultiBoltSweepVariable) = v == ENTdSweep ? "EN_Td" :
    v == TKSweep ? "T_K" : v == pTorrSweep ? "p_Torr" : v == BinFracSweep ? "bin_frac" :
    v == NuSweep ? "Nu" : "N_terms"

_multibolt_scattering_str(s::MultiBoltScatteringModel) = s == IsotropicScattering ? "Isotropic" :
    s == IdealForwardScattering ? "IdealForward" : "ScreenedCoulomb"

_multibolt_interp_str(m::MultiBoltInterpMethod) = m == LinearInterp ? "Linear" : "Logarithmic"

abstract type MultiBoltSweepStyle end

"""
    MultiBoltLinearSweep(start, stop, count)

`--sweep_style lin start stop count`.
"""
struct MultiBoltLinearSweep <: MultiBoltSweepStyle
    start::Float64
    stop::Float64
    count::Int
end

"""
    MultiBoltLogSweep(log10_start, log10_stop, count)

`--sweep_style log log10(start) log10(stop) count`.
"""
struct MultiBoltLogSweep <: MultiBoltSweepStyle
    log10_start::Float64
    log10_stop::Float64
    count::Int
end

"""
    MultiBoltRegularSweep(start, step, stop)

`--sweep_style reg start step stop`.
"""
struct MultiBoltRegularSweep <: MultiBoltSweepStyle
    start::Float64
    step::Float64
    stop::Float64
end

"""
    MultiBoltDefinedSweep(values)

`--sweep_style def v1 v2 ...` — an explicit list of values, MultiBolt's
manual/arbitrary-list scan (the analog of BOLSIG+'s `BOLSIGExplicitRun`).
"""
struct MultiBoltDefinedSweep <: MultiBoltSweepStyle
    values::Vector{Float64}
end

function _multibolt_sweep_style_args(s::MultiBoltLinearSweep)
    return ["--sweep_style", "lin", string(s.start), string(s.stop), string(s.count)]
end
function _multibolt_sweep_style_args(s::MultiBoltLogSweep)
    return ["--sweep_style", "log", string(s.log10_start), string(s.log10_stop), string(s.count)]
end
function _multibolt_sweep_style_args(s::MultiBoltRegularSweep)
    return ["--sweep_style", "reg", string(s.start), string(s.step), string(s.stop)]
end
function _multibolt_sweep_style_args(s::MultiBoltDefinedSweep)
    return ["--sweep_style", "def", string.(s.values)...]
end

"""
    MultiBoltSweep(variable, style)

`--sweep_option <variable>` + the args from `style`. `variable`'s current
value in [`MultiBoltInput`](@ref) (e.g. `EN_Td`) acts as a placeholder,
silently overridden by the sweep — no `VAR`-marker mechanism is needed here
the way BOLSIG+'s `CONDITIONS` block needs one (confirmed against MultiBolt's
own `example_sweep_calculations_EN_Td.sh`, which passes a concrete `--EN_Td
100` alongside `--sweep_option EN_Td`).
"""
struct MultiBoltSweep
    variable::MultiBoltSweepVariable
    style::MultiBoltSweepStyle
end

"""
    MultiBoltSpecies(name, fraction)

`--species name fraction`. `name` is matched against whatever species are
identifiable (as reactants/products) in the files listed in
`MultiBoltInput.cross_section_files` — there's no explicit file→species
pairing at the CLI level the way `BOLSIGReadCollisions` pairs one file with
its species.
"""
struct MultiBoltSpecies
    name::String
    fraction::Float64
end

"""
    MultiBoltXsecScale(process, scale)

`--scale_Xsec process scale` — a scalar weight applied to one named
collision process.
"""
struct MultiBoltXsecScale
    process::String
    scale::Float64
end

"""
    MultiBoltInput(; kwargs...)

Configuration for one MultiBolt run. Unlike [`BOLSIGInput`](@ref), this
isn't read from/written to a file — MultiBolt takes its configuration
entirely as command-line arguments — so `MultiBoltInput` is only ever
translated into an argument vector (via `_multibolt_args`, used by
[`run_multibolt`](@ref)), not serialized.

Required: `cross_section_files`, `species`, `export_name`. Everything else
has a default matching MultiBolt's own CLI defaults, or `nothing`/`false`
meaning "don't pass this flag, let MultiBolt use its own default".
"""
Base.@kwdef struct MultiBoltInput
    # Basic parameters — EN_Td/T_K/p_Torr/Nu/N_terms are placeholders,
    # silently overridden if `sweep` targets that same quantity.
    model::MultiBoltModelType = HDGEModel
    N_terms::Int = 6
    Nu::Int = 200
    p_Torr::Float64 = 760.0
    T_K::Float64 = 300.0
    EN_Td::Float64 = 100.0

    # Cross sections / composition
    cross_section_files::Vector{String}
    species::Vector{MultiBoltSpecies}
    keep_all_species::Bool = false
    xsec_scales::Vector{MultiBoltXsecScale} = MultiBoltXsecScale[]

    sweep::Union{MultiBoltSweep,Nothing} = nothing   # nothing = single run at the conditions above (see _multibolt_args: internally becomes a 1-point `def` sweep at EN_Td, since a truly bare invocation crashes the real binary)

    # Export. `export_location` is deliberately not a field here —
    # `run_multibolt` controls it (a fresh temp directory per run), the
    # same way BOLSIGSaveResults.file is a bare filename within a
    # run_bolsig-managed workdir rather than an absolute path.
    export_name::String
    limit_export::Bool = false
    export_xsecs::Bool = false

    # Scattering. Elastic's valid choices are really only
    # {Isotropic, ScreenedCoulomb} (MultiBolt's own docs call IdealForward
    # "non-physical" for elastic collisions specifically) — not enforced
    # here at the type level, just documented; MultiBolt will reject it
    # itself if misused.
    elastic_scattering::Union{MultiBoltScatteringModel,Nothing} = nothing
    excitation_scattering::Union{MultiBoltScatteringModel,Nothing} = nothing
    ionization_scattering::Union{MultiBoltScatteringModel,Nothing} = nothing
    superelastic_scattering::Union{MultiBoltScatteringModel,Nothing} = nothing
    sharing::Union{Float64,Nothing} = nothing

    # Advanced / convergence
    conv_err::Union{Float64,Nothing} = nothing
    iter_max::Union{Int,Nothing} = nothing
    iter_min::Union{Int,Nothing} = nothing
    initial_eV_max::Union{Float64,Nothing} = nothing
    use_ev_max_guess::Bool = false
    weight_f0::Union{Float64,Nothing} = nothing

    # Remap
    use_energy_remap::Bool = false
    remap_target_order_span::Union{Float64,Nothing} = nothing
    remap_grid_trial_max::Union{Int,Nothing} = nothing
    remap_allowance::Union{Float64,Nothing} = nothing
    use_nu_remap::Bool = false
    remap_Nu_max::Union{Int,Nothing} = nothing
    remap_Nu_increment::Union{Int,Nothing} = nothing

    # Misc
    silent::Bool = false
    num_threads::Union{Int,Nothing} = nothing
    shy::Bool = false
    dont_enforce_sum::Bool = false
    interp_method::Union{MultiBoltInterpMethod,Nothing} = nothing
end

function _validate_multibolt_input(c::MultiBoltInput)
    if !isnothing(c.sweep) && c.sweep.variable == BinFracSweep && length(c.species) != 2
        error(
            "MultiBoltInput: sweep variable is bin_frac, which requires exactly 2 species " *
            "(the first gets the swept fraction, the second gets 1-fraction, and any others " *
            "are silently removed by MultiBolt) — got $(length(c.species))."
        )
    end
    return nothing
end

"""
    _multibolt_args(config::MultiBoltInput, export_location::AbstractString) -> Vector{String}

Translate `config` into the argument vector `multibolt_linux`/
`multibolt_win64.exe` expects, with `--export_location` fixed to
`export_location` (managed by [`run_multibolt`](@ref)).
"""
function _multibolt_args(c::MultiBoltInput, export_location::AbstractString)
    _validate_multibolt_input(c)

    args = String[
        "--model", _multibolt_model_str(c.model),
        "--N_terms", string(c.N_terms),
        "--Nu", string(c.Nu),
        "--p_Torr", string(c.p_Torr),
        "--T_K", string(c.T_K),
        "--EN_Td", string(c.EN_Td),
    ]

    for f in c.cross_section_files
        append!(args, ["--LXCat_Xsec_fid", f])
    end
    for s in c.species
        append!(args, ["--species", s.name, string(s.fraction)])
    end
    c.keep_all_species && push!(args, "--KEEP_ALL_SPECIES")
    for xs in c.xsec_scales
        append!(args, ["--scale_Xsec", xs.process, string(xs.scale)])
    end

    # A *truly* bare invocation — no --sweep_option/--sweep_style at all —
    # reliably crashes multibolt_linux with std::bad_alloc (confirmed
    # directly: reproducible regardless of Nu/N_terms/convergence settings,
    # so it's a real bug in the binary's single-run code path, not a
    # settings issue). Route `sweep = nothing` through a trivial 1-point
    # `def` sweep at the current EN_Td instead — semantically identical to a
    # bare single run (EN_Td doesn't change) but uses the sweep code path,
    # which doesn't crash.
    sweep = something(c.sweep, MultiBoltSweep(ENTdSweep, MultiBoltDefinedSweep([c.EN_Td])))
    append!(args, ["--sweep_option", _multibolt_sweep_var_str(sweep.variable)])
    append!(args, _multibolt_sweep_style_args(sweep.style))

    append!(args, ["--export_location", export_location, "--export_name", c.export_name])
    c.limit_export && push!(args, "--LIMIT_EXPORT")
    c.export_xsecs && push!(args, "--EXPORT_XSECS")

    isnothing(c.elastic_scattering) || append!(args, ["--elastic_scattering", _multibolt_scattering_str(c.elastic_scattering)])
    isnothing(c.excitation_scattering) || append!(args, ["--excitation_scattering", _multibolt_scattering_str(c.excitation_scattering)])
    isnothing(c.ionization_scattering) || append!(args, ["--ionization_scattering", _multibolt_scattering_str(c.ionization_scattering)])
    isnothing(c.superelastic_scattering) || append!(args, ["--superelastic_scattering", _multibolt_scattering_str(c.superelastic_scattering)])
    isnothing(c.sharing) || append!(args, ["--sharing", string(c.sharing)])

    isnothing(c.conv_err) || append!(args, ["--conv_err", string(c.conv_err)])
    isnothing(c.iter_max) || append!(args, ["--iter_max", string(c.iter_max)])
    isnothing(c.iter_min) || append!(args, ["--iter_min", string(c.iter_min)])
    isnothing(c.initial_eV_max) || append!(args, ["--initial_eV_max", string(c.initial_eV_max)])
    c.use_ev_max_guess && push!(args, "--USE_EV_MAX_GUESS")
    isnothing(c.weight_f0) || append!(args, ["--weight_f0", string(c.weight_f0)])

    if c.use_energy_remap
        push!(args, "--USE_ENERGY_REMAP")
        isnothing(c.remap_target_order_span) || append!(args, ["--remap_target_order_span", string(c.remap_target_order_span)])
        isnothing(c.remap_grid_trial_max) || append!(args, ["--remap_grid_trial_max", string(c.remap_grid_trial_max)])
        isnothing(c.remap_allowance) || append!(args, ["--remap_allowance", string(c.remap_allowance)])
    end
    if c.use_nu_remap
        push!(args, "--USE_NU_REMAP")
        isnothing(c.remap_Nu_max) || append!(args, ["--remap_Nu_max", string(c.remap_Nu_max)])
        isnothing(c.remap_Nu_increment) || append!(args, ["--remap_Nu_increment", string(c.remap_Nu_increment)])
    end

    c.silent && push!(args, "--SILENT")
    isnothing(c.num_threads) || append!(args, ["--multibolt_num_threads", string(c.num_threads)])
    c.shy && push!(args, "--SHY")
    c.dont_enforce_sum && push!(args, "--DONT_ENFORCE_SUM")
    isnothing(c.interp_method) || append!(args, ["--interp_method", _multibolt_interp_str(c.interp_method)])

    return args
end
