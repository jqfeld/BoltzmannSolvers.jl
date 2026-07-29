using BoltzmannSolvers
using Test

const DATA = joinpath(@__DIR__, "data", "bolsig_input")

# Generic recursive structural-equality check: these config structs don't
# (and shouldn't) define a custom `Base.==` — the default falls back to
# identity for structs holding heap fields (Vector/String), so two
# separately-parsed-but-equivalent BOLSIGInputs would otherwise compare
# unequal. Used only here, for round-trip testing.
function _deep_eq(a, b)
    if a isa AbstractArray
        return b isa AbstractArray && length(a) == length(b) && all(_deep_eq(x, y) for (x, y) in zip(a, b))
    end
    (a isa Number || a isa AbstractString || a isa Symbol || a isa Nothing || a isa Enum) && return a == b
    typeof(a) == typeof(b) || return false
    return all(_deep_eq(getfield(a, f), getfield(b, f)) for f in fieldnames(typeof(a)))
end

@testset "read_bolsig_input: example1 (RUNSERIES, chained segments)" begin
    input = read_bolsig_input(joinpath(DATA, "example1_input.dat"))

    @test length(input.collisions) == 2
    @test input.collisions[1].species == ["Ar"]
    @test input.collisions[2].species == ["He"]
    @test all(c -> c.extrapolate, input.collisions)

    @test input.conditions.reduced_field == 10.0
    @test input.conditions.gas_temperature == 300.0
    @test input.conditions.gas_fractions == [0.1, 0.9]
    @test input.conditions.normalize_fractions == true

    @test length(input.runs) == 1
    run = input.runs[1]
    @test run isa BOLSIGSeriesRun
    @test length(run.segments) == 2
    @test run.segments[1] == BOLSIGSeriesSegment(ReducedFieldVar, 0.1, 100.0, 31, ExponentialSeries)
    @test run.segments[2] == BOLSIGSeriesSegment(ReducedFieldVar, 200.0, 1000.0, 9, LinearSeries)
    # total points (31+9=40) matches the already-verified example1.dat output size
    @test run.segments[1].count + run.segments[2].count == 40

    @test input.save isa BOLSIGSaveResults
    @test input.save.file == "example1.dat"
    @test input.save.format == 3   # "E/N" -> the verbose wrapped-table reader
end

@testset "read_bolsig_input: example2 (VAR + explicit RUN blocks)" begin
    input = read_bolsig_input(joinpath(DATA, "example2_input.dat"))

    @test input.conditions.reduced_field === nothing   # VAR
    @test input.conditions.angular_frequency === nothing   # VAR
    @test input.conditions.gas_fractions == [0.99, 0.01]

    @test length(input.runs) == 3
    @test all(r -> r isa BOLSIGExplicitRun, input.runs)
    @test size(input.runs[1].values) == (7, 2)
    @test size(input.runs[2].values) == (1, 2)
    @test size(input.runs[3].values) == (2, 2)
    # total rows (7+1+2=10) matches the already-verified example2.dat output size
    @test sum(size(r.values, 1) for r in input.runs) == 10
    @test input.runs[1].values[1, :] == [100.0, 0.0]

    @test input.save.format == 2   # "Combined" -> the condensed R#/A#/C# reader
end

@testset "read_bolsig_input: example3 (RUN2D, VAR gas fraction, no SAVERESULTS)" begin
    input = read_bolsig_input(joinpath(DATA, "example3_input.dat"))

    @test input.conditions.reduced_field == 10.0   # NOT VAR: RUN2D's 1st variable is chosen by code
    @test input.conditions.gas_fractions == [1.0, nothing]   # 2nd fraction is VAR: RUN2D's 2nd variable

    @test length(input.runs) == 1
    run = input.runs[1]
    @test run isa BOLSIGRun2D
    @test run.variable1 == ReducedFieldVar
    @test (run.min1, run.max1, run.num1, run.type1) == (0.1, 100.0, 20, QuadraticSeries)
    @test (run.min2, run.max2, run.num2, run.type2) == (0.0, 1.0, 5, LinearSeries)
    @test run.output_file == "example3.dat"

    @test input.save === nothing   # RUN2D writes its own file, no SAVERESULTS block
end

@testset "read_bolsig_input: example4 (multi-species, & continuation, bare RUN)" begin
    input = read_bolsig_input(joinpath(DATA, "example4_input.dat"))

    @test length(input.collisions) == 1
    @test input.collisions[1].species == ["Ar", "He", "N2", "O2"]

    # fractions split across 3 physical lines via "&" continuation
    @test input.conditions.gas_fractions == [0.7, 0.2, 0.08, 0.02]

    @test length(input.runs) == 1
    @test input.runs[1] isa BOLSIGFixedRun

    @test input.save.format == 1   # "Run by run" -> the single-run-report reader
end

@testset "round-trip: parse -> write -> reparse matches for all 4 examples" begin
    for i in 1:4
        input1 = read_bolsig_input(joinpath(DATA, "example$(i)_input.dat"))
        tmp = tempname()
        write_bolsig_input(tmp, input1)
        input2 = read_bolsig_input(tmp)
        @test _deep_eq(input1, input2)
    end
end

@testset "write_bolsig_input validates VAR/run-spec consistency" begin
    collisions = [BOLSIGReadCollisions("x.txt", ["Ar"], true)]

    # ExplicitRun column count must match the number of VAR fields (here: 0)
    bad_conditions = BOLSIGConditions(; reduced_field=10.0, gas_fractions=[1.0])
    bad_input = BOLSIGInput(;
        collisions, conditions=bad_conditions,
        runs=[BOLSIGExplicitRun([1.0 2.0])],
        save=BOLSIGSaveResults(file="out.dat", format=2),
    )
    @test_throws ErrorException write_bolsig_input(tempname(), bad_input)

    # Run2D needs exactly one VAR field (here: 0)
    bad_input2 = BOLSIGInput(;
        collisions, conditions=bad_conditions,
        runs=[BOLSIGRun2D(ReducedFieldVar, 0.1, 100.0, 10, LinearSeries, 0.0, 1.0, 5, LinearSeries, "out2d.dat")],
        save=nothing,
    )
    @test_throws ErrorException write_bolsig_input(tempname(), bad_input2)

    # non-RUN2D runs need a `save`
    bad_input3 = BOLSIGInput(;
        collisions, conditions=bad_conditions,
        runs=[BOLSIGFixedRun()],
        save=nothing,
    )
    @test_throws ErrorException write_bolsig_input(tempname(), bad_input3)

    # a consistent input writes without error
    good_conditions = BOLSIGConditions(; reduced_field=nothing, gas_fractions=[1.0])
    good_input = BOLSIGInput(;
        collisions, conditions=good_conditions,
        runs=[BOLSIGExplicitRun(reshape([1.0, 2.0], 2, 1))],
        save=BOLSIGSaveResults(file="out.dat", format=2),
    )
    tmp = tempname()
    @test write_bolsig_input(tmp, good_input) == tmp
    @test _deep_eq(read_bolsig_input(tmp), good_input)
end
