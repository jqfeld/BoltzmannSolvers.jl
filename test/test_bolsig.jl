using BoltzmannSolvers
using DataFrames
using Test

const DATA = joinpath(@__DIR__, "data", "bolsig")

@testset "BOLSIG+ verbose format (example1.dat)" begin
    df = load_dataframe(BOLSIG(), joinpath(DATA, "example1.dat"))

    @test size(df) == (40, 28)
    @test allunique(names(df))

    @test df.reduced_field[1] ≈ 0.1
    @test df.reduced_field[2] ≈ 0.125893
    @test df.mean_energy[1] ≈ 0.0923585
    @test df.reduced_mobility[1] ≈ 1.52233e25 rtol=1e-5

    # Rate coefficients, normalized reaction names.
    @test df[!, Symbol("Effective(Ar)")][1] ≈ 1.538668e-15

    # BOLSIG+ echoes "Maximum energy" as two distinct blocks under the exact
    # same title (a real quirk, see src/solvers/bolsig.jl) — both must survive
    # as separate, disambiguated columns rather than one overwriting the other.
    @test df[!, "Maximum energy"][1] ≈ 200.0
    @test df[!, "Maximum energy (dup2)"][1] ≈ 1.03464

    # BOLSIG+ drops the "E" in an exponent when a fixed-width field would
    # otherwise overflow with a 3-digit exponent, e.g. "-0.852184-208" instead
    # of "-0.852184E-208" — must still parse to the correct (very small) value.
    @test df[!, "Inelastic power loss /N (eV m3/s)"][1] ≈ -8.52184e-209
end

@testset "BOLSIG+ condensed format (example2.dat)" begin
    df = load_dataframe(BOLSIG(), joinpath(DATA, "example2.dat"))

    @test size(df) == (10, 27)
    @test allunique(names(df))

    @test df.reduced_field[1] ≈ 100.0
    @test df.reduced_field[9] ≈ 1000.00
    @test df.reduced_field[10] ≈ 10000.0
    @test df.mean_energy[1] ≈ 6.94903
    @test df.reduced_mobility[1] ≈ 7.93112e23 rtol=1e-5

    @test df[!, "Re/perp mobility *N (1/m/V/s)"][1] ≈ 0.0
    @test df[!, "Re/perp mobility *N (1/m/V/s)"][2] ≈ 7.93111e23 rtol=1e-5

    @test df[!, Symbol("Effective(Ar)")][1] ≈ 1.56811e-13 rtol=1e-5

    # Same dropped-"E" exponent quirk as example1.dat, here in a data row of
    # the condensed table instead of a wrapped verbose-format block.
    @test df[!, "Inelastic power loss /N (eV m3/s)"][7] ≈ -1.58317e-207
end

@testset "BOLSIG+ single-run report format (example4.dat)" begin
    df = load_dataframe(BOLSIG(), joinpath(DATA, "example4.dat"))

    @test size(df) == (1, 90)
    @test allunique(names(df))

    @test df.mean_energy[1] ≈ 1.35522
    @test df.reduced_mobility[1] ≈ 4.3799e24 rtol=1e-5
    @test df[!, Symbol("Effective(Ar)")][1] ≈ 1.50155e-14 rtol=1e-5
    @test df[!, Symbol("Elastic(N2)")][1] ≈ 7.92173e-14 rtol=1e-5

    # BOLSIG+ itself gives two distinct N2 collisions (C9, C10) the identical
    # truncated description "N2    Excitation    0.29 eV" — both values must
    # survive as separate, disambiguated columns, not silently collapse to one.
    @test df[!, Symbol("Excitation(N2,0.29eV)")][1] ≈ 1.57234e-16 rtol=1e-5
    @test df[!, "Excitation(N2,0.29eV) (dup2)"][1] ≈ 3.23263e-15 rtol=1e-5
end

@testset "BOLSIG+ 2D scan format (example3.dat) is not supported" begin
    # A genuine 2D parametric scan (E/N rows × gas-mixture-fraction columns) —
    # each quantity is a matrix, not a single column, so it doesn't fit the
    # "one row per condition" DataFrame shape the other three formats share.
    # Should fail clearly rather than silently misparse.
    @test_throws ErrorException load_dataframe(BOLSIG(), joinpath(DATA, "example3.dat"))
end

@testset "create_interpolation works on parsed BOLSIG+ output" begin
    df = load_dataframe(BOLSIG(), joinpath(DATA, "example1.dat"))
    itp = create_interpolation(df, :mean_energy, :reduced_field)
    @test itp(0.1) ≈ 0.0923585
    @test itp(0.125893) ≈ 0.108603
end
