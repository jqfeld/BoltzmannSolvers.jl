using BoltzmannSolvers
using DataFrames
using Test

const DATA = joinpath(@__DIR__, "data", "loki", "simulation_1")
const H5_PATH = joinpath(DATA, "simulation_1.h5")

@testset "LoKI-B lookup tables (with power balance)" begin
    df = load_dataframe(LoKI(), DATA)

    @test size(df) == (100, 131)
    @test allunique(names(df))

    @test df.reduced_field[1] ≈ 0.001
    @test df.reduced_field[2] ≈ 1.14975699539774e-03
    @test df.reduced_diffusion_coef[1] ≈ 9.68802449186642e23
    @test df.reduced_mobility[1] ≈ 3.90470412651995e25
    @test df.drift_velocity[1] ≈ 39.0470412651994
    @test df.mean_energy[1] ≈ 3.72047005392004e-02
    @test df.characteristic_energy[1] ≈ 2.48111615578434e-02
    @test df.electron_temperature[1] ≈ 2.48031336928003e-02

    # Effective collision (irreversible, "->") — normalized as "Effective(lhs)".
    @test df[!, Symbol("Effective(e+N2(X))")][1] ≈ 4.27413027133661e-15

    # Reversible ("<->") vibrational collision produces both an ine (forward)
    # and sup (backward) column.
    @test df[!, "e+N2(X,v=0)-->e+N2(X,v=1)"][1] ≈ 5.45981930854129e-23
    @test df[!, "e+N2(X,v=0)<--e+N2(X,v=1)"][1] ≈ 1.08289127873106e-17

    # Reversible rotational collisions parse the same way (a distinct "type"
    # string in the header comment, not specially handled beyond "Effective").
    @test "e+N2(X,v=0,J=0)-->e+N2(X,v=0,J=2)" in names(df)
    @test "e+N2(X,v=0,J=0)<--e+N2(X,v=0,J=2)" in names(df)

    # Irreversible ionization ("->") gives only a forward column.
    @test "e+N2(X)-->2e+N2(+,X)" in names(df)
    @test !("e+N2(X)<--2e+N2(+,X)" in names(df))

    # lookUpTablePower.txt (per-channel electron energy gain/loss balance) —
    # newer LoKI-B output, joined in when present.
    @test df[!, "PowerField(eVm^3s^-1)"][1] ≈ 3.90470412651995e-23
    @test df[!, "PwrElaGain(eVm^3s^-1)"][1] ≈ 9.98219297139508e-21

    # "RelPwrBalance" is written with a literal trailing "%" glued onto the
    # number (e.g. "6.95053583371272e-13%") — must still parse to Float64.
    @test eltype(df[!, "RelPwrBalance"]) == Float64
    @test df[!, "RelPwrBalance"][1] ≈ 6.95053583371272e-13
end

@testset "LoKI-B lookup tables (no power balance file — backward compatible)" begin
    dir = mktempdir()
    cp(joinpath(DATA, "lookUpTableSwarm.txt"), joinpath(dir, "lookUpTableSwarm.txt"))
    cp(joinpath(DATA, "lookUpTableRateCoeff.txt"), joinpath(dir, "lookUpTableRateCoeff.txt"))

    df = load_dataframe(LoKI(), dir)
    @test size(df) == (100, 110)
    @test !any(n -> occursin("Pwr", n) || occursin("PowerField", n), names(df))
end

@testset "create_interpolation works on parsed LoKI-B output" begin
    df = load_dataframe(LoKI(), DATA)
    itp = create_interpolation(df, :mean_energy, :reduced_field)
    @test itp(0.001) ≈ 3.72047005392004e-02
end

# Must run *before* HDF5 is `using`-loaded anywhere below (package extensions
# activate process-wide, not per-SafeTestset-module, so this is the only
# point in the whole suite where the un-triggered stub's behavior can still
# be observed).
@testset "LoKI-B HDF5 output without HDF5.jl loaded gives a clear error" begin
    err = try
        load_dataframe(LoKI(), H5_PATH)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("HDF5.jl", err.msg)
end

if Base.identify_package("HDF5") !== nothing
    using HDF5

    @testset "LoKI-B HDF5 output (simulation_1.h5)" begin
        df = load_dataframe(LoKI(), H5_PATH)

        @test size(df) == (100, 131)
        @test allunique(names(df))

        @test df.reduced_field[1] ≈ 0.001
        @test df.mean_energy[1] ≈ 3.72047005392004e-02 rtol=1e-10
        @test df[!, Symbol("Effective(e+N2(X))")][1] ≈ 4.27413027133661e-15 rtol=1e-10
        @test df[!, "e+N2(X,v=0)-->e+N2(X,v=1)"][1] ≈ 5.45981930854129e-23 rtol=1e-10
        @test df[!, "e+N2(X,v=0)<--e+N2(X,v=1)"][1] ≈ 1.08289127873106e-17 rtol=1e-10
        @test "e+N2(X,v=0,J=0)-->e+N2(X,v=0,J=2)" in names(df)
        @test "e+N2(X)-->2e+N2(+,X)" in names(df)
        @test !("e+N2(X)<--2e+N2(+,X)" in names(df))

        @test df[!, "PowerField(eVm^3s^-1)"][1] ≈ 3.90470412651995e-23 rtol=1e-10
        @test df[!, "PwrElaGain(eVm^3s^-1)"][1] ≈ 9.98219297139508e-21 rtol=1e-10
        @test df[!, "RelPwrBalance"][1] ≈ 6.95053583371272e-13 rtol=1e-8
    end

    @testset "LoKI-B txt and HDF5 paths agree on the same simulation" begin
        df_txt = load_dataframe(LoKI(), DATA)
        df_h5 = load_dataframe(LoKI(), H5_PATH)

        @test Set(names(df_txt)) == Set(names(df_h5))
        for c in names(df_txt)
            @test df_txt[!, c] ≈ df_h5[!, c] rtol=1e-9
        end
    end
else
    @warn "Skipping LoKI-B HDF5 tests — HDF5 not found in current environment"
end
