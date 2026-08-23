using BoltzmannSolvers
using Test

# The committed N2 swarm-analysis setup shipped with LoKI-B.
input = read_loki_input(joinpath(@__DIR__, "data", "loki_n2_swarm_setup.in"))

@testset "working conditions" begin
    @test length(input.reduced_field) == 100
    @test input.reduced_field[1] ≈ 1e-3
    @test input.reduced_field[end] ≈ 1e3
    @test length(input.electron_temperature) == 100
    @test input.electron_temperature[1] ≈ 0.03
    @test input.excitation_frequency == 0.0
    @test input.gas_pressure ≈ 133.32
    @test input.gas_temperature == 300.0
    @test input.electron_density ≈ 1e19
end

@testset "electron kinetics" begin
    @test input.eedf_type == "boltzmann"
    @test input.shape_parameter === nothing   # commented out in the file
    @test input.ionization_operator == "usingSDCS"
    @test input.growth_model == "spatial"
    @test input.include_ee == false
    @test input.lxcat_files ==
          ["Nitrogen/N2_LXCat.txt", "Nitrogen/N2_rot_LXCat.txt"]
    @test isempty(input.lxcat_extra_files)    # commented out
    @test isempty(input.car_gases)            # commented out
    @test input.fractions == ["N2" => 1.0]
end

@testset "state populations" begin
    # `N2(X) = 1.0`, a file reference, and a population function — kept verbatim
    @test length(input.populations) == 3
    @test input.populations[1] == ("N2(X)" => 1.0)
    @test input.populations[2] == "Nitrogen/N2_vibpop.txt"
    @test input.populations[3] ==
          ("N2(X,v=0,J=*)" => "boltzmannPopulation@gasTemperature")
end

@testset "numerics" begin
    @test input.max_energy == 1.0
    @test input.cell_number == 1000
    @test input.smart_grid isa LoKISmartGrid
    @test input.smart_grid.min_eedf_decay == 20.0
    @test input.smart_grid.max_eedf_decay == 25.0
    @test input.smart_grid.update_factor == 0.05
    @test input.mixing_parameter == 0.7
    @test input.max_eedf_rel_error == 1e-9
end

@testset "raw tree" begin
    @test input.raw["chemistry"]["isOn"] == false
    @test input.raw["output"]["dataFormat"] == "hdf5+txt"
    @test input.raw["workingConditions"]["totalSccmOutFlow"] == "ensureIsobaric"
    # state energies/statistical weights are not modeled but stay reachable
    energies = input.raw["electronKinetics"]["stateProperties"]["energy"]
    @test ("N2(X,v=*)" => "harmonicOscillatorEnergy") in energies
end

@testset "defaults on a minimal file" begin
    path = joinpath(mktempdir(), "minimal.in")
    write(path, """
        workingConditions:
          reducedField: 40
          gasTemperature: 400
        electronKinetics:
          isOn: true
          LXCatFiles: Ar.txt
          gasProperties:
            fraction:
              - Ar = 1
        """)
    m = read_loki_input(path)
    @test m.reduced_field == [40.0]
    @test m.gas_temperature == 400.0
    @test m.eedf_type == "boltzmann"
    @test m.growth_model == "temporal"
    @test m.lxcat_files == ["Ar.txt"]
    @test m.fractions == ["Ar" => 1.0]
    @test m.smart_grid === nothing
    @test isempty(m.populations)
end
