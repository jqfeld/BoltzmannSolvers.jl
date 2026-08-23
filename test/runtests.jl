using SafeTestsets

@safetestset "LoKI solver" begin
    include("test_loki.jl")
end

@safetestset "LoKI input files" begin
    include("test_loki_input.jl")
end

@safetestset "BOLSIG solver" begin
    include("test_bolsig.jl")
end

@safetestset "BOLSIG input files" begin
    include("test_bolsig_input.jl")
end

@safetestset "run_solver(::BOLSIGInput)" begin
    include("test_bolsig_run.jl")
end

@safetestset "run_solver(::MultiBoltInput)" begin
    include("test_multibolt_run.jl")
end
