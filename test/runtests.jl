using SafeTestsets

@safetestset "LoKI solver" begin
    include("test_loki.jl")
end

@safetestset "BOLSIG solver" begin
    include("test_bolsig.jl")
end

@safetestset "BOLSIG input files" begin
    include("test_bolsig_input.jl")
end

@safetestset "run_bolsig" begin
    include("test_bolsig_run.jl")
end

@safetestset "run_multibolt" begin
    include("test_multibolt_run.jl")
end
