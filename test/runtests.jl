using SafeTestsets

@safetestset "LoKI solver" begin
    include("test_loki.jl")
end

@safetestset "BOLSIG solver" begin
    include("test_bolsig.jl")
end
