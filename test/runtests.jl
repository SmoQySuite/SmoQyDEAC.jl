using SmoQyDEAC
using Test
using Random


@testset "SmoQyDEAC.jl" begin
    # Test gemmavx!()
    A = rand(Float64,(101,20))
    B = rand(Float64,(20,45))
    C = Array{Float64}(undef,(size(A,1),size(B,2)))
    SmoQyDEAC.gemmavx!(C,A,B)
    @test all(A*B .≈  C)

end
