using SmoQyDEAC
using Test
using Random


@testset "SmoQyDEAC.jl" begin
    # Test gemmSIMD!()
    A = rand(Float64,(101,20))
    B = rand(Float64,(20,45))
    C = Array{Float64}(undef,(size(A,1),size(B,2)))
    SmoQyDEAC.gemmSIMD!(C,A,B)
    @test all(A*B .≈  C)
    Greens = [0.5,0.4,0.3,0.2,0.1,0.05,0.1,0.2,0.3,0.4,0.5]
    Greens_std = zeros(Float64,11) .+ 0.02

    @test typeof( SmoQyDEAC.DEAC_Std(Greens,Greens_std,1.0,collect(LinRange(0.0,1.0,11)),collect(LinRange(-10.0,10.0,401)),"time_fermionic",2,10,"x.jld2",".")) == Dict{String,Any}
    rm("x.jld2")


end
