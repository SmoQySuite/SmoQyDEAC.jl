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

    # Uniform grids preserve the original rectangle-rule weights, while
    # nonuniform grids use their local energy-bin widths.
    @test SmoQyDEAC.omega_weights([0.0, 0.5, 1.0]) == [0.5, 0.5, 0.5]
    @test SmoQyDEAC.omega_weights([0.0, 0.25, 1.0]) == [0.25, 0.5, 0.75]
    @test_throws ArgumentError SmoQyDEAC.omega_weights([0.0, 1.0, 0.5])

    params = SmoQyDEAC.DEACParameters(
        1.0, [0.0, 1.0], [0.0, 0.25, 1.0], "time_bosonic",
        "", "", 2, 1, 8, 1, 0.9, 0.1, 0.9, 0.1, 0.9, 1,
    )
    K = SmoQyDEAC.generate_K(params)
    @test K[1, :] ≈ [0.25, 0.5, 0.75]

    covariance = [0.04 0.0; 0.0 0.01]
    cov_params = SmoQyDEAC.DEACParameters(
        1.0, [0.0, 1.0], [-1.0, 0.0, 1.0], "time_bosonic",
        "", "", 2, 1, 8, 1, 0.9, 0.1, 0.9, 0.1, 0.9, 1,
    )
    cov_K = SmoQyDEAC.generate_K(cov_params)
    W, Kp, corr_p, eigenvalues = SmoQyDEAC.calculate_fit_matrices(
        ([0.5, 0.25], covariance), cov_K, false, false, cov_params, 1e-8,
    )
    @test sort(W) ≈ [0.5 / (2 * 0.04), 0.5 / (2 * 0.01)]
    @test size(Kp) == size(cov_K)
    @test length(corr_p) == 2
    @test eigenvalues ≈ [0.01, 0.04]

    @test_throws DimensionMismatch SmoQyDEAC.DEAC_Cov(
        [0.5, 0.25], ones(3, 3), 1.0, [0.0, 1.0], [-1.0, 0.0, 1.0],
        "time_bosonic", 2, 1, "", ""; number_of_generations=1,
    )

    Greens = [0.5,0.4,0.3,0.2,0.1,0.05,0.1,0.2,0.3,0.4,0.5]
    Greens_std = zeros(Float64,11) .+ 0.02

    @test typeof( SmoQyDEAC.DEAC_Std(Greens,Greens_std,1.0,collect(LinRange(0.0,1.0,11)),collect(LinRange(-10.0,10.0,401)),"time_fermionic",2,10,"x.jld2","DEAC_checkpoint.jld2")) == Dict{String,Any}
    rm("x.jld2")


end
