include("../Lie.jl")

using LinearAlgebra: norm
using Test: @test, @testset, @test_throws

using .Lie

@testset "SO3" begin
    include("so3/so3.jl")
end
