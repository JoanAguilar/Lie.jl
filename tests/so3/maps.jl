@testset "Exponential Map" begin
    @test exp(zero(VectorSO3Algebra), RotationMatrix) ≈ one(RotationMatrix)
    @test exp(VectorSO3Algebra([0.5 * π; 0; 0]), RotationMatrix) ≈
        RotationMatrix([1 0 0; 0 0 -1; 0 1 0])
    @test exp(VectorSO3Algebra([0; 0.5 * π; 0]), RotationMatrix) ≈
        RotationMatrix([0 0 1; 0 1 0; -1 0 0])
    @test exp(VectorSO3Algebra([0; 0; 0.5 * π]), RotationMatrix) ≈
        RotationMatrix([0 -1 0; 1 0 0; 0 0 1])
end

@testset "Logarithmic Map" begin
    @test log(one(RotationMatrix), VectorSO3Algebra) ≈ zero(VectorSO3Algebra)
    @test log(RotationMatrix([1 0 0; 0 0 -1; 0 1 0]), VectorSO3Algebra) ≈
        VectorSO3Algebra([0.5 * π; 0; 0])
    @test log(RotationMatrix([0 0 1; 0 1 0; -1 0 0]), VectorSO3Algebra) ≈
        VectorSO3Algebra([0; 0.5 * π; 0])
    @test log(RotationMatrix([0 -1 0; 1 0 0; 0 0 1]), VectorSO3Algebra) ≈
        VectorSO3Algebra([0; 0; 0.5 * π])
end
