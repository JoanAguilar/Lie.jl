@testset "One Element" begin
    x = RotationMatrix([0 0 1; 0 1 0; -1 0 0])

    @test one(RotationMatrix) * x == x
    @test x * one(RotationMatrix) == x
    @test one(RotationMatrix) * one(RotationMatrix) == one(RotationMatrix)
    @test one(RotationMatrix) == inv(one(RotationMatrix))
    @test log(one(RotationMatrix), VectorSO3Algebra) ≈ zero(VectorSO3Algebra)

    @test one(x) * x ≈ x
    @test x * one(x) ≈ x
    @test one(x) * one(x) == one(x)
    @test one(x) == inv(one(x))
    @test log(one(x), VectorSO3Algebra) ≈ zero(VectorSO3Algebra)
end

@testset "Zero Element" begin
    x = VectorSO3Algebra([1; 2; 3])

    @test zero(VectorSO3Algebra) + x == x
    @test x + zero(VectorSO3Algebra) == x
    @test zero(VectorSO3Algebra) + zero(VectorSO3Algebra) == zero(VectorSO3Algebra)
    @test zero(VectorSO3Algebra) == -zero(VectorSO3Algebra)
    @test exp(zero(VectorSO3Algebra), RotationMatrix) ≈ one(RotationMatrix)

    @test zero(x) + x == x
    @test x + zero(x) == x
    @test zero(x) + zero(x) == zero(x)
    @test zero(x) == -zero(x)
    @test exp(zero(x), RotationMatrix) ≈ one(RotationMatrix)
end
