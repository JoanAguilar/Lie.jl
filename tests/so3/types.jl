@testset begin "Lie Group Constructor"
    x = [1 0  0;
         0 0 -1;
         0 1  0]
    @test isa(RotationMatrix(x), RotationMatrix{eltype(x)})
    @test isa(RotationMatrix{Float32}(x), RotationMatrix{Float32})
    @test isa(RotationMatrix(rand(Float32, 4, 2), checks=false), RotationMatrix{Float32})
    @test_throws ArgumentError RotationMatrix(rand(4, 2))
    @test_throws ArgumentError RotationMatrix(rand(3, 3))
    @test_throws ArgumentError RotationMatrix([0.1 0 0;
                                               0 0.1 0;
                                               0 0 0.1])
end


@testset begin "Lie Algebra Constructor"
    ω = [1; 2; 3]
    @test isa(VectorSO3Algebra(ω), VectorSO3Algebra{eltype(ω)})
    @test isa(VectorSO3Algebra{Float32}(ω), VectorSO3Algebra{Float32})
    @test isa(VectorSO3Algebra(rand(Float32, 4), checks=false), VectorSO3Algebra{Float32})
    @test_throws ArgumentError VectorSO3Algebra(rand(4))
end
