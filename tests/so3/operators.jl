@testset "Lie Group Multiplication" begin
    x = RotationMatrix([1 0  0;
                        0 0 -1;
                        0 1  0])
    y = RotationMatrix([0 0 1;
                        0 1 0;
                       -1 0 0])
    z = RotationMatrix([0 -1 0;
                        1  0 0;
                        0  0 1])
    xy = RotationMatrix([0 0 1;
                         1 0 0;
                         0 1 0])
    xz = RotationMatrix([0 -1  0;
                         0  0 -1;
                         1  0  0])
    yz = RotationMatrix([0 0 1;
                         1 0 0;
                         0 1 0])
    xyz = RotationMatrix([0  0 1;
                          0 -1 0;
                          1  0 0])

    @test one(RotationMatrix) * one(RotationMatrix) ≈ one(RotationMatrix)

    @test x * one(RotationMatrix) ≈ x
    @test y * one(RotationMatrix) ≈ y
    @test z * one(RotationMatrix) ≈ z
    @test one(RotationMatrix) * x ≈ x
    @test one(RotationMatrix) * y ≈ y
    @test one(RotationMatrix) * z ≈ z

    @test x * y ≈ xy
    @test x * z ≈ xz
    @test y * z ≈ yz
    @test x * y * z ≈ xyz
    @test *(x, y, z) ≈ xyz
end

@testset "Lie Group Inverse" begin
    x = RotationMatrix([1 0  0;
                        0 0 -1;
                        0 1  0])
    y = RotationMatrix([0 0 1;
                        0 1 0;
                       -1 0 0])
    z = RotationMatrix([0 -1 0;
                        1  0 0;
                        0  0 1])

    @test one(RotationMatrix) * inv(one(RotationMatrix)) ≈ one(RotationMatrix)
    @test x * inv(x) ≈ one(RotationMatrix)
    @test y * inv(y) ≈ one(RotationMatrix)
    @test z * inv(z) ≈ one(RotationMatrix)
    @test inv(x) * x ≈ one(RotationMatrix)
    @test inv(y) * y ≈ one(RotationMatrix)
    @test inv(z) * z ≈ one(RotationMatrix)
    @test x * y * z * inv(z) * inv(y) * inv(x) ≈ one(RotationMatrix)
end

@testset "Lie Algebra Summation" begin
    @test zero(VectorSO3Algebra) + zero(VectorSO3Algebra) == zero(VectorSO3Algebra)

    n_tests = 10

    for n = 1:n_tests
        w = rand(Float64, 3)
        x = rand(Float64, 3)
        y = rand(Float64, 3)
        ω = VectorSO3Algebra(w)
        χ = VectorSO3Algebra(x)
        η = VectorSO3Algebra(y)

        @test ω + zero(VectorSO3Algebra) ≈ ω
        @test χ + zero(VectorSO3Algebra) ≈ χ
        @test η + zero(VectorSO3Algebra) ≈ η
        @test zero(VectorSO3Algebra) + ω ≈ ω
        @test zero(VectorSO3Algebra) + χ ≈ χ
        @test zero(VectorSO3Algebra) + η ≈ η

        @test ω + χ ≈ VectorSO3Algebra(w + x)
        @test ω + η ≈ VectorSO3Algebra(w + y)
        @test χ + η ≈ VectorSO3Algebra(x + y)
        @test ω + χ + η ≈ VectorSO3Algebra(w + x + y)
        @test +(ω, χ, η) ≈ VectorSO3Algebra(w + x + y)
    end
end

@testset "Lie Algebra Subtraction" begin
    n_tests = 10

    for n = 1:n_tests
        w = rand(Float64, 3)
        x = rand(Float64, 3)
        y = rand(Float64, 3)
        ω = VectorSO3Algebra(w)
        χ = VectorSO3Algebra(x)
        η = VectorSO3Algebra(y)

        @test ω - χ ≈ VectorSO3Algebra(w - x)
        @test ω - η ≈ VectorSO3Algebra(w - y)
        @test χ - η ≈ VectorSO3Algebra(x - y)
        @test ω - ω ≈ zero(VectorSO3Algebra)
        @test η - η ≈ zero(VectorSO3Algebra)
        @test χ - χ ≈ zero(VectorSO3Algebra)
    end
end

@testset "Lie Algebra Negation" begin
    n_tests = 10

    for n = 1:n_tests
        w = rand(Float64, 3)
        x = rand(Float64, 3)
        y = rand(Float64, 3)
        ω = VectorSO3Algebra(w)
        χ = VectorSO3Algebra(x)
        η = VectorSO3Algebra(y)

        @test ω + (-ω) ≈ zero(VectorSO3Algebra)
        @test η + (-η) ≈ zero(VectorSO3Algebra)
        @test χ + (-χ) ≈ zero(VectorSO3Algebra)
    end
end

@testset "Lie Algebra Scalar Multiplication" begin
    n_tests = 10

    for n = 1:n_tests
        a = rand(Float64)
        b = rand(Float64)
        c = rand(Float64)
        w = rand(Float64, 3)
        x = rand(Float64, 3)
        y = rand(Float64, 3)
        ω = VectorSO3Algebra(w)
        χ = VectorSO3Algebra(x)
        η = VectorSO3Algebra(y)

        @test a * ω ≈ VectorSO3Algebra(a * w)
        @test b * χ ≈ VectorSO3Algebra(b * x)
        @test c * η ≈ VectorSO3Algebra(c * y)
        @test ω * a ≈ VectorSO3Algebra(a * w)
        @test χ * b ≈ VectorSO3Algebra(b * x)
        @test η * c ≈ VectorSO3Algebra(c * y)
    end
end

@testset "Lie Group-Algebra Multiplication" begin
    x = RotationMatrix([1 0  0;
                        0 0 -1;
                        0 1  0])
    y = RotationMatrix([0 0 1;
                        0 1 0;
                       -1 0 0])
    z = RotationMatrix([0 -1 0;
                        1  0 0;
                        0  0 1])
    ω = VectorSO3Algebra([1; 2; 3])
    χ = VectorSO3Algebra([2; 3; 1])
    η = VectorSO3Algebra([3; 1; 2])

    @test one(RotationMatrix) * zero(VectorSO3Algebra) ≈ zero(VectorSO3Algebra)
    @test one(RotationMatrix) * ω ≈ ω
    @test one(RotationMatrix) * χ ≈ χ
    @test one(RotationMatrix) * η ≈ η
    @test x * zero(VectorSO3Algebra) ≈ zero(VectorSO3Algebra)
    @test y * zero(VectorSO3Algebra) ≈ zero(VectorSO3Algebra)
    @test z * zero(VectorSO3Algebra) ≈ zero(VectorSO3Algebra)

    @test x * ω ≈ VectorSO3Algebra([1; -3; 2])
    @test x * χ ≈ VectorSO3Algebra([2; -1; 3])
    @test x * η ≈ VectorSO3Algebra([3; -2; 1])
    @test y * ω ≈ VectorSO3Algebra([3; 2; -1])
    @test y * χ ≈ VectorSO3Algebra([1; 3; -2])
    @test y * η ≈ VectorSO3Algebra([2; 1; -3])
    @test z * ω ≈ VectorSO3Algebra([-2; 1; 3])
    @test z * χ ≈ VectorSO3Algebra([-3; 2; 1])
    @test z * η ≈ VectorSO3Algebra([-1; 3; 2])
end
