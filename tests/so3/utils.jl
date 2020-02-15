@testset begin "Conversion"
    n_tests = 10
    for n = 1:n_tests
        x = rand(Float32, 3, 3)
        ω = rand(Float32, 3)
        T = Float64
        @test convert(Array, RotationMatrix(x; checks=false)) == x
        @test eltype(convert(Array, RotationMatrix(x; checks=false))) == eltype(x)
        @test eltype(convert(Array{T}, RotationMatrix(x; checks=false))) == T
        @test eltype(convert(Array{T, 2}, RotationMatrix(x; checks=false))) == T
        @test convert(Array, VectorSO3Algebra(ω; checks=false)) == ω
        @test eltype(convert(Array, VectorSO3Algebra(ω; checks=false))) == eltype(ω)
        @test eltype(convert(Array{T}, VectorSO3Algebra(ω; checks=false))) == T
        @test eltype(convert(Array{T, 1}, VectorSO3Algebra(ω; checks=false))) == T
    end
end


@testset begin "Angle"
    x = RotationMatrix([1 0  0;
                        0 0 -1;
                        0 1  0])
    y = RotationMatrix([0 0 1;
                        0 1 0;
                       -1 0 0])
    z = RotationMatrix([0 -1 0;
                        1  0 0;
                        0  0 1])

    ω = VectorSO3Algebra([1; 1; 1])
    χ = VectorSO3Algebra([0.5π; 0; 0])
    η = VectorSO3Algebra([0; 0.5π; 0])
    ζ = VectorSO3Algebra([0; 0; 0.5π])

    @test angle(one(RotationMatrix)) ≈ 0
    @test angle(x) ≈ 0.5π
    @test angle(y) ≈ 0.5π
    @test angle(z) ≈ 0.5π

    @test angle(zero(VectorSO3Algebra)) ≈ 0
    @test angle(ω) ≈ norm([1; 1; 1])
    @test angle(χ) ≈ 0.5π
    @test angle(η) ≈ 0.5π
    @test angle(ζ) ≈ 0.5π

    n_tests = 1
    for n = 1:n_tests
        υ = VectorSO3Algebra(rand(3))
        @test angle(exp(υ, RotationMatrix)) ≈ angle(υ)
    end
end


@testset begin "Axis"
    x = RotationMatrix([1 0  0;
                        0 0 -1;
                        0 1  0])
    y = RotationMatrix([0 0 1;
                        0 1 0;
                       -1 0 0])
    z = RotationMatrix([0 -1 0;
                        1  0 0;
                        0  0 1])
    x2 = RotationMatrix([1  0  0;
                         0 -1  0;
                         0  0 -1])
    y2 = RotationMatrix([-1 0  0;
                          0 1  0;
                          0 0 -1])
    z2 = RotationMatrix([-1  0 0;
                          0 -1 0;
                          0  0 1])

    ω = VectorSO3Algebra([1; 1; 1])
    χ = VectorSO3Algebra([0.5π; 0; 0])
    η = VectorSO3Algebra([0; 0.5π; 0])
    ζ = VectorSO3Algebra([0; 0; 0.5π])

    @test all(map(isnan, axis(one(RotationMatrix))))
    @test axis(x) ≈ [1; 0; 0]
    @test axis(y) ≈ [0; 1; 0]
    @test axis(z) ≈ [0; 0; 1]
    @test all(map(isnan, axis(x2)))
    @test all(map(isnan, axis(y2)))
    @test all(map(isnan, axis(z2)))

    @test all(map(isnan, axis(zero(VectorSO3Algebra))))
    @test axis(ω) ≈ (1 / √3) * [1; 1; 1]
    @test axis(χ) ≈ [1; 0; 0]
    @test axis(η) ≈ [0; 1; 0]
    @test axis(ζ) ≈ [0; 0; 1]

    n_tests = 1
    for n = 1:n_tests
        υ = VectorSO3Algebra(rand(3))
        @test axis(exp(υ, RotationMatrix)) ≈ axis(υ)
    end
end


@testset begin "Conversion"
    x = RotationMatrix([1 0  0;
                        0 0 -1;
                        0 1  0])
    y = RotationMatrix([0 0 1;
                        0 1 0;
                       -1 0 0])
    z = RotationMatrix([0 -1 0;
                        1  0 0;
                        0  0 1])

    χ = VectorSO3Algebra([0.5π; 0; 0])
    η = VectorSO3Algebra([0; 0.5π; 0])
    ζ = VectorSO3Algebra([0; 0; 0.5π])

    @test one(RotationMatrix) ≈ one(RotationMatrix)
    @test x ≈ x
    @test y ≈ y
    @test z ≈ z
    @test !(x ≈ y)
    @test !(x ≈ z)
    @test !(y ≈ z)

    @test one(RotationMatrix) == one(RotationMatrix)
    @test x == x
    @test y == y
    @test z == z
    @test !(x == y)
    @test !(x == z)
    @test !(y == z)

    @test isequal(one(RotationMatrix), one(RotationMatrix))
    @test isequal(x, x)
    @test isequal(y, y)
    @test isequal(z, z)
    @test !isequal(x, y)
    @test !isequal(x, z)
    @test !isequal(y, z)

    @test zero(VectorSO3Algebra) ≈ zero(VectorSO3Algebra)
    @test χ ≈ χ
    @test η ≈ η
    @test ζ ≈ ζ
    @test !(χ ≈ η)
    @test !(χ ≈ ζ)
    @test !(η ≈ ζ)

    @test zero(VectorSO3Algebra) == zero(VectorSO3Algebra)
    @test χ == χ
    @test η == η
    @test ζ == ζ
    @test !(χ == η)
    @test !(χ == ζ)
    @test !(η == ζ)

    @test isequal(zero(VectorSO3Algebra), zero(VectorSO3Algebra))
    @test isequal(χ, χ)
    @test isequal(η, η)
    @test isequal(ζ, ζ)
    @test !isequal(χ, η)
    @test !isequal(χ, ζ)
    @test !isequal(η, ζ)
end
