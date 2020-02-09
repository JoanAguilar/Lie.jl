"""
    convert(T::Type, R)

Convert rotation matrix `R` to an array of type `T`.
"""
function Base.:convert(T::Type{Array}, R::RotationMatrix)
    return convert(T, R.R)
end


"""
    angle(R)

Rotation angle for the rotation matrix `R`.
"""
function angle(R::RotationMatrix)
    return angle(log(R, VectorSO3Algebra))
end

"""
    axis(ω)

Rotation unit axis for the rotation matrix `R`. Returns an array full of `NaN` if the rotation
angle is ``0`` or a multiple of ``π``.
"""
function axis(R::RotationMatrix)
    return axis(log(R, VectorSO3Algebra))
end


"""
    isapprox(R, S; args...)

Inexact equality comparison for rotation matrices. The comparison is performed by computing the
angle `θ` between the two matrices and then calling `isapprox(0, θ, args...)`. See `Base`'s
`isapprox` for more information about `args`.
"""
function Base.:isapprox(R::RotationMatrix, S::RotationMatrix; args...)
    θ = angle(log(inv(R) * S, VectorSO3Algebra))
    return isapprox(0, θ; args...)
end


"""
    isequal(R, S)

Equality comparison for rotation matrices. The comparison is performed by applying `isequal` to
the pair of arrays representing the rotation matrices. See `Base`'s `isequal` for information on
how `isequal` differs from `==`.
"""
function Base.:isequal(R::RotationMatrix, S::RotationMatrix)
    return isequal(convert(Array, R), convert(Array, S))
end


"""
    ==(R, S)

Generic equality comparison for rotation matrices. The comparison is performed by applying
`==` to the pair of arrays representing the rotation matrices. See `Base`'s `==` for information
on how `==` differs from `isequal`.
"""
function Base.:(==)(R::RotationMatrix, S::RotationMatrix)
    return ==(convert(Array, R), convert(Array, S))
end


"""
    convert(T::Type, ω)

Convert ``so3`` element `ω` to an array of type `T`.
"""
function Base.:convert(T::Type{Array}, ω::VectorSO3Algebra)
    return convert(T, ω.ω)
end

"""
    angle(ω)

Rotation angle for the ``so3`` element `ω`.
"""
function angle(ω::VectorSO3Algebra)
    return norm(convert(Array, ω))
end

"""
    axis(ω)

Rotation unit axis for the ``so3`` element `ω`. Returns an array full of `NaN` if the rotation
angle is ``0`` or a multiple of ``π``.
"""
function axis(ω::VectorSO3Algebra)
    θ = angle(ω)
    if θ % π == 0
        return [NaN; NaN; NaN]
    else
        return convert(Array, ω) / θ
    end
end


"""
    isapprox(ω, χ; args...)

Inexact equality comparison for ``so3`` elements. The comparison is performed by computing the
angle `θ` between the two ``so3`` elements and then calling `isapprox(0, θ; args..)`. See
`Base`'s `isapprox` for more information about `args`.
"""
function Base.:isapprox(ω::VectorSO3Algebra, χ::VectorSO3Algebra; args...)
    return isapprox(0, angle(ω - χ); args...)
end


"""
    isequal(ω, χ)

Equality comparison for ``so3`` elements. The comparison is performed by applying `isequal` to
the pair of vectors representing the ``so3`` elements. See `Base`'s `isequal` for information on
how `isequal` differs from `==`.
"""
function Base.:isequal(ω::VectorSO3Algebra, χ::VectorSO3Algebra)
    return isequal(convert(Array, ω), convert(Array, χ))
end


"""
    ==(ω, χ)

Generic equality comparison for ``so3`` elements. The comparison is performed by applying `==` to
the pair of vectors representing the ``so3`` elements. See `Base`'s `==` for information on how
`==` differs from `isequal`.
"""
function Base.:(==)(ω::VectorSO3Algebra, χ::VectorSO3Algebra)
    return ==(convert(Array, ω), convert(Array, χ))
end
