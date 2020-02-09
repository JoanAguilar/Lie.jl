"""
    *(R, S, T...)

Rotation matrix multiplication.
"""
function Base.:*(R::RotationMatrix, S::RotationMatrix, T...)
    if length(T) == 0
        return RotationMatrix(convert(Array, R) * convert(Array, S), checks=false)
    else
        return *(RotationMatrix(convert(Array, R) * convert(Array, S), checks=false), T...)
    end
end


"""
   inv(R)

Rotation matrix inverse, which corresponds to the matrix transpose.
"""
function Base.:inv(R::RotationMatrix)
    return RotationMatrix(convert(Array, convert(Array, R)'), checks=false)
end


"""
    +(ω, χ, η...)

``so3`` element summation.
"""
function Base.:+(ω::VectorSO3Algebra, χ::VectorSO3Algebra, η...)
    if length(η) == 0
        return VectorSO3Algebra(convert(Array, ω) + convert(Array, χ), checks=false)
    else
        return +(VectorSO3Algebra(convert(Array, ω) + convert(Array, χ), checks=false), η...)
    end
end


"""
    -(ω, χ)

``so3`` element subtraction.
"""
function Base.:-(ω::VectorSO3Algebra, χ::VectorSO3Algebra)
    return +(ω, -χ)
end


"""
    -(ω)

``so3`` element negation.
"""
function Base.:-(ω::VectorSO3Algebra)
    return VectorSO3Algebra(-convert(Array, ω), checks=false)
end


"""
    *(x, ω, η...)

Multiplication between a scalar `x` and an ``so3`` element `ω`.
"""
function Base.:*(x::T, ω::VectorSO3Algebra, η...) where T<:Real
    if length(η) == 0
        return VectorSO3Algebra(x * convert(Array, ω), checks=false)
    else
        return *(VectorSO3Algebra(x * convert(Array, ω), checks=false), η...)
    end
end


"""
    *(ω, x, η...)

Multiplication between a scalar `x` and an ``so3`` element `ω`.
"""
function Base.:*(ω::VectorSO3Algebra, x::T, η...) where T<:Real
    return *(x, ω, a...)
end


"""
    *(R, ω)

Multiplication between a rotation matrix `R` and an ``so3`` element `ω`.
"""
function Base.:*(R::RotationMatrix, ω::VectorSO3Algebra)
    return VectorSO3Algebra(convert(Array, R) * convert(Array, ω))
end
