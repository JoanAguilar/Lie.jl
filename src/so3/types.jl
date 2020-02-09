using LinearAlgebra: det, I


"""
    AbstractSO3Group

Abstract type for ``SO(3)``.

Every implementation of ``SO(3)`` must be a subtype of this type.
"""
abstract type AbstractSO3Group <: AbstractLieGroup end


"""
    RotationMatrix{T}

Concrete type that stores an ``SO(3)`` Group element as a 3×3 rotation matrix.

The elements of the matrix are stored as scalars of type `T`.
"""
struct RotationMatrix{T<:Real} <: AbstractSO3Group
    R::Array{T, 2}


    """
        RotationMatrix(x, [checks=true])

    Construct a `RotationMatrix` instance from array `x`, if `checks` is set to `true`, correctness
    checks are performed during construction.
    """
    function RotationMatrix(R::Array{T, 2}; checks::Bool=true) where T<:Real
        if checks
            if size(R) ≠ (3, 3)
                throw(ArgumentError("Expected size (3, 3), got $(size(R)). `R` is not a " *
				    "rotation matrix."))
            elseif R * R' ≉ I
                throw(ArgumentError("Transpose is not an inverse. `R` is not a rotation  matrix."))
            elseif det(R) ≉ 1
                throw(ArgumentError("Expected determinant 1, got $(det(R)). `x` is not a " *
                                    "rotation matrix."))
            end
        end
	new{T}(R)
    end
end


"""
    AbstractSO3Algebra

Abstract type for ``so3``.

Every implementation of ``so3`` must be a subtype of this type.
"""
abstract type AbstractSO3Algebra <: AbstractLieAlgebra end


"""
    VectorSO3Algebra{T}

Concrete type that stores an ``so3`` Algebra element as 3-element vector.
"""
struct VectorSO3Algebra{T<:Real} <: AbstractLieAlgebra
    ω::Array{T, 1}

    """
        VectorSO3Algebra(ω, [checks=true])

    Construct a `VectorSO3Algebra` instance from vector `ω`, if `checks` is set to `true`,
    correctness checks are performed during construction.
    """
    function VectorSO3Algebra(ω::Array{T, 1}; checks::Bool=true) where T<:Real
        if checks
            if size(ω) ≠ (3,)
                throw(ArgumentError("Expected size (3,) got $(size(ω)). `ω` can't be used as a " *
                                    "representation of an ``so3`` algebra element."))
            end
        end
        new{T}(ω)
    end
end
