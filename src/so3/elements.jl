"""
    one(x)

Identity element.
"""
function Base.:one(x::T) where T<:RotationMatrix
    return T([1 0 0; 0 1 0; 0 0 1], checks=false)
end


"""
    one(T::Type)

Identity element.
"""
function Base.:one(T::Type{RotationMatrix})
    return T([1 0 0; 0 1 0; 0 0 1], checks=false)
end


"""
    zero(ω)

Zero element.
"""
function Base.:zero(ω::T) where T<:VectorSO3Algebra
    return T([0; 0; 0], checks=false)
end


"""
    zero(T::Type)

Zero element.
"""
function Base.:zero(T::Type{VectorSO3Algebra})
    return T([0; 0; 0], checks=false)
end
