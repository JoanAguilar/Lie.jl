# FIXME: Use dual quaternions instead of transformation matrices as the internal representation.
"""
	SpecialEuclidean

Concrete Lie type for the special Euclidean group in dimension 3, SE(3).

Use of the type constructor is discouraged, use the `from` method instead.
"""
struct SpecialEuclidean <: AbstractLieGroup
    T::Array{Float64, 2}
end


# TODO: Add support for dual quaternion representations.
"""
	SpecialEuclideanRepresentation

Enum with the different representations for SE(3).
"""
@enum SpecialEuclideanRepresentation begin
	matrix
	disp_rot
end


"""
	from(repr::SpecialEuclideanRepresentation, args...)

Create an instance of `SpecialEuclidean` from `args` assuming representation `repr`.
"""
function from(repr::SpecialEuclideanRepresentation, args...)
	if repr == SpecialEuclideanRepresentation.matrix
		return from_matrix(args...)
	elseif repr == SpecialEuclideanRepresentation.disp_rot
		from_disp_rot(args...)
	else
		error("Transformation from " * string(repr) * " not implemented.")
	end
end


"""
	from_matrix(T::Array{<:Real, 2})

Create an instance of `SpecialEuclidean` from the transformation matrix `T`.
"""
function from_matrix(T::Array{<:Real, 2})
	return SpecialEuclidean(copy(T))
end


"""
	from_disp_rot(d::Array{<:Real, 1}, r::SpecialOrthogonal)

Create an instance of `SpecialEuclidean` from the displacement vector `d` and the SO(3) element `r`.
"""
function from_disp_rot(d::Array{<:Real, 1}, r::SpecialOrthogonal)
	T = Array{Float64, 2}(undef, 4, 4)
	T[1:3, 1:3] = to_rotation_matrix(r)
	T[1:3, 4] = copy(d)
	T[4, 1:3] = 0
	T[4, 4] = 1
	return SpecialEuclidean(T)
end


"""
	to(repr::SpecialEuclideanRepresentation, q::SpecialEuclidean)

Convert an instance of `SpecialEuclidean` to an output that uses representation `repr`.

Note that certain representations can be multi-valued, that is, there might exist multiple valid
output values for a single `q`.
"""
function to(repr::SpecialEuclideanRepresentation, q::SpecialEuclidean)
	if repr == SpecialEuclideanRepresentation.matrix
		return to_matrix(q)
	elseif repr == SpecialEuclideanRepresentation.disp_rot
		return to_disp_rot(d)
	else
		error("Transformation to " * string(repr) * " not implemented.")
	end
end


"""
	to_matrix(q::SpecialEuclidean)

Create an `Array` containing `q` in the form of a 4×4 transformation matrix.
"""
function to_matrix(q::SpecialEuclidean)
	return copy(q.T)
end


"""
	to_disp_rot(q::SpecialEuclidean)

Create a pair of outputs (`d` and `r`) containing the displacement vector as an `Array` and the rotation as
a `SpecialOrthogonal` instance, .
"""
function to_disp_rot(q::SpecialEuclidean)
	return disp(q), rot(q)
end


"""
	disp(q::SpecialEuclidean)

Create an `Array` containing the displacement vector of `q`.
"""
function disp(q::SpecialEuclidean)
	return copy(q.T[1:3, 4])
end


"""
	rot(q::SpecialOrthogonal)

Create a `SpecialOrthogonal` instance with the rotation of `q`.
"""
function rot(q::SpecialOrthogonal)
	return from_rotation_matrix(q.T[1:3, 1:3])
end


"""
	exp(ω::Array{<:Real, 1}, T::Type{SpecialEuclidean})

Exponential map for `SpecialEuclidean`.

The input `ω` is an element of the Lie algebra se(3) as a 6-element `Array` (the first 3 elements
corresponding to displacement, and the remaining 3 to rotation). The output is an element of the Lie group
SE(3) as the corresponding instance of `SpecialEuclidean`.
"""
function Base.:exp(ω::Array{<:Real, 1}, T::Type{SpecialEuclidean})
	T = Array{Float64, 2}(undef, 4, 4)
	T[1:3, 1:3] = rodrigues(ω[1:3])
	T[1:3, 4] = copy(ω[4:6])
	T[4, 1:3] = 0
	T[4, 4] = 1
	return SpecialEuclidean(T)
end


"""
	log(q::SpecialEuclidean)

Logarithmic map for `SpecialEuclidean`.

The input is an element of the Lie group SE(3) as `q`, the output is an element of the Lie algebra se(3) as
a 6-element `Array` (the first 3 elements corresponding to displacement, and the remaining 3 to rotation).
Note that the logarithmic map is multivalued, the implementation here will return the `Array` with the
smallest norm.
"""
function Base.:log(q::SpecialEuclidean)
	ω = Array{Float64, 1}(undef, 6)
	ω[1:3] = copy(q.T[1:3, 4])
	ω[4:6] = log(rot(q))
	return ω
end


"""
	*(q::SpecialEuclidean, r::SpecialEuclidean, s...)

Lie-group multiplication for `SpecialEuclidean` instances.
"""
function Base.:*(q::SpecialEuclidean, r::SpecialEuclidean, s...)
	mul = SpecialOrthogonal(q.T * r.T)
	if length(s) == 0
		return mul
	else
		return *(mul, s...)
	end
end


"""
	*(q::SpecialEuclidean, v::Array{<:Real, 1})

Multiplication between the Lie group element `q` and the 3-element vector `v`.

The output is a 3-element vector which corresponds to `v` after applying the displacement and rotation
specified by `q`.
"""
function Base.:*(q::SpecialEuclidean, v::Array{<:Real, 1})
	return q.T[1:3, 1:3] * v + q.T[1:3, 4]
end


"""
	inv(q::SpecialEuclidean)

Lie-group multiplicative inverse of `q`.
"""
function Base.:inv(q::SpecialEuclidean)
	d, r = to_disp_rot(q)
	return from_disp_rot(inv(r) * (-d), inv(r))
end


"""
	one(q::SpecialEuclidean)

Lie-group multiplicative identity element corresponding to `q`.
"""
function Base.:one(q::SpecialEuclidean)
	return one(SpecialEuclidean)
end


"""
	one(T::Type{SpecialEuclidean})

Lie-group multiplicative identity for `SpecialEuclidean`.
"""
function Base.:one(T::Type{SpecialEuclidean})
	return SpecialEuclidean(convert(Array{Float64}, I(4)))
end
