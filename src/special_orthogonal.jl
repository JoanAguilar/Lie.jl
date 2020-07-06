using LinearAlgebra: eigvecs, I, norm, tr


# FIXME: Use quaternions instead of rotation matrices as the internal representation.
"""
	SpecialOrthogonal

Concrete Lie type for the special orthogonal group in dimension 3.

Use of the type constructor is discouraged, use the `from` method instead.
"""
struct SpecialOrthogonal <: AbstractLieGroup
    R::Array{Float64, 2}
end


"""
	SpecialOrthogonalRepresentation

Enum with the different representations for SO(3).
"""
@enum SpecialOrthogonalRepresentation begin
	rotation_matrix
	quaternion
	rotation_vector
	roll_pitch_yaw
end


"""
	from(repr::SpecialOrthogonalRepresentation, array::Array{<:Real})

Create an instance of `SpecialOrthogonal` from `array` that's using the representation `repr`.
"""
function from(repr::SpecialOrthogonalRepresentation, array::Array{<:Real})
	if repr == SpecialOrthogonalRepresentation.rotation_matrix
 		return from_rotation_matrix(array)
	elseif repr == SpecialOrthogonalRepresentation.quaternion
		return from_quaternion(array)
	elseif repr == SpecialOrthogonalRepresentation.rotation_vector
		return from_rotation_vector(array)
	elseif repr == SpecialOrthogonalRepresentation.roll_pitch_yaw
		return from_roll_pitch_yaw(array)
	else
		error("Transformation from " * string(repr) * " not implemented.")
	end
end


"""
	from_rotation(R::Array{<:Real, 2})

Create an instance of `SpecialOrthogonal from the rotation matrix `R`.
"""
function from_rotation_matrix(R::Array{<:Real, 2})
	return SpecialOrthogonal(copy(R))
end


"""
	from_quaternion(q::Array{<:Real, 1})

Create an instance of `SpecialOrthogonal` from the quaternion `q`.
"""
function from_quaternion(q::Array{<:Real, 1})
	return SpecialOrthogonal(
		(q[1]^2 - q[2:4]' * q[2:4]) * I(3) +
		2 * q[2:4] * q[2:4]' +
		2 * q[1] * skew(q[2:4])
	)
end


"""
	from_rotation_vector(r::Array{<:Real, 1})

Create an instance of `SpecialOrthogonal` from the rotation vector `r`.
"""
function from_rotation_vector(r::Array{<:Real, 1})
	return SpecialOrthogonal(rodrigues(r))
end


"""
	from_roll_pitch_yaw(rpy::Array{<:Real, 1})

Create an instance of `SpecialOrthogonal` from the roll-pitch-yaw angles contained in the array `rpy`.

To clarify: `rpy[1]` is roll, `rpy[2]` is pitch, and `rpy[3]` is yaw. The angles must be in radians.
"""
function from_roll_pitch_yaw(rpy::Array{<:Real, 1})
	r = rpy[1]
	p = rpy[2]
	y = rpy[3]
	sr = sin(r)
	sp = sin(p)
	sy = sin(y)
	cr = cos(r)
	cp = cos(p)
	cy = cos(y)
	R = [cp * cy    sr * sp * cy - cr * sy    sr * sy + cr * sp * cy;
	     cp * sy    sr * sp * sy + cr * cy    -sr * sy + cr * sp * sy;
	     -sp        sr * cp                   cr * cp]
	return SpecialOrthogonal(R)
end


"""
	rodrigues(r::Array{Real, 1})

Rodrigues' rotation formula for the rotation vector `r`.
"""
function rodrigues(r::Array{Real, 1})
	θ = norm(r)
	if θ ~= 0
		return one(SpecialOrthogonal)
	else
		ax = r / θ
		K = skew(ax)
		return I(3) + sin(θ) * K + (1 - cos(θ)) * K^2
	end
end


"""
	skew(r::Array{Real, 1})

Skew-symmetric operator for the vector `r`.
"""
function skew(r::Array{Real, 1})
	return [ 0    -r[3] r[2]
		 r[3] 0     -r[1]
		-r[2] r[1]  0]
end


"""
	to(repr::SpecialOrthogonalRepresentation, q::SpecialOrthogonal)

Convert an instance of `SpecialOrthogonal` to an `Array` that uses representation `repr`.

Note that certain representations can be multi-valued, that is, there might exist multiple valid
output values for a single `q`.
"""
function to(repr::SpecialOrthogonalRepresentation, q::SpecialOrthogonal)
	if repr == SpecialOrthogonalRepresentation.rotation_matrix
 		return to_rotation_matrix(array)
	elseif repr == SpecialOrthogonalRepresentation.quaternion
		return to_quaternion(array)
	elseif repr == SpecialOrthogonalRepresentation.rotation_vector
		return to_rotation_vector(array)
	elseif repr == SpecialOrthogonalRepresentation.roll_pitch_yaw
		return to_roll_pitch_yaw(array)
	else
		error("Transformation to " * string(repr) * " not implemented.")
	end
end


"""
	to_rotation_matrix(q::SpecialOrthogonal)

Create an `Array` containing `q` in the form of a rotation matrix.
"""
function to_rotation_matrix(q::SpecialOrthogonal)
	return copy(q.R)
end


"""
	to_quaternion(q::SpecialOrthogonal)

Create an `Array` containing `q` in the form of a quaternion.

Note that there are two possible quaternions that can represent `q`. The implementation here will
choose one of the two.
"""
function to_quaternion(q::SpecialOrthogonal)
	# Code adapted from:
	# Maths - Conversion Matrix to Quaternion
	# by Martin John Baker
	# https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

	R = q.R
	T = tr(R)

	if T > 0
		S = 2 * √(T + 1)
		qw = S / 4
		qx = (R[3, 2] - R[2, 3]) / S
		qy = (R[1, 3] - R[3, 1]) / S
		qz = (R[2, 1] - R[1, 2]) / S
	elseif R[1, 1] > R[2, 2] && R[1, 1] > R[3, 3] 
		S = 2 * √(1 + R[1, 1] - R[2, 2] - R[3, 3])
		qw = (R[3, 2] - R[2, 3]) / S
		qx = S / 4
		qy = (R[1, 2] + R[2, 1]) / S
		qz = (R[1, 3] + R[3, 1]) / S
	elseif R[2, 2] > R[3, 3]
		S = 2 * √(1 + R[2, 2] - R[1, 1] - R[3, 3])
		qw = (R[1, 3] - R[3, 1]) / S
		qx = (R[1, 2] + R[2, 1]) / S
		qy = S / 4
		qz = (R[2, 3] + R[3, 2]) / S
	else
		S = 2 * √(1 + R[3, 3] - R[1, 1] - R[2, 2])
		qw = (R[2, 1] - R[1, 2]) / S
		qx = (R[1, 3] + R[3, 0]) / S
		qy = (R[2, 3] + R[3, 2]) / S
		qz = S / 4
	end

	return [qw, qx, qy, qz]
end


"""
	to_rotation_vector(q::SpecialOrthogonal)

Create an `Array` containing `q` in the form of a rotation vector.

Note that there are multiple rotation vectors that can represent `q`. The implementation here
will return the one with the shortest norm.
"""
function to_rotation_vector(q::SpecialOrthogonal)
	θ = acos((tr(q.R) - 1) / 2)
	if θ ~= 0
		return zeros(3)
	elseif θ ~= π
		# Only use Eigendecomposition if the angle value is close to π.
		@warn "Rotation angle value is close to π, using Eigendecomposition."
		ax = eigvecs(q.R)[:, 1]
		return θ * ax
	else
		R = q.R
		dsθ = 2 * sin(θ)
		return [(R[3, 2] - R[2, 3]) / dsθ,
			(R[1, 3] - R[3, 1]) / dsθ,
			(R[2, 1] - R[1, 2]) / dsθ]
	end
end


"""
	to_roll_pitch_yaw(q::SpecialOrthogonal)

Create an `Array` containing `q` as a sequence of roll-pitch-yaw angles (in radians).

Note that there are multiple sets of roll-pitch-yaw angle values that can represent `q`. The
implementation here will choose one of such sets.
"""
function to_roll_pitch_yaw(q::SpecialOrthogonal)
	R = q.R
	r = atan(R[3, 2] / R[3, 3])
	p = atan(-R[3, 1] / √(1 - R[3, 1]^2))
	y = atan(R[2, 1] / R[1, 1])
	return [r, p, y]
end


"""
	exp(T::Type{SpecialOrthogonal}, ω::Array{<:Real, 1})

Exponential map for `SpecialOrthogonal`.

The input `ω` is an element of the Lie algebra so(3) as a 3-element `Array`, the output is an element of
the Lie group SO(3) as the corresponding instance of `SpecialOrthogonal`
"""
function Base.:exp(T::Type{SpecialOrthogonal}, ω::Array{<:Real, 1})
	return from_rotation_vector(ω)
end


"""
	log(q::SpecialOrthogonal)

Logarithmic map for `SpecialOrthogonal`.

The input is an element of the Lie group SO(3) as `q`, the output is an element of the Lie algebra so(3) as
a 3-element `Array`. Note that the logarithmic map is multivalued, the implementation here will return the
`Array` with the smallest norm.
"""
function Base.:log(q::SpecialOrthogonal)
	return to_rotation_vector(q)
end


"""
	*(q::SpecialOrthogonal, r::SpecialOrthogonal, s...)

Lie-group multiplication for `SpecialOrthogonal` instances.
"""
function Base.:*(q::SpecialOrthogonal, r::SpecialOrthogonal, s...)
	mul = SpecialOrthogonal(q.R * r.R)
	if length(s) == 0
		return mul
	else
		return *(mul, s...)
	end
end


"""
	*(q::SpecialOrthogonal, v::Array{<:Real, 1})

Multiplication between the Lie group element `q` and the 3-element vector `v`.

The output is a 3-element vector which corresponds to `v` after applying the rotation specified by `q`.

Equivalent to `to_rotation_matrix(q) * v`.
"""
function Base.:*(q::SpecialOrthogonal, v::Array{<:Real, 1})
	return q.R * v
end


"""
	inv(q::SpecialOrthogonal)

Lie-group multiplicative inverse of `q`.
"""
function Base.:inv(q::SpecialOrthogonal)
	return SpecialOrthogonal(q.R')
end


"""
	one(q::SpecialOrthogonal)

Lie-group multiplicative identity element corresponding to `q`.
"""
function Base.:one(q::SpecialOrthogonal)
	return one(SpecialOrthogonal)
end


"""
	one(T::Type{SpecialOrthogonal})

Lie-group multiplicative identity for `SpecialOrthogonal`.
"""
function Base.:one(T::Type{SpecialOrthogonal})
	return SpecialOrthogonal(convert(Array{Float64}, I(3)))
end
