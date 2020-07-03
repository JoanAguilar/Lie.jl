# TODO: Parametrize array element type.
"""
	LieRealArray

Concrete Lie type for real-valued arrays.

Note that the internal representation uses `float64` as the array element type.

Use of the type constructor is discouraged, use the `from` method instead.
"""
struct LieRealArray <: AbstractLieGroup
    a::Array{Float64}
end


"""
	LieRealArrayRepresentation

Enum with the different representations for a real array.

As this is a trivial case, only one representation is available (`array`), and is provided in order to
match the interface of other Lie types.
"""
@enum LieRealArrayRepresentation begin
	array
end


"""
	from(array::Array{<:Real}, repr::LieRealArrayRepresentation)

Create an instance of `LieRealArray` from `array` that's using the representation `repr`.

As this is a trivial case, only one representation is available (`array`), and the method is only provided
in order to match the interface of other Lie types.
"""
function from(array::Array{<:Real}, repr::LieRealArrayRepresentation)
	if repr == LieRealArrayRepresentation.array
		return from_array(array)
	else
		error("The only representation available for `LieRealArray` is \"array\".")
	end
end


"""
	from_array(array::Array{<:Real})

Create an instance of `LieRealArray` from `array`.

The `Array` contained in the newly created instance is a copy of `array`.
"""
function from_array(array::Array{<:Real})
	return LieRealArray(copy(array))
end


"""
	to(q::LieRealArray, repr::LieRealArrayRepresentation)

Convert an instance of `LieRealArray` to an `Array` that uses representation `repr`.

As this is a trivial case, only one representation is available (`array`), and the method is only provided
in order to match the interface of other Lie types.
"""
function to(q::LieRealArray, repr::LieRealArrayRepresentation)
	if repr == LieRealArrayRepresentation.array
		return to_array(q)
	else
		error("The only representation available for `LieRealArray` is \"array\".")
	end
end

 
"""
	to_array(q::LieRealArray)

Create an `Array` from `q`.

The output is a copy of the `Array` contained in `q`.
"""
function to_array(q::LieRealArray)
	return copy(q.array)
end


"""
	exp(ω::Array{<:Real}, T::Type{LieRealArra})

Exponential map for `LieRealArray`.

The output is an instance of `LieRealArray` containing a copy of `ω`.
"""
function Base.:exp(ω::Array{<:Real}, T::Type{LieRealArray})
	return from_array(ω)
end


"""
	log(q::LieRealArray)

Logarithmic map for `LieRealArray`.

The output is a copy of the array contained by `q`.
"""
function Base.:log(q::LieRealArray)
	return to_array(q)
end


"""
	*(q::LieRealArray, r::LieRealArray, s...)

Lie-group multiplication for `LieRealArray` instances.

Note that Lie-group multiplication corresponds to matrix addition.
"""
function Base.:*(q::LieRealArray, r::LieRealArray, s...)
	mul = LieRealArray(q.array + r.array)
	if length(s) == 0
		return mul
	else
		return *(mul, s...)
	end
end


"""
	inv(q::LieRealArray)

Lie-group multiplicative inverse of `q`.

Note that the Lie-group multiplication inverse corresponds to the matrix additive inverse (`-q`).
"""
function Base.:inv(q::LieRealArray)
	return LieRealArray(-q.array)
end


"""
	*(q::LieRealArray, ω::Array{<:Real})

Multiplication between the Lie group element `q` and the Lie algebra element `ω`.

Note that this operation for `LieRealArray` has no effect on `ω` and this method is only provided to match
the interface of other Lie types.
"""
function Base.:*(q::LieRealArray, ω::Array{<:Real})
	if size(to_array(q)) != size(ω)
		error("Both arrays must be the same size, got " * string(size(to_array(q))) * " and " *
		      string(size(ω)) * ".")
	end
	return copy(ω)
end


"""
	one(q::LieRealArray)

Lie-group multiplicative identity element corresponding to `q`.
"""
function Base.:one(q::LieRealArray)
	return LieRealArray(zeros(size(q.array)...))
end


# FIXME: Avoid the `dims` argument.
"""
	one(T::Type{LieRealArray}, dims::Tuple)

Lie-group multiplicative identity for `LieRealArray`s with dimensions `dims`.
"""
function Base.:one(T::Type{LieRealArray}, dims::Tuple)
	return LieRealArray(zeros(dims...))
end
