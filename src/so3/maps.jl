using LinearAlgebra: I, norm


"""
    exp(ω, T::Type)

Exponential map for the ``so3`` element `ω`.
"""
function Base.:exp(ω::VectorSO3Algebra, T::Type{RotationMatrix})
    θ = angle(ω)
    if θ > 0
        ax = axis(ω)
	k = [0 -ax[3] ax[2];
             ax[3] 0 -ax[1];
             -ax[2] ax[1] 0]
        R = I + sin(θ) * k + (1 - cos(θ)) * k^2
        return T(R, checks=false)
    else
        return one(T)
    end
end


"""
    log(R, T::Type)

Logarithmic map for the ``SO(3)`` element `R`.
"""
function Base.:log(R::RotationMatrix, T::Type{VectorSO3Algebra})
    θ = angle(R)
    if θ == 0
        return zero(T)
    elseif θ == π
        throw(ArgumentError("Rotation angle of input rotation matrix is π. Logarithmic map is " *
                            "not unique."))
    else
        return T(θ * axis(R), checks=false)
    end
end
