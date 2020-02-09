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
    mat = convert(Array, R)
    θ = acos((mat[1, 1] + mat[2, 2] + mat[3, 3] - 1) / 2)
    if θ > 0
        ax_den = 2 * sin(θ)
        ax_num_vec = [mat[3, 2] - mat[2, 3],
                      mat[1, 3] - mat[3, 1],
                      mat[2, 1] - mat[1, 2]]
        return T(θ * ax_num_vec / ax_den, checks=false)
    else
        return zero(T)
    end
end
