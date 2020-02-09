module Lie

include("src/abstract_types.jl")
include("src/so3/so3.jl")

# Abstract Lie types
export AbstractLieGroup
export AbstractLieAlgebra

# SO(3) types
export AbstractSO3Group
export AbstractSO3Algebra
export RotationMatrix
export VectorSO3Algebra, angle, axis

end
