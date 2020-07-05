module Lie

include("src/abstract_types.jl")
include("src/real_array.jl")
include("src/special_orthogonal.jl")
include("src/special_euclidean.jl")

# Abstract types
export AbstractLieGroup

# LieRealArray
export LieRealArray
export LieRealArrayRepresentation
export from, from_array
export to, to_array

# Special Orthogonal 3 - SO(3)
export SpecialOrthogonal
export SpecialOrthogonalRepresentation
export from, from_rotation_matrix, from_quaternion, from_rotation_vector, from_roll_pitch_yaw
export to, to_rotation_matrix, to_quaternion, to_rotation_vector, to_roll_pitch_yaw

# Special Euclidean 3 - SE(3)
export SpecialEuclidean
export SpecialEuclideanRepresentation
export from, from_matrix, from_disp_rot
export to, to_matrix, to_disp_rot
export disp, rot

end
