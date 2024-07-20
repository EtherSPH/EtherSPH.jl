#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/13 22:16:32
  @ license: MIT
  @ description:
 =#

@inline function dyad(x::Vector2D, y::Vector2D)::Matrix2D
    @inbounds return Matrix2D(x[1] * y[1], x[2] * y[1], x[1] * y[2], x[2] * y[2])
end

@inline function dyad(x::Vector3D, y::Vector3D)::Matrix3D
    @inbounds return Matrix3D(
        x[1] * y[1],
        x[2] * y[1],
        x[3] * y[1],
        x[1] * y[2],
        x[2] * y[2],
        x[3] * y[2],
        x[1] * y[3],
        x[2] * y[3],
        x[3] * y[3],
    )
end

@inline function dot(A::Matrix2D, x::Vector2D)::Vector2D
    @inbounds return Vector2D(A[1] * x[1] + A[3] * x[2], A[2] * x[1] + A[4] * x[2])
end

@inline function dot(x::Vector2D, A::Matrix2D)::Vector2D
    @inbounds return Vector2D(A[1] * x[1] + A[2] * x[2], A[3] * x[1] + A[4] * x[2])
end

@inline function dot(A::Matrix3D, x::Vector3D)::Vector3D
    @inbounds return Vector3D(
        A[1] * x[1] + A[4] * x[2] + A[7] * x[3],
        A[2] * x[1] + A[5] * x[2] + A[8] * x[3],
        A[3] * x[1] + A[6] * x[2] + A[9] * x[3],
    )
end

@inline function dot(x::Vector3D, A::Matrix3D)::Vector3D
    @inbounds return Vector3D(
        A[1] * x[1] + A[2] * x[2] + A[3] * x[3],
        A[4] * x[1] + A[5] * x[2] + A[6] * x[3],
        A[7] * x[1] + A[8] * x[2] + A[9] * x[3],
    )
end

# * ========== funny operator redefinition ========== * #

@inline ·(x, y) = dot(x, y) # \cdotp
@inline ⋯(x, y) = ddot(x, y) # \cdots
@inline ×(x, y) = cross(x, y) # times
@inline ⊗(x, y) = dyad(x, y) # otimes
@inline tr(A) = trace(A)
