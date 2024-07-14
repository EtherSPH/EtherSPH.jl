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
