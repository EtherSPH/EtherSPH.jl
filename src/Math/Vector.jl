#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 16:49:27
  @ license: MIT
  @ description:
 =#

# * ========== basic vector definitions ========== * #

const RealVector{Dimension} = StaticArrays.MArray{Tuple{Dimension}, Float64, 1, Dimension}

@inline dimension(::RealVector{Dimension}) where {Dimension} = Dimension

const Vector2D = RealVector{Dimension2}
const Vector3D = RealVector{Dimension3}

@inline function Vector0(dimension::Int64 = 2)::RealVector{dimension}
    return RealVector{dimension}(zeros(Float64, dimension))
end

@inline function VectorX(dimension::Int64 = 2)::RealVector{dimension}
    vector = Vector0(dimension)
    @inbounds vector[1] = 1.0
    return vector
end

@inline function VectorY(dimension::Int64 = 2)::RealVector{dimension}
    vector = Vector0(dimension)
    @inbounds vector[2] = 1.0
    return vector
end

@inline function VectorZ(dimension::Int64 = 3)::RealVector{dimension}
    vector = Vector0(dimension)
    @inbounds vector[3] = 1.0
    return vector
end

const Vector02D = Vector0(2)
const VectorX2D = VectorX(2)
const VectorY2D = VectorY(2)

const Vector03D = Vector0(3)
const VectorX3D = VectorX(3)
const VectorY3D = VectorY(3)
const VectorZ3D = VectorZ(3)

# * ========== basic vector operations ========== * #

@inline function dot(x::Vector2D, y::Vector2D)::Float64
    @inbounds return x[1] * y[1] + x[2] * y[2]
end

@inline function dot(x::Vector3D, y::Vector3D)::Float64
    @inbounds return x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
end

@inline function minusplus(x::Vector2D, y::Vector2D, z::Vector2D)::Vector2D
    @inbounds return Vector2D(x[1] - y[1] + z[1], x[2] - y[2] + z[2])
end

@inline function minusplus(x::Vector3D, y::Vector3D, z::Vector3D)::Vector3D
    @inbounds return Vector3D(x[1] - y[1] + z[1], x[2] - y[2] + z[2], x[3] - y[3] + z[3])
end

@inline function norm(x::RealVector{Dimension})::Float64 where {Dimension}
    @inbounds return sqrt(dot(x, x))
end

@inline function normalize(x::RealVector{Dimension})::RealVector{Dimension} where {Dimension}
    return x / norm(x)
end

@inline function cross(x::Vector2D, y::Vector2D)::Float64
    @inbounds return x[1] * y[2] - x[2] * y[1]
end

@inline function cross(x::Vector3D, y::Vector3D)::Vector3D
    @inbounds return Vector3D(x[2] * y[3] - x[3] * y[2], x[3] * y[1] - x[1] * y[3], x[1] * y[2] - x[2] * y[1])
end
