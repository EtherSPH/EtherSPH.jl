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

macro Vec(args)
    return :(StaticArrays.@MArray($args))
end

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

# * ========== basic vector operations ========== * #

@inline function dot(x::RealVector{Dimension}, y::RealVector{Dimension})::Float64 where {Dimension}
    return sum(x .* y)
end

@inline function norm(x::RealVector{Dimension})::Float64 where {Dimension}
    return sqrt(dot(x, x))
end

@inline function normalize(x::RealVector{Dimension})::RealVector{Dimension} where {Dimension}
    return x / norm(x)
end

@inline function cross(x::Vector3D, y::Vector3D)::Vector3D
    @inbounds return Vector3D(x[2] * y[3] - x[3] * y[2], x[3] * y[1] - x[1] * y[3], x[1] * y[2] - x[2] * y[1])
end
