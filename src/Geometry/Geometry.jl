#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 22:15:11
  @ license: MIT
  @ description:
 =#

abstract type Shape end

const Point{Dimension} = RealVector{Dimension}
const Point2D = Point{Dimension2}
const Point3D = Point{Dimension3}

struct Box{Dimension} <: Shape
    lower_::Point{Dimension}
    upper_::Point{Dimension}
    range_::RealVector{Dimension}
end
@inline dimension(::Box{Dimension}) where {Dimension} = Dimension

const Box2D = Box{Dimension2}
const Box3D = Box{Dimension3}

@inline function Box(lower::Point{Dimension}, upper::Point{Dimension})::Box{Dimension} where {Dimension}
    return Box(lower, upper, upper - lower)
end

@inline function isInsideShape(point::Point{Dimension}, box::Box{Dimension}) where {Dimension}
    @simd for i in 1:dimension(point)
        @inbounds if point[i] < box.lower_[i] || point[i] > box.upper_[i]
            return false
        end
    end
    return true
end

include("Rectangle.jl")
include("Ring.jl")
