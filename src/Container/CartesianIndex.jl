#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 21:11:47
  @ license: MIT
  @ description:
 =#

const CartesianIndex2D = CartesianIndex{Dimension2}
const CartesianIndex3D = CartesianIndex{Dimension3}

@inline dimension(::CartesianIndex{Dimension}) where {Dimension} = Dimension

struct CartesianRange{Dimension}
    lower_::CartesianIndex{Dimension}
    upper_::CartesianIndex{Dimension}
end

const CartesianRange2D = CartesianRange{Dimension2}
const CartesianRange3D = CartesianRange{Dimension3}

@inline function CartesianRange(upper::CartesianIndex{Dimension})::CartesianRange{Dimension} where {Dimension}
    lower = CartesianIndex{Dimension}() # usually, the lower bound is (1, 1, ...)
    return CartesianRange(lower, upper)
end

@inline function Base.show(io::IO, range::CartesianRange{Dimension}) where {Dimension}
    print(io, "CartesianRange{", Dimension, "}(lower = ")
    show(io, range.lower_)
    print(io, ", upper = ")
    show(io, range.upper_)
    print(io, ")")
    return nothing
end

@inline function isInsideCartesianRange(
    index::CartesianIndex{Dimension},
    range::CartesianRange{Dimension},
)::Bool where {Dimension}
    @simd for i in 1:dimension(index)
        @inbounds if index[i] < range.lower_[i] || index[i] > range.upper_[i]
            return false
        end
    end
    return true
end

# const CartesianGridIndex{Dimension} = StaticArrays.MArray{Tuple{Dimension}, Int64, 1, Dimension}

# const CartesianGridIndex2D = CartesianGridIndex{Dimension2}
# const CartesianGridIndex3D = CartesianGridIndex{Dimension3}
# @inline dimension(::CartesianGridIndex{Dimension}) where {Dimension} = Dimension

# @inline function CartesianGridIndex(dimension::Int64)::CartesianGridIndex{dimension}
#     return CartesianGridIndex{dimension}(zeros(Int64, dimension))
# end

# struct CartesianGridRange{Dimension}
#     lower_::CartesianGridIndex{Dimension}
#     upper_::CartesianGridIndex{Dimension}
# end

# const CartesianGridRange2D = CartesianGridRange{Dimension2}
# const CartesianGridRange3D = CartesianGridRange{Dimension3}
# @inline dimension(::CartesianGridRange{Dimension}) where {Dimension} = Dimension

# @inline function CartesianGridRange(
#     upper::CartesianGridIndex{Dimension},
# )::CartesianGridRange{Dimension} where {Dimension}
#     lower = CartesianGridIndex{Dimension}(ones(Int64, dimension(upper))) # usually, the lower bound is (1, 1, ...)
#     return CartesianGridRange(lower, upper)
# end

# @inline function Base.show(io::IO, grid_range::CartesianGridRange{Dimension}) where {Dimension}
#     print(io, "CartesianGridRange{", Dimension, "}(lower = ")
#     show(io, grid_range.lower_)
#     print(io, ", upper = ")
#     show(io, grid_range.upper_)
#     print(io, ")")
#     return nothing
# end

# @inline function Base.CartesianIndex(index::CartesianGridIndex{Dimension})::CartesianIndex{Dimension} where {Dimension}
#     return CartesianIndex(index...)
# end

# @inline function isInsideCartesianGridRange(
#     index::CartesianGridIndex{Dimension},
#     grid_range::CartesianGridRange{Dimension},
# )::Bool where {Dimension}
#     @simd for i in 1:dimension(index)
#         @inbounds if index[i] < grid_range.lower_[i] || index[i] > grid_range.upper_[i]
#             return false
#         end
#     end
#     return true
# end
