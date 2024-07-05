#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 17:32:16
  @ license: MIT
  @ description:
 =#

const RealMatrix{Row, Column, Number} = StaticArrays.MArray{Tuple{Row, Column}, Float64, 2, Number}

const SquareMatrix{Dimension, Number} = RealMatrix{Dimension, Dimension, Number}

const Matrix2D = SquareMatrix{Dimension2, Dimension2 * Dimension2}
const Matrix3D = SquareMatrix{Dimension3, Dimension3 * Dimension3}

@inline function Matrix0(dimension::Int64 = 2)::SquareMatrix{dimension, dimension * dimension}
    return SquareMatrix{dimension, dimension * dimension}(zeros(Float64, dimension, dimension))
end

@inline function MatrixI(dimension::Int64 = 2)::SquareMatrix{dimension, dimension * dimension}
    matrix = Matrix0(dimension)
    @simd for i in 1:dimension
        @inbounds matrix[i, i] = 1.0
    end
    return matrix
end

# * ========== basic vector operations ========== * #

@inline dot(A::Matrix2D, B::Matrix2D)::Float64 = LinearAlgebra.dot(A, B)
