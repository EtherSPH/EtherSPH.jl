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

const Matrix02D = Matrix0(2)
const MatrixI2D = MatrixI(2)

const Matrix03D = Matrix0(3)
const MatrixI3D = MatrixI(3)

# * ========== basic vector operations ========== * #

@inline dot(A::Matrix2D, B::Matrix2D)::Float64 = LinearAlgebra.dot(A, B)

@inline dot(A::Matrix3D, B::Matrix3D)::Float64 = LinearAlgebra.dot(A, B)

@inline function trace(A::Matrix2D)::Float64
    @inbounds return A[1] + A[4]
end

@inline function trace(A::Matrix3D)::Float64
    @inbounds return A[1] + A[5] + A[9]
end

@inline function det(A::Matrix2D)::Float64
    @inbounds return A[1] * A[4] - A[2] * A[3]
end

@inline function det(A::Matrix3D)::Float64
    @inbounds return A[1] * (A[5] * A[9] - A[6] * A[8]) - A[2] * (A[4] * A[9] - A[6] * A[7]) +
                     A[3] * (A[4] * A[8] - A[5] * A[7])
end
