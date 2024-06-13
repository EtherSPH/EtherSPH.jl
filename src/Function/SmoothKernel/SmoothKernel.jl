#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 20:02:12
  @ license: MIT
  @ description:
 =#

abstract type SmoothKernel{Dimension} end

# ! @code_warntype -> not stable
@inline function asCallable(smooth_kernel::SmoothKernel{Dimension})::Tuple{Function, Function} where {Dimension}
    @inline function value(r::Float64)::Float64
        return kernelValue(r, smooth_kernel)
    end
    @inline function gradient(r::Float64)::Float64
        return kernelGradient(r, smooth_kernel)
    end
    return value, gradient
end

include("CubicSpline.jl")
include("Gaussian.jl")
include("WendlandC2.jl")
include("WendlandC4.jl")
