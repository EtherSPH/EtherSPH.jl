#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 21:05:23
  @ license: MIT
  @ description:
 =#

const kGaussianRadiusRatio = 3.0
const kGaussianSigmaList = [1.0 / sqrt(pi), 1.0 / pi, 1.0 / sqrt(pi^3)]

struct Gaussian{Dimension} <: SmoothKernel{Dimension}
    h_::Float64
    h_inv_::Float64
    radius_::Float64
    sigma_::Float64
    kernel_value_0_::Float64
end

@inline function Gaussian{Dimension}(radius::Float64) where {Dimension}
    radius_ratio = kGaussianRadiusRatio
    h = radius / radius_ratio
    h_inv = 1.0 / h
    sigma = kGaussianSigmaList[Dimension] / h^Dimension
    kernel_value_0 = sigma
    return Gaussian{Dimension}(h, h_inv, radius, sigma, kernel_value_0)
end

@inline @fastmath function kernelValue(r::Float64, kernel::Gaussian{Dimension})::Float64 where {Dimension}
    q::Float64 = r * kernel.h_inv_
    if q < 3.0
        return kernel.sigma_ * exp(-q * q)
    else
        return 0.0
    end
end

@inline @fastmath function kernelGradient(r::Float64, kernel::Gaussian{Dimension})::Float64 where {Dimension}
    q::Float64 = r * kernel.h_inv_
    if q < 3.0
        return -2.0 * kernel.sigma_ * kernel.h_inv_ * q * exp(-q * q)
    else
        return 0.0
    end
end
