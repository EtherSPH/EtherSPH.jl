#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 21:03:14
  @ license: MIT
  @ description:
 =#

const kCubicSplineRadiusRatio = 2.0
const kCubicSplineSigmaList = [2.0 / 3.0, 10.0 / 7.0 / pi, 1.0 / pi]

struct CubicSpline{Dimension} <: SmoothKernel{Dimension}
    h_::Float64
    h_inv_::Float64
    radius_::Float64
    sigma_::Float64
    kernel_value_0_::Float64
end

@inline function CubicSpline{Dimension}(radius::Float64) where {Dimension}
    radius_ratio = kCubicSplineRadiusRatio
    h = radius / radius_ratio
    h_inv = 1.0 / h
    sigma = kCubicSplineSigmaList[Dimension] / h^Dimension
    kernel_value_0 = sigma
    return CubicSpline{Dimension}(h, h_inv, radius, sigma, kernel_value_0)
end

@inline @fastmath function kernelValue(r::Float64, kernel::CubicSpline{Dimension})::Float64 where {Dimension}
    q::Float64 = r * kernel.h_inv_
    if q < 1.0
        return kernel.sigma_ * (3 * q * q * (q - 2.0) + 4.0) * 0.25
    elseif q < 2.0
        return kernel.sigma_ * (2.0 - q)^3 * 0.25
    else
        return 0.0
    end
end

@inline @fastmath function kernelGradient(r::Float64, kernel::CubicSpline{Dimension})::Float64 where {Dimension}
    q::Float64 = r * kernel.h_inv_
    if q < 1.0
        return kernel.sigma_ * kernel.h_inv_ * 0.75 * q * (3 * q - 4.0)
    elseif q < 2.0
        return -kernel.sigma_ * kernel.h_inv_ * 0.75 * (2.0 - q)^2
    else
        return 0.0
    end
end
