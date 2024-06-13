#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 21:14:41
  @ license: MIT
  @ description:
 =#

kWendlandC2RadiusRatio = 2.0
kWendlandC2SigmaList = [0.0, 7.0 / 4.0 / pi, 21.0 / 16.0 / pi]

struct WendlandC2{Dimension} <: SmoothKernel{Dimension}
    h_::Float64
    radius_::Float64
    sigma_::Float64
    kernel_value_0_::Float64
end

@inline function WendlandC2{Dimension}(radius::Float64) where {Dimension}
    radius_ratio = kWendlandC2RadiusRatio
    h = radius / radius_ratio
    sigma = kWendlandC2SigmaList[Dimension] / h^Dimension
    kernel_value_0 = sigma
    return WendlandC2{Dimension}(h, radius, sigma, kernel_value_0)
end

@inline function kernelValue(r::Float64, kernel::WendlandC2{Dimension})::Float64 where {Dimension}
    q::Float64 = r / kernel.h_
    if q < 2.0
        return kernel.sigma_ * (2.0 - q)^4 * (1.0 + 2.0 * q) / 16.0
    else
        return 0.0
    end
end

@inline function kernelGradient(r::Float64, kernel::WendlandC2{Dimension})::Float64 where {Dimension}
    q::Float64 = r / kernel.h_
    if q < 2.0
        return -kernel.sigma_ / kernel.h_ * 5.0 / 8.0 * q * (2.0 - q)^3
    else
        return 0.0
    end
end
