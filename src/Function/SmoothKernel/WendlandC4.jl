#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 21:17:58
  @ license: MIT
  @ description:
 =#

const kWendlandC4RadiusRatio = 2.0
const kWendlandC4SigmaList = [5.0 / 8.0, 9.0 / 4.0 / pi, 495.0 / 256.0 / pi]

struct WendlandC4{Dimension} <: SmoothKernel{Dimension}
    h_::Float64
    radius_::Float64
    sigma_::Float64
    kernel_value_0_::Float64
end

@inline function WendlandC4{Dimension}(radius::Float64) where {Dimension}
    radius_ratio = kWendlandC4RadiusRatio
    h = radius / radius_ratio
    sigma = kWendlandC4SigmaList[Dimension] / h^Dimension
    kernel_value_0 = sigma
    return WendlandC4{Dimension}(h, radius, sigma, kernel_value_0)
end

@inline function kernelValue(r::Float64, kernel::WendlandC4{Dimension})::Float64 where {Dimension}
    q::Float64 = r / kernel.h_
    if q < 2.0
        return kernel.sigma_ * (2.0 - q)^6 * (35.0 * q^2 + 36.0 * q + 12.0) / 768.0
    else
        return 0.0
    end
end

@inline function kernelGradient(r::Float64, kernel::WendlandC4{Dimension})::Float64 where {Dimension}
    q::Float64 = r / kernel.h_
    if q < 2.0
        return -kernel.sigma_ / kernel.h_ * (35.0 / 96.0 * q^2 + 7.0 / 48.0 * q) * (2.0 - q)^5
    else
        return 0.0
    end
end
