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
    h_inv_::Float64
    radius_::Float64
    sigma_::Float64
    kernel_value_0_::Float64
end

@inline function WendlandC4{Dimension}(radius::Float64) where {Dimension}
    radius_ratio = kWendlandC4RadiusRatio
    h = radius / radius_ratio
    h_inv = 1.0 / h
    sigma = kWendlandC4SigmaList[Dimension] / h^Dimension
    kernel_value_0 = sigma
    return WendlandC4{Dimension}(h, h_inv, radius, sigma, kernel_value_0)
end

@inline @fastmath function kernelValue(r::Float64, kernel::WendlandC4{Dimension})::Float64 where {Dimension}
    q::Float64 = r * kernel.h_inv_
    if q < 2.0
        return kernel.sigma_ * (2.0 - q)^6 * (35 * q * q + 36 * q + 12.0) * 0.0013020833333333333
    else
        return 0.0
    end
end

@inline @fastmath function kernelGradient(r::Float64, kernel::WendlandC4{Dimension})::Float64 where {Dimension}
    q::Float64 = r * kernel.h_inv_
    if q < 2.0
        return -kernel.sigma_ * kernel.h_inv_ * (0.3645833333333333 * q * q + 0.14583333333333334 * q) * (2.0 - q)^5
    else
        return 0.0
    end
end
