#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 21:22:17
  @ license: MIT
  @ description:
 =#

@inline function libKernelAverageDensityFilterInteraction!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_value::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.sum_kernel_weight_ += kernel_value * q.mass_ / q.rho_
    p.sum_kernel_weighted_value_ += kernel_value * q.mass_
    return nothing
end

@inline function libKernelAverageDensityFilterSelfaction!(
    p::ParticleType;
    kernel_value::Float64,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.sum_kernel_weight_ += kernel_value * p.mass_ / p.rho_
    p.sum_kernel_weighted_value_ += kernel_value * p.mass_
    p.rho_ = p.sum_kernel_weighted_value_ / p.sum_kernel_weight_
    p.sum_kernel_weight_ = 0.0
    p.sum_kernel_weighted_value_ = 0.0
    return nothing
end
