#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 17:58:48
  @ license: MIT
  @ description:
 =#

@inline function libTraditionalPressureForce!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
    kernel_value::Float64 = 0.0,
    reference_kernel_value::Float64 = 1.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p_rho_2 = p.p_ / (p.rho_ * p.rho_) + q.p_ / (q.rho_ * q.rho_)
    p_rho_2 += abs(p_rho_2) * kernel_value / reference_kernel_value * 0.01
    p.dv_vec_ .-= q.mass_ * p_rho_2 * kernel_gradient / r * rpq
    return nothing
end

@inline function libBalancedPressureForce!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
    kernel_value::Float64 = 0.0,
    reference_kernel_value::Float64 = 1.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    pressure_force = p.p_ + q.p_
    pressure_force += abs(pressure_force) * kernel_value / reference_kernel_value * 0.01
    p.dv_vec_ .-= q.mass_ * pressure_force * kernel_gradient / (r * p.rho_ * q.rho_) * rpq
    return nothing
end

@inline function libDensityWeightedPressureForce!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
    kernel_value::Float64 = 0.0,
    reference_kernel_value::Float64 = 1.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    pressure_force = 2 * (p.p_ * p.rho_ + q.p_ * q.rho_) / (p.rho_ + q.rho_)
    pressure_force += abs(pressure_force) * kernel_value / reference_kernel_value * 0.01
    p.dv_vec_ .-= q.mass_ * pressure_force * kernel_gradient / (r * p.rho_ * q.rho_) * rpq
    return nothing
end

@inline function libPressureExtrapolationInteraction!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_value::Float64 = 0.0,
    body_force_vec::RealVector{Dimension} = ConstVec0{Dimension}(), # this is not exactly the body force, it depends on `p`'s acceleration as well
    background_pressure::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    kernel_weight = kernel_value * q.mass_ / q.rho_
    p.sum_kernel_weight_ += kernel_weight
    p.sum_kernel_weighted_p_ +=
        kernel_weight * (max(q.p_, background_pressure) + q.rho_ * max(dot(body_force_vec, rpq), background_pressure))
    return nothing
end

@inline function libPressureExtrapolationSelfaction!(
    p::ParticleType;
    background_pressure::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}} # Dimension is not needed here
    if p.sum_kernel_weight_ > 0.0
        p.p_ = p.sum_kernel_weighted_p_ / p.sum_kernel_weight_
        p.sum_kernel_weight_ = 0.0
        p.sum_kernel_weighted_p_ = 0.0
        return nothing
    else
        p.p_ = background_pressure
        return nothing
    end
    return nothing
end
