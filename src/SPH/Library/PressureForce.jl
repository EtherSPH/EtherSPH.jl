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
    p.dv_vec_ .+= -q.mass_ * p_rho_2 * kernel_gradient / r * rpq
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
    pressure_force /= p.rho_ * q.rho_
    p.dv_vec_ .+= -q.mass_ * pressure_force * kernel_gradient / r * rpq
    return nothing
end
