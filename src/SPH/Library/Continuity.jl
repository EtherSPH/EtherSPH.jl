#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 16:50:56
  @ license: MIT
  @ description:
 =#

@inline function libTraditionalContinuity!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.drho_ += q.mass_ * dot(p.v_vec_ .- q.v_vec_, rpq) / r * kernel_gradient
    return nothing
end

@inline function libBalancedContinuity!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.drho_ += p.rho_ * q.mass_ * dot(p.v_vec_ .- q.v_vec_, rpq) / (r * q.rho_) * kernel_gradient
    return nothing
end
