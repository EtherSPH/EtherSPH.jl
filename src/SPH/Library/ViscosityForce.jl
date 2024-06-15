#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 18:07:32
  @ license: MIT
  @ description:
 =#

@inline function libTraditionalViscosityForce!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
    h::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    mu_pq = harmonicMean(p.mu_, q.mu_)
    p.dv_vec_ .+=
        2 * mu_pq * q.mass_ * r * kernel_gradient / (p.rho_ * q.rho_ * (r^2 + 0.01 * h^2)) * (p.v_vec_ .- q.v_vec_)
    return nothing
end
