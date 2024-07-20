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
        2 * mu_pq * q.mass_ * r * kernel_gradient / (p.rho_ * q.rho_ * (r * r + 0.01 * h * h)) * (p.v_vec_ .- q.v_vec_)
    return nothing
end

@inline function libArtificialViscosityForce!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
    alpha::Float64 = 0.1,
    beta::Float64 = 0.1,
    h::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    v_dot_x = dot(p.v_vec_ .- q.v_vec_, rpq)
    if v_dot_x > 0
        return nothing
    end
    mean_rho = 0.5 * (p.rho_ + q.rho_)
    mean_c = 0.5 * (p.c_ + q.c_)
    phi = h * v_dot_x / (r * r + 0.01 * h * h)
    artificial_stress = (-alpha * mean_c + beta * phi) * phi / mean_rho
    p.dv_vec_ .-= q.mass_ / (p.rho_ * q.rho_ * r) * artificial_stress * kernel_gradient * rpq
    return nothing
end
