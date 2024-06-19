#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/18 18:03:48
  @ license: MIT
  @ description:
 =#

@inline function libTraditionalThermalConduction!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
    h::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    kappa_pq = harmonicMean(p.kappa_, q.kappa_)
    heat = 2 * kappa_pq * r * kernel_gradient * (p.t_ - q.t_) / (p.rho_ * q.rho_ * (r * r + 0.01 * h * h))
    p.dt_ += q.mass_ * heat / p.cp_
    return nothing
end
