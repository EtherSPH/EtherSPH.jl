#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 21:14:26
  @ license: MIT
  @ description:
 =#

@inline function libUpdateDensity!(
    p::ParticleType;
    dt::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.rho_ += dt * p.drho_
    p.drho_ = 0.0
    return nothing
end
