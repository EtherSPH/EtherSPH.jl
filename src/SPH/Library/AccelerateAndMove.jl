#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 18:13:44
  @ license: MIT
  @ description:
 =#

@inline function libAccelerate!(
    p::ParticleType;
    dt::Float64 = 0.0,
    body_force_vec = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.v_vec_ .+= dt * (p.dv_vec_ .+ body_force_vec)
    return nothing
end

@inline function libMove!(
    p::ParticleType;
    dt::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.x_vec_ .+= dt * p.v_vec_
    p.dv_vec_ .= 0.0
    return nothing
end

@inline function libAccelerateAndMove!(
    p::ParticleType;
    dt::Float64 = 0.0,
    body_force_vec = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.dv_vec_ .+= body_force_vec
    p.x_vec_ .+= dt * p.v_vec_ .+ 0.5 * dt * dt * p.dv_vec_
    p.v_vec_ .+= dt * p.dv_vec_
    p.dv_vec_ .= 0.0
    return nothing
end
