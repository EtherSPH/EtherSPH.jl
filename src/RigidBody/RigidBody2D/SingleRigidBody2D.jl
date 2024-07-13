#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/13 16:29:52
  @ license: MIT
  @ description:
 =#

@kwdef mutable struct SingleRigidBody2D{ParticleType <: AbstractParticle2D}
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    omega_::Float64 = 0.0
    alpha_::Float64 = 0.0
    force_vec_::Vector2D = Vector2D(0.0, 0.0)
    torque_::Float64 = 0.0
    mass_::Float64 = 0.0
    inertia_::Float64 = 0.0
    particles_::Vector{ParticleType} = ParticleType[]
end

@inline function SingleRigidBody2D(ParticleType::DataType)::SingleRigidBody2D{ParticleType}
    return SingleRigidBody2D{ParticleType}()
end

@inline function computeCenter!(
    rigid_body_2d::SingleRigidBody2D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle2D}
    sum_mass = 0.0
    rigid_body_2d.x_vec_ .= 0.0
    for p in rigid_body_2d.particles_
        sum_mass += p.mass_
        rigid_body_2d.x_vec_ .+= p.mass_ .* p.x_vec_
    end
    rigid_body_2d.x_vec_ ./= sum_mass
    return nothing
end

@inline function computeMass!(
    rigid_body_2d::SingleRigidBody2D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle2D}
    rigid_body_2d.mass_ = 0.0
    for p in rigid_body_2d.particles_
        rigid_body_2d.mass_ += p.mass_
    end
    return nothing
end

@inline function computeInertia!(
    rigid_body_2d::SingleRigidBody2D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle2D}
    rigid_body_2d.inertia_ = 0.0
    for p in rigid_body_2d.particles_
        @inbounds rx = p.x_vec_[1] - rigid_body_2d.x_vec_[1]
        @inbounds ry = p.x_vec_[2] - rigid_body_2d.x_vec_[2]
        rigid_body_2d.inertia_ += p.mass_ * (rx * rx + ry * ry)
    end
    return nothing
end

@inline function computeForceAndTorque!(
    rigid_body_2d::SingleRigidBody2D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle2D}
    rigid_body_2d.force_vec_ .= 0.0
    rigid_body_2d.torque_ = 0.0
    for p in rigid_body_2d.particles_
        rigid_body_2d.force_vec_ .+= p.mass_ * p.dv_vec_
        @inbounds rx = p.x_vec_[1] - rigid_body_2d.x_vec_[1]
        @inbounds ry = p.x_vec_[2] - rigid_body_2d.x_vec_[2]
        @inbounds rigid_body_2d.torque_ += p.mass_ * (rx * p.dv_vec_[2] - ry * p.dv_vec_[1])
    end
    return nothing
end

@inline function centerMotionLaw!(
    rigid_body_2d::SingleRigidBody2D{ParticleType};
    body_force_vec = 0.0,
)::Nothing where {ParticleType <: AbstractParticle2D}
    rigid_body_2d.dv_vec_ .= rigid_body_2d.force_vec_ / rigid_body_2d.mass_ .+ body_force_vec
    rigid_body_2d.alpha_ = rigid_body_2d.torque_ / rigid_body_2d.inertia_
    return nothing
end

@inline function accelerateAndMoveRigidBodyParticles!(
    rigid_body_2d::SingleRigidBody2D{ParticleType};
    dt::Float64 = 0.0,
)::Nothing where {ParticleType <: AbstractParticle2D}
    Threads.@threads for p in rigid_body_2d.particles_
        p.dv_vec_ .= 0.0
        @inbounds direction_vec =
            Vector2D(-(p.x_vec_[2] - rigid_body_2d.x_vec_[2]), (p.x_vec_[1] - rigid_body_2d.x_vec_[1]))
        p.v_vec_ .=
            rigid_body_2d.v_vec_ .+ rigid_body_2d.omega_ * direction_vec .+ rigid_body_2d.alpha_ * dt * direction_vec
        p.x_vec_ .+=
            rigid_body_2d.v_vec_ * dt .+ rigid_body_2d.omega_ * dt * direction_vec .+
            0.5 * rigid_body_2d.alpha_ * dt * dt * direction_vec
    end
    return nothing
end

@inline function accelerateAndMoveRigidBody!(
    rigid_body_2d::SingleRigidBody2D{ParticleType};
    dt::Float64 = 0.0,
)::Nothing where {ParticleType <: AbstractParticle2D}
    rigid_body_2d.x_vec_ .+= rigid_body_2d.v_vec_ * dt .+ 0.5 * rigid_body_2d.dv_vec_ * dt * dt
    rigid_body_2d.v_vec_ .+= rigid_body_2d.dv_vec_ * dt
    rigid_body_2d.omega_ += rigid_body_2d.alpha_ * dt
    return nothing
end

@inline function motionRigidBody!(
    rigid_body_2d::SingleRigidBody2D{ParticleType};
    dt::Float64 = 0.0,
    body_force_vec = 0.0,
)::Nothing where {ParticleType <: AbstractParticle2D}
    computeForceAndTorque!(rigid_body_2d)
    centerMotionLaw!(rigid_body_2d; body_force_vec = body_force_vec)
    accelerateAndMoveRigidBodyParticles!(rigid_body_2d; dt = dt)
    accelerateAndMoveRigidBody!(rigid_body_2d; dt = dt)
    return nothing
end
