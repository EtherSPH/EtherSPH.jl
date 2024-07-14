#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/11 17:12:12
  @ license: MIT
  @ description:
 =#

# * default AbstractParticle must at least have fields:
# * - x_vec_::RealVector
# * - rho_::Float64
# * - mass_::Float64
# * - type_::Int64
# * - otherwise you may need override the following functions:
# * - getPosition(particle::AbstractParticle)::RealVector
# * - getDensity(particle::AbstractParticle)::Float64
# * - getMass(particle::AbstractParticle)::Float64
# * - getType(particle::AbstractParticle)::Int64

abstract type AbstractParticle{Dimension} end

const AbstractParticle2D = AbstractParticle{Dimension2}
const AbstractParticle3D = AbstractParticle{Dimension3}

@inline dimension(::Type{<:AbstractParticle{Dimension}}) where {Dimension} = Dimension

@inline function getPosition(particle::AbstractParticle{Dimension})::RealVector{Dimension} where {Dimension}
    return particle.x_vec_
end
@inline function getDensity(particle::AbstractParticle{Dimension})::Float64 where {Dimension}
    return particle.rho_
end
@inline function getMass(particle::AbstractParticle{Dimension})::Float64 where {Dimension}
    return particle.mass_
end
@inline function getType(particle::AbstractParticle{Dimension})::Int64 where {Dimension}
    return particle.type_
end

mutable struct ParticleSystem{Dimension, ParticleType <: AbstractParticle{Dimension}}
    particles_::Vector{ParticleType}
    cell_link_list_::CellLinkList{Dimension}
end

const ParticleSystem2D{ParticleType <: AbstractParticle2D} = ParticleSystem{Dimension2, ParticleType}
const ParticleSystem3D{ParticleType <: AbstractParticle3D} = ParticleSystem{Dimension3, ParticleType}

Base.length(
    particle_system::ParticleSystem{Dimension, ParticleType},
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = length(particle_system.particles_)
Base.size(
    particle_system::ParticleSystem{Dimension, ParticleType},
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = size(particle_system.particles_)
Base.eltype(::ParticleSystem{Dimension, ParticleType}) where {Dimension, ParticleType <: AbstractParticle{Dimension}} =
    ParticleType
Base.getindex(
    particle_system::ParticleSystem{Dimension, ParticleType},
    i::Int64,
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = particle_system.particles_[i]
Base.setindex!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    value::ParticleType,
    i::Int64,
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = particle_system.particles_[i] = value
Base.push!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    value::ParticleType,
) where {Dimension, ParticleType} = push!(particle_system.particles_, value)
Base.pop!(
    particle_system::ParticleSystem{Dimension, ParticleType},
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = pop!(particle_system.particles_)
Base.eachindex(
    particle_system::ParticleSystem{Dimension, ParticleType},
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = Base.OneTo(length(particle_system.particles_))
Base.append!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    particles::AbstractVector{ParticleType},
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = append!(particle_system.particles_, particles)
Base.resize!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    new_size::Int64,
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = resize!(particle_system.particles_, new_size)
Base.deleteat!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    i::Int64,
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = deleteat!(particle_system.particles_, i)
# ! warning: deleteat!'s range should be unique and sorted
Base.deleteat!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    range::Vector{Int64},
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = deleteat!(particle_system.particles_, range)
Base.show(
    io::IO,
    particle_system::ParticleSystem{Dimension, ParticleType},
) where {Dimension, ParticleType <: AbstractParticle{Dimension}} = print(
    io,
    "ParticleSystem{",
    Dimension,
    "}:\n",
    "  particles: $(length(particle_system.particles_)) $(eltype(particle_system.particles_)) particles\n",
    "  cell link list: $(particle_system.cell_link_list_)\n",
)

@inline function ParticleSystem(
    ParticleType::DataType,
    cell_link_list::CellLinkList{Dimension},
)::ParticleSystem{Dimension} where {Dimension}
    particles = Vector(ParticleType[])
    return ParticleSystem{Dimension, ParticleType}(particles, cell_link_list)
end

@inline function ParticleSystem(
    ParticleType::DataType,
    reference_gap::Float64,
    lower::Point{Dimension},
    upper::Point{Dimension},
)::ParticleSystem{Dimension} where {Dimension}
    cell_link_list = CellLinkList(reference_gap, lower, upper)
    return ParticleSystem(ParticleType, cell_link_list)
end
