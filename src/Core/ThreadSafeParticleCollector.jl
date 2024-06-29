#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/28 16:04:20
  @ license: MIT
  @ description:
 =#

# * ==================== ParticleContainer struct =================== * #
mutable struct ParticleContainer{Dimension, ParticleType <: AbstractParticle{Dimension}}
    particles_::Vector{ParticleType}
    length_::Int64
end

@inline function ParticleContainer(ParticleType::DataType, dimension::Int64)::ParticleContainer{dimension, ParticleType}
    return ParticleContainer{dimension, ParticleType}(Vector{ParticleType}(), 0)
end

@inline function ParticleContainer(ParticleType::DataType)::ParticleContainer{dimension(ParticleType), ParticleType}
    return ParticleContainer{dimension(ParticleType), ParticleType}(Vector{ParticleType}(), 0)
end

@inline function Base.length(
    particle_container::ParticleContainer{Dimension, ParticleType},
)::Int64 where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    return particle_container.length_
end

@inline function Base.getindex(
    particle_container::ParticleContainer{Dimension, ParticleType},
    i::Int64,
)::ParticleType where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @inbounds return particle_container.particles_[i]
end

@inline function Base.getindex(
    particle_container::ParticleContainer{Dimension, ParticleType},
    range::UnitRange{Int64},
)::Vector{ParticleType} where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @inbounds return particle_container.particles_[range]
end

@inline function Base.setindex!(
    particle_container::ParticleContainer{Dimension, ParticleType},
    value::ParticleType,
    i::Int64,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @inbounds particle_container.particles_[i] = value
    return nothing
end

@inline function Base.eachindex(
    particle_container::ParticleContainer{Dimension, ParticleType},
)::UnitRange{Int64} where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    return Base.OneTo(length(particle_container))
end

@inline function capacity(
    particle_container::ParticleContainer{Dimension, ParticleType},
)::Int64 where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    return length(particle_container.particles_)
end

@inline function setlength!(
    particle_container::ParticleContainer{Dimension, ParticleType},
    length::Int64,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    particle_container.length_ = length
    return nothing
end

@inline function reset!(
    particle_container::ParticleContainer{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    setlength!(particle_container, 0)
    return nothing
end

@inline function Base.resize!(
    particle_container::ParticleContainer{Dimension, ParticleType},
    new_capacity::Int64,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    if new_capacity > capacity(particle_container)
        resize!(particle_container.particles_, new_capacity)
    end
    return nothing
end

@inline function Base.push!(
    particle_container::ParticleContainer{Dimension, ParticleType},
    value::ParticleType,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    if length(particle_container) == capacity(particle_container)
        resize!(particle_container, capacityExpandPolicy(capacity(particle_container)))
    end
    setlength!(particle_container, length(particle_container) + 1)
    setindex!(particle_container, value, length(particle_container))
    return nothing
end

# * ==================== ThreadSafeParticleCollector struct =================== * #
mutable struct ThreadSafeParticleCollector{Dimension, ParticleType <: AbstractParticle{Dimension}}
    particles_containers_::Vector{ParticleContainer{Dimension, ParticleType}} # nthread * particle number each thread
end

@inline function ThreadSafeParticleCollector(
    ParticleType::DataType,
    dimension::Int64,
)::ThreadSafeParticleCollector{dimension, ParticleType}
    n_threads = Threads.nthreads()
    return ThreadSafeParticleCollector{dimension, ParticleType}([
        ParticleContainer{dimension, ParticleType}(Vector{ParticleType}(), 0) for _ in 1:n_threads
    ],)
end

@inline function ThreadSafeParticleCollector(
    ParticleType::DataType,
)::ThreadSafeParticleCollector{dimension(ParticleType), ParticleType}
    n_threads = Threads.nthreads()
    return ThreadSafeParticleCollector{dimension(ParticleType), ParticleType}([
        ParticleContainer{dimension(ParticleType), ParticleType}(Vector{ParticleType}(), 0) for _ in 1:n_threads
    ],)
end

@inline function Base.length(
    thread_safe_particle_collector::ThreadSafeParticleCollector{Dimension, ParticleType},
)::Int64 where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    return sum(length(container) for container in thread_safe_particle_collector.particles_containers_)
end

@inline function Base.show(
    io::IO,
    thread_safe_particle_collector::ThreadSafeParticleCollector{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    print(io, "ThreadSafeParticleCollector{")
    print(io, Dimension)
    print(io, ", ")
    print(io, ParticleType)
    print(io, "}: \n")
    for i in 1:length(thread_safe_particle_collector.particles_containers_)
        print(
            io,
            "    thread id $(i): number of particles $(length(thread_safe_particle_collector.particles_containers_[i]))\n",
        )
    end
    return nothing
end

@inline function Base.push!(
    thread_safe_particle_collector::ThreadSafeParticleCollector{Dimension, ParticleType},
    value::ParticleType,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    thread_id = Threads.threadid()
    push!(thread_safe_particle_collector.particles_containers_[thread_id], value)
    return nothing
end

@inline function reset!(
    thread_safe_particle_collector::ThreadSafeParticleCollector{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    Threads.@threads for container in thread_safe_particle_collector.particles_containers_
        reset!(container)
    end
    return nothing
end

@inline function Base.append!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    thread_safe_particle_collector::ThreadSafeParticleCollector{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    for i in 1:length(thread_safe_particle_collector.particles_containers_)
        @inbounds length_of_container = length(thread_safe_particle_collector.particles_containers_[i])
        if length_of_container > 0
            @inbounds append!(
                particle_system.particles_,
                thread_safe_particle_collector.particles_containers_[i][1:length_of_container],
            )
        end
        @inbounds reset!(thread_safe_particle_collector.particles_containers_[i])
    end
    return nothing
end
