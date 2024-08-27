#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/08/16 21:40:52
  @ license: MIT
  @ description:
 =#

@inline function collectmin(
    udf::Function,
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Float64 where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @floop for i in eachindex(particle_system.particles_)
        @inbounds value = udf(particle_system.particles_[i])
        @reduce() do (min_value = Inf64; value)
            return min_value = min(min_value, value)
        end
    end
    return min_value
end

@inline function collectmax(
    udf::Function,
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Float64 where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @floop for i in eachindex(particle_system.particles_)
        @inbounds value = udf(particle_system.particles_[i])
        @reduce() do (max_value = -Inf64; value)
            return max_value = max(max_value, value)
        end
    end
    return max_value
end
