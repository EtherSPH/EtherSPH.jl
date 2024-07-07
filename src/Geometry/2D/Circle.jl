#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/04 21:02:57
  @ license: MIT
  @ description:
 =#

struct Circle <: Shape
    center_::Point2D
    radius_::Float64
end

@inline dimension(::Circle)::Int64 = Dimension2

@inline function isInsideShape(point::Point2D, circle::Circle)::Bool
    return norm(point .- circle.center_) < circle.radius_
end

@inline function createParticles(
    ParticleType::DataType,
    reference_gap::Float64,
    circle::Circle;
    modify!::Function = p -> nothing,
    parallel::Bool = false,
)::Vector{ParticleType}
    ring = Ring(circle.center_, reference_gap * 0.5, circle.radius_)
    particles = createParticles(ParticleType, reference_gap, ring; modify! = modify!, parallel = parallel)
    center_particle = ParticleType()
    center_particle.x_vec_ .= circle.center_
    modify!(center_particle)
    center_particle.mass_ = pi * (0.5 * reference_gap)^2 * center_particle.rho_
    push!(particles, center_particle)
    return particles
end
