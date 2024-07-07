#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/21 00:35:18
  @ license: MIT
  @ description:
 =#

struct Ring <: Shape
    center_::Point2D
    inner_radius_::Float64
    outer_radius_::Float64
end

@inline dimension(::Ring)::Int64 = Dimension2

@inline function isInsideShape(point::Point2D, ring::Ring)::Bool
    r = norm(point .- ring.center_)
    if r < ring.inner_radius_ || r > ring.outer_radius_
        return false
    end
    return true
end

@inline function createParticles(
    ParticleType::DataType,
    reference_gap::Float64,
    ring::Ring;
    modify!::Function = p -> nothing,
    parallel::Bool = false
)::Vector{ParticleType}
    n_r = round(Int, (ring.outer_radius_ - ring.inner_radius_) / reference_gap)
    gap_r = (ring.outer_radius_ - ring.inner_radius_) / n_r
    particles = Vector{ParticleType}()
    for i in 1:n_r
        r_i = ring.inner_radius_ + (i - 0.5) * gap_r
        ring_length = 2 * pi * r_i
        n_theta = round(Int, ring_length / reference_gap)
        gap_theta = ring_length / n_theta
        delta_theta = 2 * pi / n_theta
        for j in 1:n_theta
            theta_j = (j - 0.5) * delta_theta
            x_p = ring.center_[1] + r_i * cos(theta_j)
            y_p = ring.center_[2] + r_i * sin(theta_j)
            push!(particles, ParticleType())
            particles[end].x_vec_ = Point2D(x_p, y_p)
            modify!(particles[end])
            particles[end].mass_ = gap_r * gap_theta * particles[end].rho_
        end
    end
    return particles
end
