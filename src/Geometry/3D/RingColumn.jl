#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/04 21:18:10
  @ license: MIT
  @ description:
 =#

struct RingColumn{Dimension} <: Shape
    center_::Point3D
    inner_radius_::Float64
    outer_radius_::Float64
    height_::Float64
end

const RingColumnX = RingColumn{Dimension1}
const RingColumnY = RingColumn{Dimension2}
const RingColumnZ = RingColumn{Dimension3}

@inline function dimension(::RingColumn{Dimension})::Int64 where {Dimension}
    return Dimension3
end
@inline function axis(::RingColumn{Dimension})::Int64 where {Dimension}
    return Dimension
end

@inline function isInsideShape(point::Point3D, ring_column::RingColumnX)::Bool
    if ring_column.center_[1] > point[1] || ring_column.center_[1] + ring_column.height_ < point[1]
        return false
    end
    r = sqrt((point[2] - ring_column.center_[2])^2 + (point[3] - ring_column.center_[3])^2)
    if r < ring_column.inner_radius_ || r > ring_column.outer_radius_
        return false
    end
    return true
end

@inline function createParticles(
    ParticleType::DataType,
    reference_gap::Float64,
    ring_column::RingColumnX;
    modify!::Function = p -> nothing,
)::Vector{ParticleType}
    n_axis = round(Int, ring_column.height_ / reference_gap)
    gap_axis = ring_column.height_ / n_axis
    n_r = round(Int, (ring_column.outer_radius_ - ring_column.inner_radius_) / reference_gap)
    gap_r = (ring_column.outer_radius_ - ring_column.inner_radius_) / n_r
    particles_along_axis = Vector{Vector{ParticleType}}(undef, n_axis)
    Threads.@threads for i_axis in 1:n_axis
        particles_along_axis[i_axis] = Vector{ParticleType}()
        x_p = ring_column.center_[1] + (i_axis - 0.5) * gap_axis
        for i in 1:n_r
            r_i = ring_column.inner_radius_ + (i - 0.5) * gap_r
            ring_length = 2 * pi * r_i
            n_theta = round(Int, ring_length / reference_gap)
            gap_theta = ring_length / n_theta
            delta_theta = 2 * pi / n_theta
            for j in 1:n_theta
                theta_j = (j - 0.5) * delta_theta
                y_p = ring_column.center_[2] + r_i * cos(theta_j)
                z_p = ring_column.center_[3] + r_i * sin(theta_j)
                push!(particles_along_axis[i_axis], ParticleType())
                particles_along_axis[i_axis][end].x_vec_ = Point3D(x_p, y_p, z_p)
                modify!(particles_along_axis[i_axis][end])
                particles_along_axis[i_axis][end].mass_ =
                    gap_axis * gap_r * gap_theta * particles_along_axis[i_axis][end].rho_
            end
        end
    end
    return vcat(particles_along_axis...)
end
