#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/04 20:49:29
  @ license: MIT
  @ description:
 =#

struct Cylinder{Dimension} <: Shape
    center_::Point3D
    radius_::Float64
    height_::Float64
end

const CylinderX = Cylinder{Dimension1}
const CylinderY = Cylinder{Dimension2}
const CylinderZ = Cylinder{Dimension3}

@inline function dimension(::Cylinder{Dimension})::Int64 where {Dimension}
    return Dimension3
end
@inline function axis(::Cylinder{Dimension})::Int64 where {Dimension}
    return Dimension
end

@inline function isInsideShape(point::Point3D, cylinder::CylinderX)::Bool
    if cylinder.center_[1] > point[1] || cylinder.center_[1] + height < point[1]
        return false
    end
    r = sqrt((point[2] - cylinder.center_[2])^2 + (point[3] - cylinder.center_[3])^2)
    if r > cylinder.radius_
        return false
    end
    return true
end

@inline function createParticles(
    ParticleType::DataType,
    reference_gap::Float64,
    cylinder::CylinderX;
    modify!::Function = p -> nothing,
)::Vector{ParticleType}
    ring_column_x = RingColumn{1}(cylinder.center_, reference_gap * 0.5, cylinder.radius_, cylinder.height_)
    particles = createParticles(ParticleType, reference_gap, ring_column_x; modify! = modify!)
    n_axis = round(Int, cylinder.height_ / reference_gap)
    gap_axis = cylinder.height_ / n_axis
    for i in 1:n_axis
        x_p = cylinder.center_[1] + (i - 0.5) * gap_axis
        particle = ParticleType()
        particle.x_vec_ = Point3D(x_p, cylinder.center_[2], cylinder.center_[3])
        modify!(particle)
        particle.mass_ = particle.rho_ * gap_axis * pi * (0.5 * reference_gap)^2
        push!(particles, particle)
    end
    return particles
end
