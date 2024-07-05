#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/29 20:27:04
  @ license: MIT
  @ description:
 =#

struct Cuboid <: Shape
    lower_::Point3D
    upper_::Point3D
end

@inline dimension(::Cuboid)::Int64 = Dimension3

@inline function isInsideShape(point::Point3D, cuboid::Cuboid)::Bool
    for i in 1:Dimension3
        if point[i] < cuboid.lower_[i] || point[i] > cuboid.upper_[i]
            return false
        end
    end
    return true
end

@inline function createParticles(
    ParticleType::DataType,
    reference_gap::Float64,
    cuboid::Cuboid;
    modify!::Function = p -> nothing,
)::Vector{ParticleType}
    n_x = round(Int, (cuboid.upper_[1] - cuboid.lower_[1]) / reference_gap)
    gap_x = (cuboid.upper_[1] - cuboid.lower_[1]) / n_x
    n_y = round(Int, (cuboid.upper_[2] - cuboid.lower_[2]) / reference_gap)
    gap_y = (cuboid.upper_[2] - cuboid.lower_[2]) / n_y
    n_z = round(Int, (cuboid.upper_[3] - cuboid.lower_[3]) / reference_gap)
    gap_z = (cuboid.upper_[3] - cuboid.lower_[3]) / n_z
    particles = Vector{ParticleType}(undef, n_x * n_y * n_z)
    Threads.@threads for index in 1:(n_x * n_y * n_z)
        i = div(index - 1, n_y * n_z) + 1
        j = div(mod(index - 1, n_y * n_z), n_z) + 1
        k = mod(index - 1, n_z) + 1
        x_p = cuboid.lower_[1] + (i - 0.5) * gap_x
        y_p = cuboid.lower_[2] + (j - 0.5) * gap_y
        z_p = cuboid.lower_[3] + (k - 0.5) * gap_z
        particles[index] = ParticleType()
        particles[index].x_vec_ = Point3D(x_p, y_p, z_p)
        modify!(particles[index])
        particles[index].mass_ = gap_x * gap_y * gap_z * particles[index].rho_
    end
    return particles
end
