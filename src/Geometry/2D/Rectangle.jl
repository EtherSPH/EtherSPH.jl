#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 21:44:53
  @ license: MIT
  @ description:
 =#

struct Rectangle <: Shape
    lower_::Point2D
    upper_::Point2D
end

@inline dimension(::Rectangle)::Int64 = Dimension2

@inline function isInsideShape(point::Point2D, rectangle::Rectangle)::Bool
    for i in 1:Dimension2
        if point[i] < rectangle.lower_[i] || point[i] > rectangle.upper_[i]
            return false
        end
    end
    return true
end

@inline function createParticles(
    ParticleType::DataType,
    reference_gap::Float64,
    rectangle::Rectangle;
    modify!::Function = p -> nothing,
    parallel::Bool = true,
)::Vector{ParticleType}
    n_x = round(Int, (rectangle.upper_[1] - rectangle.lower_[1]) / reference_gap)
    gap_x = (rectangle.upper_[1] - rectangle.lower_[1]) / n_x
    n_y = round(Int, (rectangle.upper_[2] - rectangle.lower_[2]) / reference_gap)
    gap_y = (rectangle.upper_[2] - rectangle.lower_[2]) / n_y
    particles = Vector{ParticleType}(undef, n_x * n_y)
    if parallel == true
        Threads.@threads for index in 1:(n_x * n_y)
            i = div(index - 1, n_y) + 1
            j = mod(index - 1, n_y) + 1
            x_p = rectangle.lower_[1] + (i - 0.5) * gap_x
            y_p = rectangle.lower_[2] + (j - 0.5) * gap_y
            particles[index] = ParticleType()
            particles[index].x_vec_ = Point2D(x_p, y_p)
            modify!(particles[index])
            particles[index].mass_ = gap_x * gap_y * particles[index].rho_
        end
    else
        for index in 1:(n_x * n_y)
            i = div(index - 1, n_y) + 1
            j = mod(index - 1, n_y) + 1
            x_p = rectangle.lower_[1] + (i - 0.5) * gap_x
            y_p = rectangle.lower_[2] + (j - 0.5) * gap_y
            particles[index] = ParticleType()
            particles[index].x_vec_ = Point2D(x_p, y_p)
            modify!(particles[index])
            particles[index].mass_ = gap_x * gap_y * particles[index].rho_
        end
    end
    return particles
end
