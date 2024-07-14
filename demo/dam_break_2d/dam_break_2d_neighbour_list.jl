#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/14 20:48:55
  @ license: MIT
  @ description:
 =#

#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 22:08:32
  @ license: MIT
  @ description:
 =#

using EtherSPH
using Parameters
using ProgressBars

const dim = 2
const dr = 0.02
const h = 3 * dr
const gap = dr

const kernel = CubicSpline{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end
@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const water_width = 1.0
const water_height = 2.0
const box_width = 4.0
const box_height = 3.0
const wall_width = h

const rho_0 = 1e3
const mass = rho_0 * dr^dim
const gravity = 9.8
const g = Vector2D(0.0, -gravity)
const c = 10 * sqrt(2 * water_height * gravity)
const mu = 1e-3

const dt = 0.1 * h / c
const t_end = 4.0
const output_dt = 100 * dt
const density_filter_dt = 10 * dt

const FLUID_TAG = 1
const WALL_TAG = 2

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # additonal properties
    p_::Float64 = 0.0
    drho_::Float64 = 0.0
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    c_::Float64 = c
    mu_::Float64 = mu
    gap_::Float64 = dr
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
    # neighbour
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer()
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer()
end

# for better written code, real 'fortran'
@inline x(p::Particle)::Vector2D = p.x_vec_
@inline rho(p::Particle)::Float64 = p.rho_
@inline m(p::Particle)::Float64 = p.mass_
@inline v(p::Particle)::Vector2D = p.v_vec_
@inline P(p::Particle)::Float64 = p.p_

@inline function continuity!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    p.drho_ += m(q) * dot(v(p) - v(q), rpq) / r * DW(r)
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    p.rho_ += p.drho_ * dt
    p.drho_ = 0.0
    p.p_ = c^2 * (rho(p) - rho_0)
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG
        dW = DW(r)
        p.dv_vec_ .+= -m(q) * (P(p) / rho(p)^2 + P(q) / rho(q)^2) * dW / r * rpq
        mu_pq = 2 * p.mu_ * q.mu_ / (p.mu_ + q.mu_)
        p.dv_vec_ .+= 2 * mu_pq * m(q) * r * dW / rho(q) / rho(p) * (r^2 + 0.01 * (h / 2)^2) * (v(p) - v(q))
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        p.x_vec_ .+= v(p) * dt .+ 0.5 * (p.dv_vec_ .+ g) * dt^2
        p.v_vec_ .+= (p.dv_vec_ .+ g) * dt
        p.dv_vec_ .= 0.0
    end
    return nothing
end

@inline function densityFilter!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        p.sum_kernel_weight_ += W(r) * m(q) / rho(q)
        p.sum_kernel_weighted_rho_ += W(r) * m(q)
    end
    return nothing
end

@inline function densityFilter!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        p.sum_kernel_weight_ += kernel.kernel_value_0_ * m(p) / rho(p)
        p.sum_kernel_weighted_rho_ += kernel.kernel_value_0_ * m(p)
        p.rho_ = p.sum_kernel_weighted_rho_ / p.sum_kernel_weight_
        p.sum_kernel_weight_ = 0.0
        p.sum_kernel_weighted_rho_ = 0.0
    end
    return nothing
end

function createRectangleParticles(
    ParticleType::DataType,
    x0::Float64,
    y0::Float64,
    width::Float64,
    height::Float64,
    reference_dr::Float64;
    modifyOnParticle!::Function = (p) -> nothing,
)::Vector{ParticleType}
    particles = Vector{ParticleType}()                 # create particles, first an empty array
    n_along_x = Int64(width / reference_dr |> round)
    n_along_y = Int64(height / reference_dr |> round)
    dx = width / n_along_x
    dy = height / n_along_y
    for i in 1:n_along_x, j in 1:n_along_y
        particle = ParticleType()
        x = x0 + (i - 0.5) * dx # calculate the x position, in the center of each square
        y = y0 + (j - 0.5) * dy # calculate the y position, in the center of each square
        particle.x_vec_ = Vector2D(x, y)
        modifyOnParticle!(particle) # modify the particle as you wish
        particle.mass_ = particle.rho_ * dx * dy # calculate the mass of the particle
        push!(particles, particle) # add the particle to the array
    end
    return particles
end

const x0 = 0.0
const y0 = 0.0

particles = Particle[]

fluid_particles = createRectangleParticles(Particle, x0, y0, water_width, water_height, dr)
append!(particles, fluid_particles)

@inline function initialWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1000
    return nothing
end

bottom_wall_particles = createRectangleParticles(
    Particle,
    x0 - wall_width,
    y0 - wall_width,
    box_width + 2 * wall_width,
    wall_width,
    dr;
    modifyOnParticle! = initialWall!,
)
append!(particles, bottom_wall_particles)

left_wall_particles = createRectangleParticles(
    Particle,
    x0 - wall_width,
    y0,
    wall_width,
    box_height,
    dr;
    modifyOnParticle! = initialWall!,
)
append!(particles, left_wall_particles)

right_wall_particles =
    createRectangleParticles(Particle, x0 + box_width, y0, wall_width, box_height, dr; modifyOnParticle! = initialWall!)
append!(particles, right_wall_particles)

lower = Vector2D(x0 - wall_width, y0 - wall_width)
upper = Vector2D(x0 + box_width + wall_width, y0 + box_height + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "dam_break_2d_neighbour_list"
vtp_writer.output_path_ = "demo/results/dam_break_2d_neighbour_list"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createNeighbourIndexList!(system)
    for step in ProgressBar(1:round(Int, t_end / dt))
        applyInteractionWithNeighbours!(system, continuity!) # continuity interaction
        applySelfaction!(system, updateDensityAndPressure!) # update the density and pressure
        applyInteractionWithNeighbours!(system, momentum!) # momentum interaction
        applySelfaction!(system, accelerateAndMove!) # accelerate and move
        createNeighbourIndexList!(system) # create neighbour index list
        if step % round(Int, density_filter_dt / dt) == 0
            applyInteractionWithNeighbours!(system, densityFilter!) # density filter
            applySelfaction!(system, densityFilter!) # density filter
        end
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end
