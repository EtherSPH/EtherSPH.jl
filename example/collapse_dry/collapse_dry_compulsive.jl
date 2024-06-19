#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/15 18:14:50
  @ license: MIT
  @ description:
 =#

using EtherSPH
using Parameters
using ProgressBars
using PyCall

const dim = 2
const dr = 0.02
const h = 3 * dr
const gap = dr

const kernel = WendlandC2{dim}(h)

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
const density_filter_dt = 30 * dt

const FLUID_TAG = 1
const WALL_TAG = 2

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c^2 * (rho - rho_0)
end

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
    sum_kernel_weighted_value_::Float64 = 0.0
    # * wall normal vector
    normal_vec_::Vector2D = Vector2D(0.0, 0.0)
end

@inline function continuity!(p::Particle, q::Particle, rpq::RealVector{dim}, r::Float64)::Nothing where {dim}
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = DW(r))
    end
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
    end
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libTraditionalPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ == WALL_TAG
        dw = DW(r)
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        EtherSPH.libCompulsiveForce!(p, q, rpq, r; h = 0.5 * h)
        return nothing
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
    end
    return nothing
end

@inline function densityFilterInteraction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libKernelAverageDensityFilterInteraction!(p, q, rpq, r; kernel_value = W(r))
    end
    return nothing
end

@inline function densityFilterSelfaction!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libKernelAverageDensityFilterSelfaction!(p; kernel_value = kernel.kernel_value_0_)
    end
    return nothing
end

const x0 = 0.0
const y0 = 0.0

@inline function modifyFluid!(p::Particle)::Nothing
    depth = water_height - p.x_vec_[2]
    p.p_ = rho_0 * gravity * depth
    p.rho_ = p.p_ / c^2 + rho_0
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1000
    return nothing
end

const fluid_column = Rectangle(Vector2D(x0, y0), Vector2D(x0 + water_width, y0 + water_height))

const bottom_wall_column = Rectangle(Vector2D(x0, y0 - wall_width), Vector2D(x0 + box_width, y0))

const left_wall_column = Rectangle(Vector2D(x0 - wall_width, y0), Vector2D(x0, y0 + box_height))
const right_wall_column =
    Rectangle(Vector2D(x0 + box_width, y0), Vector2D(x0 + box_width + wall_width, y0 + box_height))

const left_corner_wall_column = Rectangle(Vector2D(x0 - wall_width, y0 - wall_width), Vector2D(x0, y0))
const right_corner_wall_column =
    Rectangle(Vector2D(x0 + box_width, y0 - wall_width), Vector2D(x0 + box_width + wall_width, y0))

@inline function modifyBottomWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ .= VectorY(2)
    return nothing
end

@inline function modifyLeftWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ .= VectorX(2)
    return nothing
end

@inline function modifyRightWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ .= -VectorX(2)
    return nothing
end

@inline function modifyLeftCornerWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ .= (VectorY(2) + VectorX(2)) .* (1 / sqrt(2))
    return nothing
end

@inline function modifyRightCornerWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ .= (-VectorX(2) + VectorY(2)) .* (1 / sqrt(2))
    return nothing
end

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_column; modify! = modifyFluid!)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyBottomWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyLeftWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; modify! = modifyRightWall!)
left_corner_wall_particles = createParticles(Particle, gap, left_corner_wall_column; modify! = modifyLeftCornerWall!)
right_corner_wall_particles = createParticles(Particle, gap, right_corner_wall_column; modify! = modifyRightCornerWall!)

append!(
    particles,
    vcat(
        fluid_particles,
        bottom_wall_particles,
        left_wall_particles,
        right_wall_particles,
        left_corner_wall_particles,
        right_corner_wall_particles,
    ),
)

const lower = Vector2D(x0 - wall_width, y0 - wall_width)
const upper = Vector2D(x0 + box_width + wall_width, y0 + box_height + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
@inline getNormal(p::Particle)::Vector2D = p.normal_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
addVector!(vtp_writer, "Normal", getNormal)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "collapse_dry_compulsive_"
vtp_writer.output_path_ = "example/results/collapse_dry/collapse_dry_compulsive"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:round(Int, t_end / dt))
        applyInteraction!(system, continuity!)
        applySelfaction!(system, updateDensityAndPressure!)
        applyInteraction!(system, momentum!)
        applySelfaction!(system, accelerateAndMove!)
        createCellLinkList!(system)
        if step % round(Int, density_filter_dt / dt) == 0
            applyInteraction!(system, densityFilterInteraction!)
            applySelfaction!(system, densityFilterSelfaction!)
        end
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end

function post()::Nothing
    PyCall.@pyinclude "example/collapse_dry/collapse_dry.py"
    CollapseDryPostProcess = py"CollapseDryPostProcess"
    post_process = CollapseDryPostProcess(key_word = "compulsive")
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
