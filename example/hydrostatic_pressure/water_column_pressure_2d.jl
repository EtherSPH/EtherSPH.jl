#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/10 21:26:20
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const dim = 2

const water_width = 2.0
const water_height = 1.0
const box_width = 2.0
const box_height = 2.0
const dr = water_height / 40
const gap = dr
const h = 3 * dr
const wall_width = h

const kernel = WendlandC2{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const rho_0 = 1e3
const mass = rho_0 * dr^dim
const mu_0 = 1e-3
const gravity = 10.0
const g = Vector2D(0.0, -gravity)
const u_max = sqrt(2 * gravity * water_height)
const c_0 = 15 * u_max
const p_0 = 0.0

const FLUID_TAG = 1
const WALL_TAG = 2

const dt = 0.1 * h / c_0
const total_time = 5.0
const output_dt = 100 * dt
const density_filter_dt = 5 * dt
const total_steps = round(Int, total_time / dt)
@info "total steps: $(total_steps)"

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c_0 * c_0 * (rho - rho_0)
end

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int = FLUID_TAG
    # additional
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    drho_::Float64 = 0.0
    p_::Float64 = p_0
    c_::Float64 = c_0
    mu_::Float64 = mu_0
    gap_::Float64 = gap
    # kernel density filter
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
end

@inline function continuity!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == WALL_TAG && q.type_ == WALL_TAG
        return nothing
    end
    EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = DW(r))
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    EtherSPH.libUpdateDensity!(p; dt = dt)
    p.p_ = getPressureFromDensity(p.rho_)
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG
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
    p.type_ = FLUID_TAG
    depth = water_height - (p.x_vec_.y - y0)
    p.p_ = rho_0 * gravity * depth
    p.rho_ = p.p_ / c_0 / c_0 + rho_0
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    depth = water_height - (p.x_vec_.y - y0)
    if depth > 0.0
        p.p_ = rho_0 * gravity * depth
        p.rho_ = p.p_ / c_0 / c_0 + rho_0
    end
    return nothing
end

const fluid_column = Rectangle(Vector2D(x0, y0), Vector2D(x0 + water_width, y0 + water_height))
const bottom_wall_column =
    Rectangle(Vector2D(x0 - wall_width, y0 - wall_width), Vector2D(x0 + box_width + wall_width, y0))
const left_wall_column = Rectangle(Vector2D(x0 - wall_width, y0), Vector2D(x0, y0 + box_height))
const right_wall_column =
    Rectangle(Vector2D(x0 + box_width, y0), Vector2D(x0 + box_width + wall_width, y0 + box_height))

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_column; modify! = modifyFluid!)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; modify! = modifyWall!)

append!(particles, vcat(fluid_particles, bottom_wall_particles, left_wall_particles, right_wall_particles))

const lower = Vector2D(x0 - wall_width, y0 - wall_width)
const upper = Vector2D(x0 + box_width + wall_width, y0 + box_height + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline function getPressure(p::Particle)::Float64
    if p.type_ == FLUID_TAG
        return p.p_
    else
        return NaN
    end
end
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "water_column_pressure_2d_"
vtp_writer.output_path_ = "example/results/hydrostatic_pressure/water_column_pressure_2d"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:total_steps)
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
    PyCall.@pyinclude "example/hydrostatic_pressure/water_column_pressure_2d.py"
    WaterColumnPressurePostProcess = py"WaterColumnPressurePostProcess"
    post_process = WaterColumnPressurePostProcess()
    post_process.viewPlot()
    return nothing
end
