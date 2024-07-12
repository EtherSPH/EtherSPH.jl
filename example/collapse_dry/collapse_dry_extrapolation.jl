#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/11 16:02:46
  @ license: MIT
  @ description:
 =#

using EtherSPH
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
const density_filter_dt = 10 * dt

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
    sum_kernel_weighted_rho_::Float64 = 0.0
    sum_kernel_weighted_p_::Float64 = 0.0
end

@inline function continuityAndPressureExpolation!(
    p::Particle,
    q::Particle,
    rpq::RealVector{dim},
    r::Float64,
)::Nothing where {dim}
    if p.type_ == FLUID_TAG
        dw = DW(r)
        EtherSPH.libBalancedContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ == WALL_TAG && q.type_ == FLUID_TAG
        w = W(r)
        kernel_weight = w * q.mass_ / q.rho_
        p.sum_kernel_weight_ += kernel_weight
        p.sum_kernel_weighted_p_ += kernel_weight * (max(q.p_, 0.0) + q.rho_ * max(dot(g, rpq), 0.0))
        return nothing
    end
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
        return nothing
    elseif p.type_ == WALL_TAG
        if p.sum_kernel_weight_ > 0.0
            p.p_ = p.sum_kernel_weighted_p_ / p.sum_kernel_weight_
        else
            p.p_ = 0.0
        end
        p.sum_kernel_weight_ = 0.0
        p.sum_kernel_weighted_p_ = 0.0
        return nothing
    end
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libBalancedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = h / 2)
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

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1000
    return nothing
end

const fluid_column = Rectangle(Vector2D(x0, y0), Vector2D(x0 + water_width, y0 + water_height))

const bottom_wall_column =
    Rectangle(Vector2D(x0 - wall_width, y0 - wall_width), Vector2D(x0 + box_width + wall_width, y0))

const left_wall_column = Rectangle(Vector2D(x0 - wall_width, y0), Vector2D(x0, y0 + box_height))

const right_wall_column =
    Rectangle(Vector2D(x0 + box_width, y0), Vector2D(x0 + box_width + wall_width, y0 + box_height))

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_column)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; modify! = modifyWall!)

append!(particles, fluid_particles)
append!(particles, bottom_wall_particles)
append!(particles, left_wall_particles)
append!(particles, right_wall_particles)

const lower = Vector2D(x0 - wall_width, y0 - wall_width)
const upper = Vector2D(x0 + box_width + wall_width, y0 + box_height + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "collapse_dry_extropolation_"
vtp_writer.output_path_ = "example/results/collapse_dry/collapse_dry_extrapolation"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:round(Int, t_end / dt))
        applyInteraction!(system, continuityAndPressureExpolation!)
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
    post_process = CollapseDryPostProcess(key_word = "extrapolation")
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
