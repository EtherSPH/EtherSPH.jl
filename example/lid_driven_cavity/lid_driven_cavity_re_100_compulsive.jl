#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/17 17:38:32
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const reynolds_number = 100.0

const dim = 2
const dr = 0.01
const h = 3 * dr
const gap = dr

const kernel = CubicSpline{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const lid_length = 1.0
const wall_width = h
const lid_velocity = 1.0
const u_max = lid_velocity

const rho_0 = 1.0
const mass = rho_0 * dr^dim
const c = 15 * u_max
const mu = rho_0 * u_max * lid_length / reynolds_number
const p_0 = c^2 * rho_0 * 0.002 # background pressure

const dt = 0.1 * h / c
const t_end = 16.0
const output_dt = 100 * dt
const density_filter_dt = 10 * dt

const FLUID_TAG = 1
const WALL_TAG = 2

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c * c * (rho - rho_0) + p_0
end

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # additonal properties
    p_::Float64 = p_0
    drho_::Float64 = 0.0
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    c_::Float64 = c
    mu_::Float64 = mu
    gap_::Float64 = gap
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_value_::Float64 = 0.0
    # used for compulsive boundary
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
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ == WALL_TAG
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
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw)
        EtherSPH.libCompulsiveForce!(p, q, rpq, r; h = 0.5 * h)
        return nothing
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt)
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

const x_0 = 0.0
const y_0 = 0.0

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1
    return nothing
end

@inline function modifyBottomWall!(p::Particle)::Nothing
    modifyWall!(p)
    if p.x_vec_[1] < x_0
        p.normal_vec_ = Vector2D(1.0, 1.0) / sqrt(2)
    elseif p.x_vec_[1] > x_0 + lid_length
        p.normal_vec_ = Vector2D(-1.0, 1.0) / sqrt(2)
    else
        p.normal_vec_ = VectorY(2)
    end
    return nothing
end

@inline function modifyTopWall!(p::Particle)::Nothing
    modifyWall!(p)
    if p.x_vec_[1] < x_0
        p.normal_vec_ = Vector2D(1.0, -1.0) / sqrt(2)
    elseif p.x_vec_[1] > x_0 + lid_length
        p.normal_vec_ = Vector2D(-1.0, -1.0) / sqrt(2)
    else
        p.normal_vec_ = -VectorY(2)
    end
    p.v_vec_ = Vector2D(lid_velocity, 0.0) # ! treat as lid
    return nothing
end

@inline function modifyLeftWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = VectorX(2)
    return nothing
end

@inline function modifyRightWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = -VectorX(2)
    return nothing
end

const fluid_column = Rectangle(Vector2D(x_0, y_0), Vector2D(x_0 + lid_length, y_0 + lid_length))

const bottom_wall_column =
    Rectangle(Vector2D(x_0 - wall_width, y_0 - wall_width), Vector2D(x_0 + lid_length + wall_width, y_0))

const top_wall_column = Rectangle(
    Vector2D(x_0 - wall_width, y_0 + lid_length),
    Vector2D(x_0 + lid_length + wall_width, y_0 + lid_length + wall_width),
)

const left_wall_column = Rectangle(Vector2D(x_0 - wall_width, y_0), Vector2D(x_0, y_0 + lid_length))

const right_wall_column =
    Rectangle(Vector2D(x_0 + lid_length, y_0), Vector2D(x_0 + lid_length + wall_width, y_0 + lid_length))

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_column)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyBottomWall!)
top_wall_particles = createParticles(Particle, gap, top_wall_column; modify! = modifyTopWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyLeftWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; modify! = modifyRightWall!)

append!(particles, fluid_particles)
append!(particles, bottom_wall_particles)
append!(particles, top_wall_particles)
append!(particles, left_wall_particles)
append!(particles, right_wall_particles)

const lower = Vector2D(x_0 - wall_width, y_0 - wall_width)
const upper = Vector2D(x_0 + lid_length + wall_width, y_0 + lid_length + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle) = p.p_
@inline getVelocity(p::Particle) = p.v_vec_
@inline getNormal(p::Particle) = p.normal_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
addVector!(vtp_writer, "Normal", getNormal)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "lid_driven_cavity_re_100_compulsive_"
vtp_writer.output_path_ = "example/results/lid_driven_cavity/lid_driven_cavity_re_100_compulsive/"

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
    PyCall.@pyinclude "example/lid_driven_cavity/lid_driven_cavity.py"
    LidDrivenCavityPostProcess = py"LidDrivenCavityPostProcess"
    post_process = LidDrivenCavityPostProcess(reynolds_number = Int(reynolds_number), reference_gap = gap)
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
