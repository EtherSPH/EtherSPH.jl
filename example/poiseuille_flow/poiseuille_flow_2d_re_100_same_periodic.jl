#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/26 21:01:00
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const x0 = 0.0
const y0 = 0.0

const dim = 2
const dr = 1e-2 / 50
const gap = dr
const h = 3 * dr

const kernel = WendlandC2{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const reynolds_number = 100
const pipe_width = 1e-2
const pipe_length = 4 * pipe_width
const wall_width = 3 * dr
const reference_length = pipe_width

const rho_0 = 1e3
const mass = rho_0 * dr^dim
const mu_0 = 1e-3
const u_mean = reynolds_number * mu_0 / (rho_0 * pipe_width)
const fx = 12 * mu_0 * u_mean / (pipe_width^2)
const ax = fx / rho_0
const g = Vector2D(ax, 0.0)

@inline function getUxFromY(y::Float64)::Float64
    return fx * 0.5 / mu_0 * y * (pipe_width - y)
end

const u_max = getUxFromY(0.5 * pipe_width)
const c_0 = 10 * u_max
const p_0 = rho_0 * c_0^2 * 0.01
@info "u max: $u_max, c_0: $c_0"

@inline function getPressureFromDensity(rho::Float64)::Float64
    return p_0 + c_0 * c_0 * (rho - rho_0)
end

const dt = 0.1 * h / c_0
const total_time = 50.0
const output_dt = 200 * dt
const density_filter_dt = 20 * dt
@info "total step: $(round(Int, total_time / dt))"

const FLUID_TAG = 1
const WALL_TAG = 2

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
    c_::Float64 = c_0
    mu_::Float64 = mu_0
    gap_::Float64 = gap
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_value_::Float64 = 0.0
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
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = h / 2)
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
        # ! periodic boundary condition here !
        if p.x_vec_[1] > x0 + pipe_length
            p.x_vec_[1] -= pipe_length
        elseif p.x_vec_[1] < x0
            p.x_vec_[1] += pipe_length
        end
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

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1000
    return nothing
end

const fluid_column = Rectangle(Vector2D(x0, y0), Vector2D(x0 + pipe_length, y0 + pipe_width))
const bottom_wall_column = Rectangle(Vector2D(x0, y0 - wall_width), Vector2D(x0 + pipe_length, y0))
const top_wall_column =
    Rectangle(Vector2D(x0, y0 + pipe_width), Vector2D(x0 + pipe_length, y0 + pipe_width + wall_width))

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_column)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyWall!)
top_wall_particles = createParticles(Particle, gap, top_wall_column; modify! = modifyWall!)

append!(particles, vcat(fluid_particles, bottom_wall_particles, top_wall_particles))

const lower = Vector2D(x0, y0 - wall_width)
const upper = Vector2D(x0 + pipe_length, y0 + pipe_width + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

# ! set the periodic boundary along x direction
setPeriodicBoundaryAlongX!(system)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "poiseuille_flow_2d_re_100_same_periodic_"
vtp_writer.output_path_ = "example/results/poiseuille_flow/poiseuille_flow_2d_re_100_same_periodic"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:round(Int, total_time / dt))
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
    PyCall.@pyinclude "example/poiseuille_flow/poiseuille_flow_2d.py"
    PoiseuilleFlow2DPostProcess = py"PoiseuilleFlow2DPostProcess"
    post_process = PoiseuilleFlow2DPostProcess(
        reference_gap = gap,
        key_word = "2d_re_$(reynolds_number)_same_periodic",
        ax = ax,
        analytical_order = 20,
    )
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
