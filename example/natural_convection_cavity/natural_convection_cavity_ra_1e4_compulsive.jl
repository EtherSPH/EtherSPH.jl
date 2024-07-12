#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/18 17:38:47
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const prandtl_number = 0.71
const rayleigh_number = 1e4

const dim = 2
const dr = 0.01
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

const cavity_length = 1.0
const rho_0 = 1.0
const mass = rho_0 * dr^dim
const c = 3.5
const nu_0 = prandtl_number / sqrt(rayleigh_number)
const mu_0 = nu_0 * rho_0
const p_0 = 0.02 * rho_0 * c^2
const alpha = nu_0 / prandtl_number

const gravity = 1.0
const g = Vector2D(0.0, -gravity)
const t_left = 1.0
const t_right = 0.0
const delta_t = t_left - t_right
const mean_t = 0.5 * (t_left + t_right)
const t_0 = mean_t
const kappa = 1.0
const cp = prandtl_number * kappa / mu_0
const beta = rayleigh_number * nu_0 * alpha / gravity / delta_t / cavity_length^3

const dt = 0.1 * h / c
const t_end = 50.0
const output_dt = 100 * dt
const density_filter_dt = 30 * dt

const FLUID_TAG = 1
const WALL_TAG = 2 # adiabatic wall
const THERMOSTATIC_WALL_TAG = 3

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c * c * (rho - rho_0) + p_0
end

@inline function bodyForceVectorByBoussiensqApproximation(t::Float64)::Vector2D
    return beta * (t_0 - t) * g
end

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have properties:
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # additional properties:
    # * first group, flow field properties
    p_::Float64 = p_0
    drho_::Float64 = 0.0
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    c_::Float64 = c
    mu_::Float64 = mu_0
    gap_::Float64 = dr
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
    # * second group, wall properties
    normal_vec_::Vector2D = Vector2D(0.0, 0.0)
    # * third group, temperature field properties
    # ρcₚ∂T/∂t = ∇·(k∇T)
    t_::Float64 = t_0
    dt_::Float64 = 0.0
    kappa_::Float64 = kappa
    cp_::Float64 = cp
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
        return nothing
    else
        return nothing
    end
    return nothing
end

@inline function continuity!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = DW(r))
        return nothing
    end
    return nothing
end

@inline function momentumAndThermalConduction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libTraditionalPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_value = w,
            kernel_gradient = dw,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        EtherSPH.libTraditionalThermalConduction!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
    elseif p.type_ == FLUID_TAG && q.type_ == WALL_TAG
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libTraditionalPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_value = w,
            kernel_gradient = dw,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        EtherSPH.libCompulsiveForce!(p, q, rpq, r; h = 0.5 * h)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ == THERMOSTATIC_WALL_TAG
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libTraditionalPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_value = w,
            kernel_gradient = dw,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        EtherSPH.libTraditionalThermalConduction!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        EtherSPH.libCompulsiveForce!(p, q, rpq, r; h = 0.5 * h)
        return nothing
    end
    return nothing
end

@inline function accelerateAndMoveAndHeated!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = bodyForceVectorByBoussiensqApproximation(p.t_))
        EtherSPH.libUpdateTemperature!(p; dt = dt)
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

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    p.t_ = t_0
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.mu_ *= 1000
    p.kappa_ *= 0.48 # i have to admit that this is a magic number, direct harmonic mean is not good
    p.t_ = NaN
    return nothing
end

@inline function modifyBottomWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.type_ = WALL_TAG
    p.normal_vec_ = Vector2D(0.0, 1.0)
    return nothing
end

@inline function modifyTopWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.type_ = WALL_TAG
    p.normal_vec_ = Vector2D(0.0, -1.0)
    return nothing
end

@inline function modifyLeftWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.type_ = THERMOSTATIC_WALL_TAG
    p.t_ = t_left
    if p.x_vec_[2] < y_0
        p.normal_vec_ = Vector2D(1.0, 1.0) / sqrt(2)
    elseif p.x_vec_[2] > y_0 + cavity_length
        p.normal_vec_ = Vector2D(1.0, -1.0) / sqrt(2)
    else
        p.normal_vec_ = VectorX(2)
    end
    return nothing
end

@inline function modifyRightWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.type_ = THERMOSTATIC_WALL_TAG
    p.t_ = t_right
    if p.x_vec_[2] < y_0
        p.normal_vec_ = Vector2D(-1.0, 1.0) / sqrt(2)
    elseif p.x_vec_[2] > y_0 + cavity_length
        p.normal_vec_ = Vector2D(-1.0, -1.0) / sqrt(2)
    else
        p.normal_vec_ = -VectorX(2)
    end
    return nothing
end

const fluid_column = Rectangle(Vector2D(x_0, y_0), Vector2D(x_0 + cavity_length, y_0 + cavity_length))

const bottom_wall_column = Rectangle(Vector2D(x_0, y_0 - wall_width), Vector2D(x_0 + cavity_length, y_0))
const top_wall_column =
    Rectangle(Vector2D(x_0, y_0 + cavity_length), Vector2D(x_0 + cavity_length, y_0 + cavity_length + wall_width))
const left_wall_column =
    Rectangle(Vector2D(x_0 - wall_width, y_0 - wall_width), Vector2D(x_0, y_0 + cavity_length + wall_width))
const right_wall_column = Rectangle(
    Vector2D(x_0 + cavity_length, y_0 - wall_width),
    Vector2D(x_0 + cavity_length + wall_width, y_0 + cavity_length + wall_width),
)

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_column; (modify!) = modifyFluid!)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; (modify!) = modifyBottomWall!)
top_wall_particles = createParticles(Particle, gap, top_wall_column; (modify!) = modifyTopWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; (modify!) = modifyLeftWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; (modify!) = modifyRightWall!)

append!(
    particles,
    vcat(fluid_particles, bottom_wall_particles, top_wall_particles, left_wall_particles, right_wall_particles),
)

const lower = Vector2D(x_0 - wall_width, y_0 - wall_width)
const upper = Vector2D(x_0 + cavity_length + wall_width, y_0 + cavity_length + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle) = p.p_
@inline getTemperature(p::Particle) = p.t_
@inline getVelocity(p::Particle) = p.v_vec_
@inline getNormal(p::Particle) = p.normal_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addScalar!(vtp_writer, "Temperature", getTemperature)
addVector!(vtp_writer, "Velocity", getVelocity)
addVector!(vtp_writer, "Normal", getNormal)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "natural_convection_cavity_ra_1e4_compulsive_"
vtp_writer.output_path_ = "example/results/natural_convection_cavity/natural_convection_cavity_ra_1e4_compulsive/"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:round(Int, t_end / dt))
        applyInteraction!(system, continuity!)
        applySelfaction!(system, updateDensityAndPressure!)
        applyInteraction!(system, momentumAndThermalConduction!)
        applySelfaction!(system, accelerateAndMoveAndHeated!)
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
    PyCall.@pyinclude "example/natural_convection_cavity/natural_convection_cavity.py"
    NaturalConvectionCavityPostProcess = py"NaturalConvectionCavityPostProcess"
    post_process = NaturalConvectionCavityPostProcess(rayleigh_number = "1e4", reference_gap = dr)
    post_process.viewPlot()
    post_process.referencePlot()
    nu = post_process.calculateNusseltNumber()
    @info "Nusselt number: $nu at Ra = 1e4, however, this post process is not accurate"
    return nothing
end
