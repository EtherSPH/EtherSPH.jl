#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/21 00:43:55
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const prandtl_number = 0.706
const rayleigh_number = 4.7e4

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

const r_outer = 1.0
const ratio_io = 2.6
const r_inner = r_outer / ratio_io
const reference_length = r_outer - r_inner

const gravity = 1.0
const g = Vector2D(0.0, -gravity)
const beta = 0.05
const mu_0 = 1e-3
const kappa = 1.0

const t_outer = 0.0
const t_inner = 1.0
const delta_t = t_inner - t_outer
const mean_t = 0.5 * (t_outer + t_inner)
const t_0 = 0.2 * delta_t + t_outer # magic number here

const cp = prandtl_number * kappa / mu_0
const rho_0 = sqrt(rayleigh_number * mu_0 * kappa / gravity / beta / (reference_length^3) / delta_t / cp)
const mass = rho_0 * dr^dim
const c_0 = 2.0
const p_0 = 0.025 * rho_0 * c_0^2
const alpha = mu_0 / rho_0 / prandtl_number

const dt = 0.1 * h / c_0
const t_end = 100.0
const output_dt = 100 * dt
const density_filter_dt = 5 * dt

const FLUID_TAG = 1
const THERMOSTATIC_WALL_TAG = 2

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c_0^2 * (rho - rho_0) + p_0
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
    c_::Float64 = c_0
    mu_::Float64 = mu_0
    gap_::Float64 = gap
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_value_::Float64 = 0.0
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
const center = Vector2D(x_0, y_0)

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    p.t_ = t_0
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = THERMOSTATIC_WALL_TAG
    p.mu_ *= 1.0
    p.kappa_ *= 1.0 # i have to admit that this is a magic number, direct harmonic mean is not good
    return nothing
end

@inline function modifyInnerWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = normalize(Vector2D(p.x_vec_.x - x_0, p.x_vec_.y - y_0))
    p.t_ = t_inner
    return nothing
end

@inline function modifyOuterWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = normalize(Vector2D(x_0 - p.x_vec_.x, y_0 - p.x_vec_.y))
    p.t_ = t_outer
    return nothing
end

const fluid_ring = Ring(center, r_inner, r_outer)
const inner_wall_ring = Ring(center, r_inner - wall_width, r_inner)
const outer_wall_ring = Ring(center, r_outer, r_outer + wall_width)

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_ring; modify! = modifyFluid!)
inner_wall_particles = createParticles(Particle, gap, inner_wall_ring; modify! = modifyInnerWall!)
outer_wall_particles = createParticles(Particle, gap, outer_wall_ring; modify! = modifyOuterWall!)

append!(particles, vcat(fluid_particles, inner_wall_particles, outer_wall_particles))

const lower = Vector2D(x_0 - r_outer - h, y_0 - r_outer - h)
const upper = Vector2D(x_0 + r_outer + h, y_0 + r_outer + h)
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
vtp_writer.file_name_ = "natural_convection_ring_ratio_2.6_ra_4.7e4_"
vtp_writer.output_path_ = "example/results/natural_convection_ring/natural_convection_ring_ratio_2.6_ra_4.7e4"

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
    PyCall.@pyinclude "example/natural_convection_ring/natural_convection_ring.py"
    NaturalConvectionRingPostProcess = py"NaturalConvectionRingPostProcess"
    post_process = NaturalConvectionRingPostProcess(key_word = "ratio_2.6_ra_4.7e4", rayleigh_number = "4.7e4")
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end