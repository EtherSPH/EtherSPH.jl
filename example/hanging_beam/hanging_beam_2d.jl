#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/25 20:10:09
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const dim = 2

const beam_length = 0.2
const beam_width = 0.02
const beam_inside_length = 0.06
const dr = beam_width / 20
const gap = dr
const h = 2.6 * dr
const wall_width = 3 * dr

const zeta = 3.8 # ζ, a magic number which is used to control the numerical stability, proposed by Zhang et al.

const kernel = CubicSpline{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const rho_0 = 1e3
const mass = rho_0 * dr^dim
const young_modulus = 2e6
const poisson_ratio = 0.3975
const shear_modulus = young_modulus / 2 / (1 + poisson_ratio)
const bulk_modulus = young_modulus / 3 / (1 - 2 * poisson_ratio)
const c = sqrt(bulk_modulus / rho_0)

const l = beam_length
const k = 1.875 / l
const kl = k * l
const hh = beam_width
const omega = sqrt(young_modulus * hh^2 * k^4 / 12 / rho_0 / (1 - poisson_ratio^2))
@info "Period time 𝒯 = $(2 * pi / omega)"

@inline function f(x::Float64)::Float64
    return (sin(kl) + sinh(kl)) * (cos(k * x) - cosh(k * x)) - (cos(kl) + cosh(kl)) * (sin(k * x) - sinh(k * x))
end

const fl = f(l)

@inline function verticalVelocity(x::Float64)::Vector2D
    if x < 0.0
        return Vector0(2)
    else
        return c * f(x) / fl * VectorY(2) * 0.01
    end
end

const dt = 0.1 * h / c
const total_time = 1.0
const output_dt = 500 * dt
@info "total steps: $(round(Int, total_time / dt))"

const MATERIAL_MOVABLE_TAG = 1
const MATERIAL_FIXED_TAG = 2

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have properties:
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = MATERIAL_MOVABLE_TAG
    # additional properties:
    p_::Float64 = 0.0
    drho_::Float64 = 0.0
    v_vec_::Vector2D = Vector0(2)
    dv_vec_::Vector2D = Vector0(2)
    gap_::Float64 = gap
    # * new, dv_shear_vec_
    dv_shear_vec_::Vector2D = Vector0(2)
end

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c^2 * (rho - rho_0)
end

@inline function continuity!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    dw = DW(r)
    libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = dw)
    p.dv_shear_vec_ .+=
        2 * zeta * shear_modulus * dt * q.mass_ * dw * dot(p.v_vec_ .- q.v_vec_, rpq) / (p.rho_ * q.rho_ * r * r * r) *
        rpq
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    libUpdateDensity!(p; dt = dt)
    p.p_ = getPressureFromDensity(p.rho_)
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == MATERIAL_MOVABLE_TAG
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
        return nothing
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == MATERIAL_MOVABLE_TAG
        p.dv_vec_ .+= p.dv_shear_vec_
        libAccelerateAndMove!(p; dt = dt)
    end
    return nothing
end

const x0 = 0.0
const y0 = 0.0

@inline function initializeBeam!(p::Particle)::Nothing
    p.type_ = MATERIAL_MOVABLE_TAG
    p.v_vec_ .= verticalVelocity(p.x_vec_[1] - x0)
    return nothing
end

@inline function initializeWall!(p::Particle)::Nothing
    p.type_ = MATERIAL_FIXED_TAG
    return nothing
end

const beam_column = Rectangle(Vector2D(x0 - beam_inside_length, y0), Vector2D(x0 + beam_length, y0 + beam_width))
const bottom_wall_column = Rectangle(Vector2D(x0 - beam_inside_length, y0 - wall_width), Vector2D(x0, y0))
const top_wall_column =
    Rectangle(Vector2D(x0 - beam_inside_length, y0 + beam_width), Vector2D(x0, y0 + beam_width + wall_width))
const left_wall_column = Rectangle(
    Vector2D(x0 - beam_inside_length - wall_width, y0 - wall_width),
    Vector2D(x0 - beam_inside_length, y0 + beam_width + wall_width),
)

particles = Particle[]
beam_particles = createParticles(Particle, gap, beam_column; modify! = initializeBeam!)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = initializeWall!)
top_wall_particles = createParticles(Particle, gap, top_wall_column; modify! = initializeWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = initializeWall!)

append!(particles, vcat(beam_particles, bottom_wall_particles, top_wall_particles, left_wall_particles))

const lower = Vector2D(x0 - beam_inside_length - wall_width, y0 - wall_width - beam_inside_length)
const upper = Vector2D(x0 + beam_length + beam_inside_length, y0 + beam_width + beam_inside_length + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle) = p.p_
@inline getVelocity(p::Particle) = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "hanging_beam_2d_"
vtp_writer.output_path_ = "example/results/hanging_beam/hanging_beam_2d"

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
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end

function post()::Nothing
    PyCall.@pyinclude "example/hanging_beam/hanging_beam_2d.py"
    HangingBeam2DPostProcess = py"HangingBeam2DPostProcess"
    post_process = HangingBeam2DPostProcess(reference_gap = gap)
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
