#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/09 14:26:04
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars

const dim = 2
const beam_width = 0.2
const beam_length = 5 * beam_width
const beam_inside_length = 2 * beam_width
const dr = beam_width / 20
const gap = dr
const h = 3 * dr
const wall_width = 3 * dr

const kernel = WendlandC4{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const young_modulus = 100e6
const poisson_ratio = 0.3
const shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))
const bulk_modulus = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio))
const rho_0 = 1e3
const c_0 = sqrt(bulk_modulus / rho_0)
const mass = rho_0 * dr^dim
const g = Vector2D(0.0, -9.8)

const artificial_alpha = 2.5
const artificial_beta = 2.5

const MOVABLE_MATERIAL_TAG = 1
const FIXED_MATERIAL_TAG = 2

const dt = 0.1 * h / c_0
const total_time = 1.0
const output_dt = 100 * dt
@info "total steps: $(round(Int, total_time / dt))"

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = MOVABLE_MATERIAL_TAG
    # user defined
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    sigma_mat_::Matrix2D = Matrix0(dim) # σₐᵦ
    depsilon_mat_::Matrix2D = Matrix0(dim) # ̇εₐᵦ
    omega_mat_::Matrix2D = Matrix0(dim) # ωₐᵦ
    c_::Float64 = c_0 # speed of sound
    p_::Float64 = 0.0 # pressure
    drho_::Float64 = 0.0
end

@inline function continuityWithDeformAndRotate!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    dw = DW(r)
    EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = dw)
    r_inv = 1 / r
    m_rho = q.mass_ / q.rho_
    # deformation
    @inbounds @simd for alpha in 1:dim
        @simd for beta in 1:dim
            dv_alpha_dx_beta = m_rho * (q.v_vec_[alpha] - p.v_vec_[alpha]) * rpq[beta] * r_inv * dw
            dv_beta_dx_alpha = m_rho * (q.v_vec_[beta] - p.v_vec_[beta]) * rpq[alpha] * r_inv * dw
            p.depsilon_mat_[alpha, beta] += 0.5 * (dv_alpha_dx_beta + dv_beta_dx_alpha)
            p.omega_mat_[alpha, beta] += 0.5 * (dv_alpha_dx_beta - dv_beta_dx_alpha)
        end
    end
    return nothing
end

@inline function dirac(i::Int64, j::Int64)::Int64
    return i == j ? 1 : 0
end

@inline function updateDensityAndStress!(p::Particle)::Nothing
    EtherSPH.libUpdateDensity!(p; dt = dt)
    principle_strain = trace(p.depsilon_mat_) # ̇εᵧᵧ
    @inbounds @simd for alpha in 1:dim
        @simd for beta in 1:dim
            sigma_omega_part_1 = 0.0
            sigma_omega_part_2 = 0.0
            @simd for gamma in 1:dim
                sigma_omega_part_1 += p.sigma_mat_[alpha, gamma] * p.omega_mat_[beta, gamma]
                sigma_omega_part_2 += p.sigma_mat_[gamma, beta] * p.omega_mat_[alpha, gamma]
            end
            double_de_alpha_beta = p.depsilon_mat_[alpha, beta] + p.depsilon_mat_[beta, alpha]
            dsigma_alpha_beta =
                sigma_omega_part_1 +
                sigma_omega_part_2 +
                shear_modulus * double_de_alpha_beta +
                bulk_modulus * principle_strain * dirac(alpha, beta)
            p.sigma_mat_[alpha, beta] += dsigma_alpha_beta * dt
        end
    end
    p.depsilon_mat_ .= 0.0
    p.omega_mat_ .= 0.0
    p.p_ = -trace(p.sigma_mat_) / dim
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == MOVABLE_MATERIAL_TAG
        dw = DW(r)
        r_inv = 1 / r
        rho_p_2_inv = 1 / (p.rho_ * p.rho_)
        rho_q_2_inv = 1 / (q.rho_ * q.rho_)
        v_dot_x = dot(p.v_vec_ .- q.v_vec_, rpq)
        if v_dot_x > 0
            artificial_stress = 0.0
        else
            mean_rho = 0.5 * (p.rho_ + q.rho_)
            mean_c = 0.5 * (p.c_ + q.c_)
            phi = h * v_dot_x / (r * r + 0.01 * h * h)
            artificial_stress = (-artificial_alpha * mean_c * phi + artificial_beta * phi * phi) / mean_rho
        end
        @inbounds @simd for alpha in 1:dim
            @simd for beta in 1:dim
                p.dv_vec_[alpha] +=
                    q.mass_ *
                    (
                        rho_p_2_inv * p.sigma_mat_[alpha, beta] + rho_q_2_inv * q.sigma_mat_[alpha, beta] -
                        artificial_stress * dirac(alpha, beta)
                    ) *
                    dw *
                    rpq[beta] *
                    r_inv
            end
        end
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == MOVABLE_MATERIAL_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
    end
end

@inline function modifyBeam!(p::Particle)::Nothing
    p.type_ = MOVABLE_MATERIAL_TAG
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = FIXED_MATERIAL_TAG
    return nothing
end

const x_0 = 0.0
const y_0 = 0.0

const beam_column = Rectangle(Point2D(x_0 - beam_inside_length, y_0), Point2D(x_0 + beam_length, y_0 + beam_width))
const bottom_wall_column = Rectangle(Point2D(x_0 - beam_inside_length, y_0 - wall_width), Point2D(x_0, y_0))
const top_wall_column =
    Rectangle(Point2D(x_0 - beam_inside_length, y_0 + beam_width), Point2D(x_0, y_0 + beam_width + wall_width))

particles = Particle[]
beam_particles = createParticles(Particle, gap, beam_column; modify! = modifyBeam!, parallel = true)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyWall!, parallel = true)
top_wall_particles = createParticles(Particle, gap, top_wall_column; modify! = modifyWall!, parallel = true)
append!(particles, vcat(beam_particles, bottom_wall_particles, top_wall_particles))

const lower = Point2D(x_0 - beam_inside_length, y_0 - wall_width - beam_inside_length)
const upper = Point2D(x_0 + beam_length + beam_inside_length, y_0 + beam_width + wall_width + beam_inside_length)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
@inline getSigmaX(p::Particle)::Vector2D = Vector2D(p.sigma_mat_[1, 1], p.sigma_mat_[1, 2])
@inline getSigmaY(p::Particle)::Vector2D = Vector2D(p.sigma_mat_[2, 1], p.sigma_mat_[2, 2])
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
addVector!(vtp_writer, "SigmaX", getSigmaX)
addVector!(vtp_writer, "SigmaY", getSigmaY)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "elastic_beam_2d"
vtp_writer.output_path_ = "demo/results/elastic_beam_2d"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:round(Int, total_time / dt))
        applyInteraction!(system, continuityWithDeformAndRotate!) # continuity
        applySelfaction!(system, updateDensityAndStress!) # update density and stress
        applyInteraction!(system, momentum!) # momentum interaction
        applySelfaction!(system, accelerateAndMove!) # accelerate and move
        createCellLinkList!(system)
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end
