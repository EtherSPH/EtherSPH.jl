#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/14 16:47:45
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

const lambda = young_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
const mu = shear_modulus

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
    x_0_vec_::Vector2D = Vector2D(0.0, 0.0) # initial position
    rho_0_::Float64 = rho_0
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    kernel_gradient_0_vec_::Vector2D = Vector2D(0.0, 0.0) # ∇₀Wᵢⱼ
    correction_mat_::Matrix2D = Matrix0(dim) # correction matrix B⁰
    deformation_mat_::Matrix2D = MatrixI(dim) # deformation matrix F
    deformation_mat_det_::Float64 = 1.0 # determinant of deformation matrix J, to apply density evolution
    deformation_mat_trace_::Float64 = 2.0 # trace of deformation matrix
    green_lagrange_strain_mat_::Matrix2D = Matrix0(dim) # Green-Lagrange strain tensor E
    piola_kirchhoff_1st_stress_mat_::Matrix2D = Matrix0(dim) # 1st Piola-Kirchhoff stress tensor P
end

@inline function initialElasticInteraction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    dw = DW(r)
    r_inv = 1 / r
    p.kernel_gradient_0_vec_ .+= dw * r_inv * q.mass_ / q.rho_ * rpq
    return nothing
end

# TODO: implement the function to calculate the deformation matrix F
