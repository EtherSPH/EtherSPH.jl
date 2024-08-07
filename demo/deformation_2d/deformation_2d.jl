#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/14 16:47:45
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars

const dim = 2
const beam_width = 0.02
const beam_length = 10 * beam_width
const beam_inside_length = 3 * beam_width
const dr = beam_width / 20
const gap = dr
const h = 1.6 * dr
const wall_width = 3 * dr

const kernel = WendlandC4{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const young_modulus = 2e6
const poisson_ratio = 0.3975
const shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))
const bulk_modulus = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio))
const rho_0 = 1e3
const c_0 = sqrt(bulk_modulus / rho_0)
const mass = rho_0 * dr^dim

const artificial_alpha = 0.1
const artificial_beta = 0.1

const lambda = young_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
const mu = shear_modulus

const MOVABLE_MATERIAL_TAG = 1
const FIXED_MATERIAL_TAG = 2

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
        return c_0 * f(x) / fl * VectorY(2) * 0.01
    end
end

const dt = 0.2 * h / c_0
const total_time = 1.0
const output_dt = 1000 * dt
@info "total steps: $(round(Int, total_time / dt))"

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = MOVABLE_MATERIAL_TAG
    # neighbour
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer()
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer()
    # user defined
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    c_::Float64 = c_0
    # need for such model ↓
    initial_neighbour_index_list_::IndexContainer = IndexContainer()
    rho_0_::Float64 = rho_
    volume_0_::Float64 = mass_ / rho_
    corrected_mat_::Matrix2D = Matrix0(dim) # corrected matrix B⁰
    kernel_gradient_0_vec_list_::Vector2DContainer = Vector2DContainer() # ∇⁰Wᵢⱼ list
    corrected_kernel_gradient_0_vec_list_::Vector2DContainer = Vector2DContainer() # B⁰∇⁰Wᵢⱼ list
    deformation_gradient_mat_::Matrix2D = MatrixI(dim) # deformation gradient matrix F
    green_lagrange_strain_mat_::Matrix2D = Matrix0(dim) # Green-Lagrange strain tensor E
    piola_kirchhoff_1st_stress_mat_::Matrix2D = Matrix0(dim) # 1st Piola-Kirchhoff stress tensor P
    piola_kirchhoff_2nd_stress_mat_::Matrix2D = Matrix0(dim) # 2nd Piola-Kirchhoff stress tensor S
    corrected_piola_kirchhoff_1st_stress_mat_::Matrix2D = Matrix0(dim) # Pᶜ
end

# after create neighbour index list ...

@inline function initialElasticParticles!(particle_system::ParticleSystem{Dimension2, Particle}, p::Particle)::Nothing
    p.rho_0_ = p.rho_
    p.volume_0_ = p.mass_ / p.rho_
    p.initial_neighbour_index_list_ = deepcopy(p.neighbour_index_list_)
    @inbounds @simd for i in eachindex(p.initial_neighbour_index_list_)
        q_index = p.initial_neighbour_index_list_[i]
        q = particle_system[q_index]
        rpq = p.neighbour_position_list_[i]
        r = p.neighbour_distance_list_[i]
        dw = DW(r)
        r_inv = 1 / r
        push!(p.kernel_gradient_0_vec_list_, dw * r_inv * rpq)
        p.corrected_mat_ .+= -q.volume_0_ * dw * r_inv * dyad(rpq, rpq)
    end
    p.corrected_mat_ .= inv(p.corrected_mat_) # B⁰
    @inbounds @simd for i in eachindex(p.initial_neighbour_index_list_)
        push!(p.corrected_kernel_gradient_0_vec_list_, dot(p.kernel_gradient_0_vec_list_[i], p.corrected_mat_))
    end
    return nothing
end

@inline function calculateDeformationMat!(particle_system::ParticleSystem{Dimension2, Particle}, p::Particle)::Nothing
    @inbounds @simd for i in eachindex(p.kernel_gradient_0_vec_list_)
        q = particle_system[p.initial_neighbour_index_list_[i]]
        p.deformation_gradient_mat_ .+=
            q.volume_0_ * dt * dyad((q.v_vec_ .- p.v_vec_), p.corrected_kernel_gradient_0_vec_list_[i])
    end
    p.rho_ = p.rho_0_ / det(p.deformation_gradient_mat_)
    p.green_lagrange_strain_mat_ .=
        0.5 * (p.deformation_gradient_mat_' * p.deformation_gradient_mat_ .- ConstMatI{dim}())
    p.piola_kirchhoff_2nd_stress_mat_ .=
        lambda * trace(p.green_lagrange_strain_mat_) * ConstMatI{dim}() .+ 2 * mu * p.green_lagrange_strain_mat_
    p.piola_kirchhoff_1st_stress_mat_ .= p.deformation_gradient_mat_ * p.piola_kirchhoff_2nd_stress_mat_
    p.corrected_piola_kirchhoff_1st_stress_mat_ .= p.piola_kirchhoff_1st_stress_mat_ * p.corrected_mat_
    return nothing
end

@inline function momentum!(particle_system::ParticleSystem{Dimension2, Particle}, p::Particle)::Nothing
    if p.type_ == MOVABLE_MATERIAL_TAG
        @inbounds @simd for i in eachindex(p.kernel_gradient_0_vec_list_)
            q = particle_system[p.initial_neighbour_index_list_[i]]
            rpq = p.x_vec_ .- q.x_vec_
            r = norm(rpq)
            v_dot_x = dot(p.v_vec_ .- q.v_vec_, rpq)
            if v_dot_x > 0.0
                artificial_stress = 0.0
            else
                mean_rho = 0.5 * (p.rho_ + q.rho_)
                mean_c = 0.5 * (p.c_ + q.c_)
                phi = h * v_dot_x / (r * r + 0.01 * h * h)
                artificial_stress = (-artificial_alpha * mean_c * phi + artificial_beta * phi * phi) / mean_rho
            end
            p.dv_vec_ .+=
                1 / p.mass_ *
                p.volume_0_ *
                q.volume_0_ *
                dot(
                    p.corrected_piola_kirchhoff_1st_stress_mat_ .+ q.corrected_piola_kirchhoff_1st_stress_mat_,
                    p.kernel_gradient_0_vec_list_[i],
                )
            p.dv_vec_ .-= q.mass_ * artificial_stress * DW(r) / r * rpq
        end
        return nothing
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == MOVABLE_MATERIAL_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt)
        return nothing
    end
    return nothing
end

@inline function modifyBeam!(p::Particle)::Nothing
    p.type_ = MOVABLE_MATERIAL_TAG
    p.v_vec_ .= verticalVelocity(p.x_vec_[1] - x_0)
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
const left_wall_column = Rectangle(
    Point2D(x_0 - beam_inside_length - wall_width, y_0 - wall_width),
    Point2D(x_0 - beam_inside_length, y_0 + beam_width + wall_width),
)

particles = Particle[]
beam_particles = createParticles(Particle, gap, beam_column; modify! = modifyBeam!, parallel = true)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyWall!, parallel = true)
top_wall_particles = createParticles(Particle, gap, top_wall_column; modify! = modifyWall!, parallel = true)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyWall!, parallel = true)
append!(particles, vcat(beam_particles, bottom_wall_particles, top_wall_particles, left_wall_particles))

const lower = Point2D(x_0 - beam_inside_length - wall_width, y_0 - wall_width - beam_inside_length)
const upper = Point2D(x_0 + beam_length + beam_inside_length, y_0 + beam_width + wall_width + beam_inside_length)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addVector!(vtp_writer, "Velocity", getVelocity)

@inline getDeformationMatDet(p::Particle)::Float64 = det(p.deformation_gradient_mat_)
addScalar!(vtp_writer, "DeformationMatDet", getDeformationMatDet)
@inline getDeformationX(p::Particle)::Vector2D =
    Vector2D(p.deformation_gradient_mat_[1], p.deformation_gradient_mat_[2])
addVector!(vtp_writer, "DeformationX", getDeformationX)
@inline getDeformationY(p::Particle)::Vector2D =
    Vector2D(p.deformation_gradient_mat_[3], p.deformation_gradient_mat_[4])
addVector!(vtp_writer, "DeformationY", getDeformationY)

@inline getVonMisesStress(p::Particle)::Float64 =
    sqrt(1.5 * ddot(p.piola_kirchhoff_2nd_stress_mat_, p.green_lagrange_strain_mat_))
addScalar!(vtp_writer, "VonMisesStress", getVonMisesStress)

vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "deformation_2d_"
vtp_writer.output_path_ = "demo/results/deformation_2d"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    createNeighbourIndexList!(system)
    applyReflection!(system, initialElasticParticles!)
    saveVTP(vtp_writer, system, 0, t)
    for step in ProgressBar(1:round(Int, total_time / dt))
        applyReflection!(system, calculateDeformationMat!)
        applyReflection!(system, momentum!)
        applySelfaction!(system, accelerateAndMove!)
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end
