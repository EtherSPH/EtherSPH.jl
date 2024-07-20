#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/19 20:26:57
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

# see `SPH/Cart/LargeDeformationContinuumModel.jl`

const dim = 2
const beam_width = 0.02
const beam_length = 10 * beam_width
const beam_inside_length = 3 * beam_width
const dr = beam_width / 20
const gap = dr
const h = 1.5 * dr
const wall_width = 3 * dr

const kernel = CubicSpline{dim}(h)

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

const dt = 0.1 * h / c_0
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

@inline function initializeContinuum!(particle_system::ParticleSystem{dim, Particle}, p::Particle)::Nothing
    EtherSPH.libLDCMInitializeContinuum!(particle_system, p; smooth_kernel = kernel)
    return nothing
end

@inline function calculateDeformationGradient!(particle_system::ParticleSystem{dim, Particle}, p::Particle;)::Nothing
    EtherSPH.libLDCMDeformationGradientEvolution!(particle_system, p; dt = dt)
    EtherSPH.libLDCMGreenLagrangeStrain!(p)
    EtherSPH.libLDCMLinearElasticityStress!(p; lambda = lambda, mu = mu)
    EtherSPH.libLDCMPioalKirchhoffStress!(p)
    return nothing
end

@inline function momentum!(particle_system::ParticleSystem{dim, Particle}, p::Particle;)::Nothing
    if p.type_ == MOVABLE_MATERIAL_TAG
        EtherSPH.libLDCMContinuumMomentum!(particle_system, p)
        EtherSPH.libLDCMContinuumArtificialStress!(
            particle_system,
            p;
            smooth_kernel = kernel,
            alpha = artificial_alpha,
            beta = artificial_beta,
            h = h,
        )
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == MOVABLE_MATERIAL_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt)
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

@info "$vtp_writer"

vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "hanging_beam_continuum_2d_"
vtp_writer.output_path_ = "example/results/hanging_beam/hanging_beam_continuum_2d"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    createNeighbourIndexList!(system)
    applyReflection!(system, initializeContinuum!)
    saveVTP(vtp_writer, system, 0, t)
    for step in ProgressBar(1:round(Int, total_time / dt))
        applyReflection!(system, calculateDeformationGradient!)
        applyReflection!(system, momentum!)
        applySelfaction!(system, accelerateAndMove!)
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end

function post()::Nothing
    PyCall.@pyinclude "example/hanging_beam/hanging_beam_continuum_2d.py"
    HangingBeam2DPostProcess = py"HangingBeam2DPostProcess"
    post_process = HangingBeam2DPostProcess(reference_gap = gap)
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
