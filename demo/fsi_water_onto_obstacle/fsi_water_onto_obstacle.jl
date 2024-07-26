#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/24 17:31:12
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars

const x_0 = 0.0
const y_0 = 0.0

const dim = 2

const water_width = 1.46
const water_height = 2.92
const box_width = 4 * water_width
const box_height = 1.2 * water_height
const beam_width = 0.12
const beam_height = 0.8
const dr = 0.02
const gap = dr
const gap_s = gap / 2
const water_h = 3 * dr
const beam_h = 1.5 * gap_s
const wall_width = 6 * dr
const beam_start_x = 2 * water_width + x_0
const beam_start_y = 0.0 + y_0

const ker_f = WendlandC2{dim}(water_h)
const ker_s = WendlandC2{dim}(beam_h)

const w0_f = ker_f.kernel_value_0_
const w0_s = ker_s.kernel_value_0_

@inline function WF(r::Float64)::Float64
    return kernelValue(r, ker_f)
end
@inline function DWF(r::Float64)::Float64
    return kernelGradient(r, ker_f)
end

@inline function WS(r::Float64)::Float64
    return kernelValue(r, ker_s)
end
@inline function DWS(r::Float64)::Float64
    return kernelGradient(r, ker_s)
end

const rho_f = 1e3
const c_f = 200.0
const mu_f = 1e-3
const rho_s = 2.5e3
const young_modulus = 20e6
const poisson_ratio = 0.2
const shear_modulus = young_modulus / (2 * (1 + poisson_ratio))
const bulk_modulus = young_modulus / (3 * (1 - 2 * poisson_ratio))
const c_s = sqrt(shear_modulus / rho_s)
const lambda_s = bulk_modulus - 2 / 3 * shear_modulus
const mu_s = shear_modulus
const artificial_alpha = 2.5
const artificial_beta = 2.5
const gravity = 9.81
const g = Vector2D(0.0, -gravity)
@info "fluid sound speed: $c_f, solid sound speed: $c_s"

const dt = min(0.1 * water_h / c_f, 0.2 * beam_h / c_s) * 0.4
const total_time = 5.0
const output_dt = 800 * dt
const density_filter_dt = 30 * dt
const total_step = round(Int, total_time / dt)
@info "time step size: $dt, total step: $total_step, output time:$(round(Int, total_time / output_dt))"

const WATER_TAG = 1
const WALL_TAG = 2
const BEAM_MOVABLE_TAG = 3
const BEAM_FIXED_TAG = 4

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    mass_::Float64 = 0.0
    rho_::Float64 = rho_f
    type_::Int64 = WATER_TAG
    # neighbour info
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer()
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer()
    # velocity and acceleration
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    drho_::Float64 = 0.0
    # shared properties
    p_::Float64 = 0.0
    gap_::Float64 = gap
    sum_kernel_weight_::Float64 = 0.0
    # fluid properties
    c_::Float64 = c_f
    mu_::Float64 = mu_f
    sum_kernel_weighted_rho_::Float64 = 0.0
    # solid properties
    sum_kernel_weighted_p_::Float64 = 0.0
    dv_vec_copy_::Vector2D = Vector2D(0.0, 0.0)
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

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c_f * c_f * (rho - rho_f)
end

@inline function continuityAndPressureExtrapolation!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == WATER_TAG && q.type_ == WATER_TAG
        dw = DWF(r)
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ == WATER_TAG && q.type_ in (WALL_TAG, BEAM_MOVABLE_TAG, BEAM_FIXED_TAG)
        dw = DWF(r)
        EtherSPH.libBalancedContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ in (WALL_TAG, BEAM_MOVABLE_TAG, BEAM_FIXED_TAG) && q.type_ == WATER_TAG
        kernel_weight = WF(r) * q.mass_ / q.rho_
        p.sum_kernel_weight_ += kernel_weight
        p.sum_kernel_weighted_p_ += kernel_weight * (max(q.p_, 0.0) + q.rho_ * max(dot(g .- p.dv_vec_copy_, rpq), 0.0))
        return nothing
    end
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    if p.type_ == WATER_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
        return nothing
    else
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

@inline function momentumTotal!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == WATER_TAG && q.type_ in (WATER_TAG, WALL_TAG)
        w = WF(r)
        dw = DWF(r)
        rw = WF(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libBalancedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * water_h)
        return nothing
    elseif p.type_ == WATER_TAG && q.type_ in (BEAM_MOVABLE_TAG, BEAM_FIXED_TAG)
        dw = DWF(r)
        w = WF(r)
        rw = WF(0.5 * (p.gap_ + q.gap_))
        pressure = (p.p_ * p.rho_ + q.p_ * q.rho_) / (p.rho_ + q.rho_)
        pressure += abs(pressure) * 0.01 * w / rw
        p.dv_vec_ .-= 2 * q.mass_ / (p.rho_ * q.rho_ * r) * pressure * dw * rpq
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * water_h)
        return nothing
    elseif p.type_ == BEAM_MOVABLE_TAG && q.type_ == WATER_TAG
        dw = DWF(r)
        w = WF(r)
        rw = WF(0.5 * (p.gap_ + q.gap_))
        pressure = (p.p_ * p.rho_ + q.p_ * q.rho_) / (p.rho_ + q.rho_)
        pressure += abs(pressure) * 0.01 * w / rw
        p.dv_vec_ .-= 2 * q.mass_ / (p.rho_ * q.rho_ * r) * pressure * dw * rpq
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * water_h)
        return nothing
    end
    return nothing
end

@inline function densityFilterInteraction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == WATER_TAG && q.type_ == WATER_TAG
        EtherSPH.libKernelAverageDensityFilterInteraction!(p, q, rpq, r; kernel_value = WF(r))
    end
    return nothing
end

@inline function densityFilterSelfaction!(p::Particle)::Nothing
    if p.type_ == WATER_TAG
        EtherSPH.libKernelAverageDensityFilterSelfaction!(p; kernel_value = w0_f)
    end
    return nothing
end

@inline function initializeBeam!(particle_system::ParticleSystem{dim, Particle}, p::Particle)::Nothing
    EtherSPH.libLDCMInitializeContinuum!(particle_system, p; smooth_kernel = ker_s)
    return nothing
end

@inline function calculateDeformationGradient!(particle_system::ParticleSystem{dim, Particle}, p::Particle;)::Nothing
    EtherSPH.libLDCMDeformationGradientEvolution!(particle_system, p; dt = dt)
    EtherSPH.libLDCMGreenLagrangeStrain!(p)
    EtherSPH.libLDCMLinearElasticityStress!(p; lambda = lambda_s, mu = mu_s)
    EtherSPH.libLDCMPioalKirchhoffStress!(p)
    return nothing
end

@inline function momentumBeam!(particle_system::ParticleSystem{dim, Particle}, p::Particle;)::Nothing
    if p.type_ == BEAM_MOVABLE_TAG
        EtherSPH.libLDCMContinuumMomentum!(particle_system, p)
        EtherSPH.libLDCMContinuumArtificialStress!(
            particle_system,
            p;
            smooth_kernel = ker_s,
            alpha = artificial_alpha,
            beta = artificial_beta,
            h = beam_h,
        )
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == WATER_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
        return nothing
    elseif p.type_ == BEAM_MOVABLE_TAG
        p.dv_vec_copy_ .= p.dv_vec_
        EtherSPH.libAccelerateAndMove!(p; dt = dt)
        return nothing
    end
    return nothing
end

@inline function modifyWater!(p::Particle)::Nothing
    p.type_ = WATER_TAG
    p.rho_ = rho_f
    p.c_ = c_f
    p.mu_ = mu_f
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.rho_ = rho_f
    p.c_ = c_f
    p.mu_ *= 1000
    return nothing
end

@inline function modifyBeamMovable!(p::Particle)::Nothing
    p.type_ = BEAM_MOVABLE_TAG
    return nothing
end

@inline function modifyBeamFixed!(p::Particle)::Nothing
    p.type_ = BEAM_FIXED_TAG
    return nothing
end

@inline function modifyBeam!(p::Particle)::Nothing
    y = p.x_vec_[2]
    p.rho_ = rho_s
    p.c_ = c_s
    p.mu_ = mu_s
    p.gap_ = gap_s
    if y < beam_start_y
        modifyBeamFixed!(p)
    else
        modifyBeamMovable!(p)
    end
    return nothing
end

const water_column = Rectangle(Point2D(x_0, y_0), Point2D(x_0 + water_width, y_0 + water_height))
const bottom_wall_column_1 = Rectangle(Point2D(x_0 - wall_width, y_0 - wall_width), Point2D(x_0 + beam_start_x, y_0))
const bottom_wall_column_2 =
    Rectangle(Point2D(x_0 + beam_start_x + beam_width, y_0 - wall_width), Point2D(x_0 + box_width + wall_width, y_0))
const left_wall_column = Rectangle(Point2D(x_0 - wall_width, y_0), Point2D(x_0, y_0 + box_height))
const right_wall_column =
    Rectangle(Point2D(x_0 + box_width, y_0), Point2D(x_0 + box_width + wall_width, y_0 + box_height))
const beam_column = Rectangle(
    Point2D(beam_start_x, beam_start_y - wall_width),
    Point2D(beam_start_x + beam_width, beam_start_y + beam_height),
)

particles = Particle[]

water_particles = createParticles(Particle, gap, water_column; modify! = modifyWater!)
bottom_wall_particles_1 = createParticles(Particle, gap, bottom_wall_column_1; modify! = modifyWall!)
bottom_wall_particles_2 = createParticles(Particle, gap, bottom_wall_column_2; modify! = modifyWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; modify! = modifyWall!)
beam_particles = createParticles(Particle, 0.5 * gap, beam_column; modify! = modifyBeam!)

append!(
    particles,
    vcat(
        water_particles,
        bottom_wall_particles_1,
        bottom_wall_particles_2,
        left_wall_particles,
        right_wall_particles,
        beam_particles,
    ),
)
total_system = ParticleSystem(
    Particle,
    water_h,
    Vector2D(x_0 - wall_width, y_0 - wall_width),
    Vector2D(x_0 + box_width + wall_width, y_0 + box_height + wall_width),
)
beam_system = ParticleSystem(
    Particle,
    beam_h,
    Vector2D(beam_start_x - wall_width, beam_start_y - wall_width),
    Vector2D(beam_start_x + beam_width + wall_width, beam_start_y + beam_height + wall_width),
)

append!(total_system, particles)
append!(beam_system, beam_particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
addScalar!(vtp_writer, "Pressure", getPressure)
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addVector!(vtp_writer, "Velocity", getVelocity)
@inline function getVonMisesStress(p::Particle)::Float64
    if p.type_ < 3
        return NaN
    else
        return sqrt(1.5 * ddot(p.piola_kirchhoff_2nd_stress_mat_, p.green_lagrange_strain_mat_))
    end
    return nothing
end
addScalar!(vtp_writer, "VonMisesStress", getVonMisesStress)
@inline function getAcceleration(p::Particle)::Vector2D
    if p.type_ == BEAM_MOVABLE_TAG
        return p.dv_vec_copy_
    else
        return Vector2D(NaN, NaN)
    end
end
addVector!(vtp_writer, "Acceleration", getAcceleration)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "fsi_water_onto_obstacle_"
vtp_writer.output_path_ = "demo/results/fsi_water_onto_obstacle"

@info "$vtp_writer"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    createNeighbourIndexList!(beam_system)
    applyReflection!(beam_system, initializeBeam!)
    createNeighbourIndexList!(total_system)
    saveVTP(vtp_writer, total_system, 0, t)
    for step in ProgressBar(1:total_step)
        applyInteractionWithNeighbours!(total_system, continuityAndPressureExtrapolation!)
        applySelfaction!(total_system, updateDensityAndPressure!)
        applyInteractionWithNeighbours!(total_system, momentumTotal!)
        applyReflection!(beam_system, calculateDeformationGradient!)
        applyReflection!(beam_system, momentumBeam!)
        applySelfaction!(total_system, accelerateAndMove!)
        createNeighbourIndexList!(total_system)
        if step % round(Int, density_filter_dt / dt) == 0
            applyInteractionWithNeighbours!(total_system, densityFilterInteraction!)
            applySelfaction!(total_system, densityFilterSelfaction!)
        end
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, total_system, step, t)
        end
        t += dt
    end
    return nothing
end

function post()::Nothing
    splitParticlesType(
        vtp_writer;
        type_name_dict = Dict(1 => "Water", 2 => "Wall", 3 => "Beam_Move", 4 => "Beam_Fixed"),
    )
    return nothing
end
