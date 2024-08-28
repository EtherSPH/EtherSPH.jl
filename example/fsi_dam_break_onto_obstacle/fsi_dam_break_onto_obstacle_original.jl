#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/08/28 14:06:07
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressMeter
using PyCall
using CSV
using DataFrames

const x_0 = 0.0
const y_0 = 0.0

const dim = 2

const water_width = 0.146
const water_height = 0.292
const box_width = 4 * water_width
const box_height = 1.2 * water_height
const beam_width = 0.012
const beam_height = 0.08
const dr_f = 0.002
const gap_f = dr_f
const dr_s = dr_f / 2
const gap_s = dr_s
const h_f = 3 * dr_f
const h_s = 2.6 * dr_s
const wall_width = 6 * dr_f
const beam_start_x = 2 * water_width + x_0
const beam_start_y = 0.0 + y_0

const ker_f = WendlandC2{dim}(h_f)
const ker_s = WendlandC2{dim}(h_s)

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
const c_f = 60.0
const u_max = c_f / 10
const mu_f = 1e-3
const nu_f = mu_f / rho_f
const rho_s = 2.5e3
const young_modulus = 1e6
const poisson_ratio = 0.2
const shear_modulus = young_modulus / (2 * (1 + poisson_ratio))
const bulk_modulus = young_modulus / (3 * (1 - 2 * poisson_ratio))
const c_s = sqrt(shear_modulus / rho_s)
const lambda_s = bulk_modulus - 2 / 3 * shear_modulus
const mu_s = shear_modulus
const artificial_alpha_f = 0.1
const artificial_beta_f = 0.1
const artificial_alpha_s = 2.5
const artificial_beta_s = 2.5
const gravity = 9.81
const g = Vector2D(0.0, -gravity)
@info "fluid sound speed: $c_f, solid sound speed: $c_s"

const reference_dt = min(0.1 * h_f / c_f, 0.1 * h_s / c_s) * 2
const total_time = 2.0
const output_step_interval = 400
const density_filter_step_interval = 30
const reference_total_step = round(Int, total_time / reference_dt)
@info "reference_dt: $reference_dt, reference_total_step: $reference_total_step"

const FLUID_TAG = 1
const WALL_TAG = 2
const BEAM_MOVABLE_TAG = 3
const BEAM_FIXED_TAG = 4

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    mass_::Float64 = 0.0
    rho_::Float64 = rho_f
    type_::Int64 = FLUID_TAG
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
    gap_::Float64 = gap_f
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
    # time step
    step_time_::Float64 = Inf64
end

@inline function timeStep(p::Particle)::Float64
    return p.step_time_
end

@inline function calculateStepTime!(p::Particle)::Nothing
    if p.type_ in (BEAM_MOVABLE_TAG, BEAM_FIXED_TAG)
        dv_dt = norm(p.dv_vec_)
        if dv_dt > 0
            p.step_time_ = 0.2 * min(sqrt(h_s / dv_dt), h_s / (norm(p.v_vec_) + c_s))
        end
        return nothing
    elseif p.type_ == FLUID_TAG
        p.step_time_ = 0.3 * min(0.6 * h_f / (norm(p.v_vec_) + c_f), 0.25 * min(h_f / u_max, h_f * h_f / nu_f))
        return nothing
    end
    return nothing
end

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c_f * c_f * (rho - rho_f)
end

# * ============================ SPH Equation Begin ============================

@inline function continuityAndPressureExtrapolation!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        dw = DWF(r)
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, BEAM_MOVABLE_TAG, BEAM_FIXED_TAG)
        dw = DWF(r)
        EtherSPH.libBalancedContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ in (WALL_TAG, BEAM_MOVABLE_TAG, BEAM_FIXED_TAG) && q.type_ == FLUID_TAG
        EtherSPH.libPressureExtrapolationInteraction!(
            p,
            q,
            rpq,
            r;
            kernel_value = WF(r),
            body_force_vec = g .- p.dv_vec_copy_,
        )
        return nothing
    end
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle; dt::Float64 = 0.0)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
        return nothing
    elseif p.type_ in (WALL_TAG, BEAM_MOVABLE_TAG, BEAM_FIXED_TAG)
        EtherSPH.libPressureExtrapolationSelfaction!(p;)
        return nothing
    end
    return nothing
end

@inline function momentumTotal!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        w = WF(r)
        dw = DWF(r)
        rw = WF(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libTraditionalPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_value = w,
            kernel_gradient = dw,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h_f)
        EtherSPH.libArtificialViscosityForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            alpha = artificial_alpha_f,
            beta = artificial_beta_f,
            h = h_f,
        )
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, BEAM_MOVABLE_TAG, BEAM_FIXED_TAG)
        w = WF(r)
        dw = DWF(r)
        rw = WF(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libBalancedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_value = w,
            kernel_gradient = dw,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h_f)
        return nothing
    elseif p.type_ == BEAM_MOVABLE_TAG && q.type_ == FLUID_TAG
        w = WF(r)
        dw = DWF(r)
        rw = WF(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libBalancedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_value = w,
            kernel_gradient = dw,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h_f)
        return nothing
    end
    return nothing
end

@inline function densityFilterInteraction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libKernelAverageDensityFilterInteraction!(p, q, rpq, r; kernel_value = WF(r))
    end
    return nothing
end

@inline function densityFilterSelfaction!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libKernelAverageDensityFilterSelfaction!(p; kernel_value = w0_f)
    end
    return nothing
end

@inline function initializeSolid!(particle_system::ParticleSystem{dim, Particle}, p::Particle)::Nothing
    EtherSPH.libLDCMInitializeContinuum!(particle_system, p; smooth_kernel = ker_s)
    return nothing
end

@inline function calculateDeformationGradient!(
    particle_system::ParticleSystem{dim, Particle},
    p::Particle;
    dt::Float64 = 0.0,
)::Nothing
    EtherSPH.libLDCMDeformationGradientEvolution!(particle_system, p; dt = dt)
    EtherSPH.libLDCMGreenLagrangeStrain!(p)
    EtherSPH.libLDCMLinearElasticityStress!(p; lambda = lambda_s, mu = mu_s)
    EtherSPH.libLDCMPioalKirchhoffStress!(p)
    return nothing
end

@inline function momentumSolid!(particle_system::ParticleSystem{dim, Particle}, p::Particle;)::Nothing
    if p.type_ == BEAM_MOVABLE_TAG
        EtherSPH.libLDCMContinuumMomentum!(particle_system, p)
        EtherSPH.libLDCMContinuumArtificialStress!(
            particle_system,
            p;
            smooth_kernel = ker_s,
            alpha = artificial_alpha_s,
            beta = artificial_beta_s,
            h = h_s,
        )
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle; dt::Float64 = 0.0, t::Float64 = 0.0)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
        return nothing
    elseif p.type_ == BEAM_MOVABLE_TAG
        p.dv_vec_copy_ .= p.dv_vec_
        EtherSPH.libAccelerateAndMove!(p; dt = dt)
        return nothing
    elseif p.type_ == BEAM_FIXED_TAG
        p.dv_vec_ .= 0.0
        return nothing
    end
    return nothing
end

# * ============================ SPH Equation End ============================

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    p.rho_ = rho_f
    p.c_ = c_f
    p.mu_ = mu_f
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.rho_ = rho_f
    p.c_ = c_f
    p.mu_ *= 100000
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
    p.mu_ *= 100000
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

water_particles = createParticles(Particle, gap_f, water_column; modify! = modifyFluid!)
bottom_wall_particles_1 = createParticles(Particle, gap_f, bottom_wall_column_1; modify! = modifyWall!)
bottom_wall_particles_2 = createParticles(Particle, gap_f, bottom_wall_column_2; modify! = modifyWall!)
left_wall_particles = createParticles(Particle, gap_f, left_wall_column; modify! = modifyWall!)
right_wall_particles = createParticles(Particle, gap_f, right_wall_column; modify! = modifyWall!)
beam_particles = createParticles(Particle, gap_s, beam_column; modify! = modifyBeam!)

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
    h_f,
    Vector2D(x_0 - wall_width, y_0 - wall_width),
    Vector2D(x_0 + box_width + wall_width, y_0 + box_height + wall_width),
)
beam_system = ParticleSystem(
    Particle,
    h_s,
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
setStepDigit!(vtp_writer, 4)
setFileName!(vtp_writer, "fsi_dam_break_onto_obstacle_original")
setOutputPath!(vtp_writer, "example/results/fsi_dam_break_onto_obstacle/fsi_dam_break_onto_obstacle_original")

@info "$vtp_writer"

@inline function findMonitorPointIndex()::Int64
    index = 1
    for i in eachindex(beam_particles)
        p = beam_particles[i]
        x = p.x_vec_[1]
        y = p.x_vec_[2]
        if abs(x - beam_start_x) <= 1.01 * 0.5 * gap_s && abs(y - beam_start_y - beam_height) <= 1.01 * 0.5 * gap_s
            index = i
        end
    end
    return index
end

const monitor_point_index = findMonitorPointIndex()
@info "monitor_point_index: $monitor_point_index"
const csv_file_name = joinpath(vtp_writer.output_path_, "monitor_point.csv")

function main()::Nothing
    t = 0.0
    dt = 0.0
    step = 0
    ts = Float64[]
    xs = Float64[]
    ys = Float64[]
    assurePathExist(vtp_writer)
    createNeighbourIndexList!(beam_system)
    applyReflection!(beam_system, initializeSolid!)
    createNeighbourIndexList!(total_system)
    saveVTP(vtp_writer, total_system, 0, t)
    progress = ProgressThresh(0.0; desc = "Computing...", showspeed = true)
    while t < total_time
        push!(ts, t)
        push!(xs, beam_particles[monitor_point_index].x_vec_[1])
        push!(ys, beam_particles[monitor_point_index].x_vec_[2])
        applyInteractionWithNeighbours!(total_system, continuityAndPressureExtrapolation!)
        applySelfaction!(total_system, updateDensityAndPressure!; dt = dt)
        applyInteractionWithNeighbours!(total_system, momentumTotal!)
        applyReflection!(beam_system, calculateDeformationGradient!; dt = dt)
        applyReflection!(beam_system, momentumSolid!)
        applySelfaction!(total_system, calculateStepTime!)
        dt = collectmin(timeStep, total_system)
        dt = min(dt, reference_dt)
        applySelfaction!(total_system, accelerateAndMove!; dt = dt, t = t)
        createNeighbourIndexList!(total_system)
        if step % density_filter_step_interval == 0
            applyInteractionWithNeighbours!(total_system, densityFilterInteraction!)
            applySelfaction!(total_system, densityFilterSelfaction!)
        end
        if step % output_step_interval == 0
            saveVTP(vtp_writer, total_system, step, t)
        end
        step += 1
        t += dt
        percent = t / total_time
        update!(progress, 1 - percent; showvalues = [("progress percentage %", percent * 100), ("dt", dt)])
    end
    df = DataFrame(Time = ts, X = xs, Y = ys)
    CSV.write(csv_file_name, df)
    return nothing
end

function post()::Nothing
    splitParticlesType(
        vtp_writer;
        type_name_dict = Dict(1 => "fluid", 2 => "wall", 3 => "beam_movable", 4 => "beam_fixed"),
    )
    PyCall.@pyinclude "example/fsi_dam_break_onto_obstacle/fsi_dam_break_onto_obstacle.py"
    FSIDamBreakOntoObstaclePostProcess = py"FSIDamBreakOntoObstaclePostProcess"
    post_process = FSIDamBreakOntoObstaclePostProcess(key_word = "original")
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
