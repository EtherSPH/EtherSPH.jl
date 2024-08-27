#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/08/18 19:13:28
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressMeter
using PyCall
using CSV
using DataFrames

const D = 0.1
const pipe_width = 4.1 * D
const buffer_length = 1 * D
const pipe_length = 11 * D + buffer_length
const cylinder_r = 0.5 * D
const cylinder_d = D
const plate_length = 3.5 * D
const plate_width = 0.2 * D
const cylinder_x = 2.0 * D
const cylinder_y = 2.0 * D
const reference_length = D

const x_0 = 0.0
const y_0 = 0.0

@kwdef struct CylinderWithPlate <: Shape
    x_::Float64 = x_0 + cylinder_x
    y_::Float64 = y_0 + cylinder_y
    r_::Float64 = cylinder_r
    l_::Float64 = plate_length
    w_::Float64 = plate_width
    plate_start_x_::Float64 = x_ + sqrt(r_^2 - (0.5 * w_)^2)
    plate_start_y_::Float64 = y_ - 0.5 * w_
end

const cylinder_with_plate = CylinderWithPlate()

@inline function isInsideShape(p::Vector2D, cp::CylinderWithPlate)::Bool
    r = norm(p .- Point2D(cp.x_, cp.y_))
    if r <= cp.r_
        return true
    elseif cp.plate_start_x_ <= p[1] <= cp.plate_start_x_ + cp.l_ &&
           cp.plate_start_y_ <= p[2] <= cp.plate_start_y_ + cp.w_
        return true
    else
        return false
    end
    return false
end

const dim = 2
const dr_f = D / 20
const dr_s = dr_f / 2
const gap_f = dr_f
const gap_s = dr_s
const h_f = 3 * dr_f
const h_s = 2.6 * dr_s
const wall_width = 3 * dr_f

const ker_f = WendlandC2{dim}(h_f)
const ker_s = WendlandC2{dim}(h_s)

@inline function WF(r::Float64)::Float64
    return kernelValue(r, ker_f)
end
@inline function WS(r::Float64)::Float64
    return kernelValue(r, ker_s)
end
@inline function DWF(r::Float64)::Float64
    return kernelGradient(r, ker_f)
end
@inline function DWS(r::Float64)::Float64
    return kernelGradient(r, ker_s)
end

const w0_f = WF(0.0)
const w0_s = WS(0.0)

const reynolds_number = 20.0
const density_ratio = 10.0
const rho_f = 1.0
const rho_s = rho_f * density_ratio
const u_mean = 1.0
const u_max = 1.5 * u_mean
const mu_f = rho_f * u_mean * reference_length / reynolds_number
const nu_f = mu_f / rho_f
const fx = 12 * mu_f * u_mean / pipe_width^2
const ax = fx / rho_f
const ax_nu = ax / nu_f # aₓ / ν
const g = Vector2D(ax, 0.0)
const c_f = 20 * u_mean
const p_0 = 0.005 * rho_f * c_f^2
const dimensionless_young_modulus = 1.4e3
const young_modulus = dimensionless_young_modulus * rho_f * u_mean^2
const poisson_ratio = 0.4
const shear_modulus = young_modulus / (2 * (1 + poisson_ratio))
const bulk_modulus = young_modulus / (3 * (1 - 2 * poisson_ratio))
const c_s = sqrt(shear_modulus / rho_s)
const lambda_s = bulk_modulus - 2 * shear_modulus / 3
const mu_s = shear_modulus
const artificial_alpha = 0.1
const artificial_beta = 0.1
@info "fluid sound speed = $c_f, solid sound speed = $c_s"

const dimensionless_time = 200.0
const total_time = dimensionless_time * D / u_mean
const start_time = 0.1 * total_time
@info "total time = $total_time"
const output_step_interval = 500
const density_filter_interval = 30
const reference_dt = min(0.1 * h_f / c_f, 0.1 * h_s / c_s) * 2
@info "total steps around $(round(total_time / reference_dt))"

@inline function getUFromY(y::Float64; t::Float64 = 0.0)::Float64
    y -= y_0
    if t <= start_time
        return 0.5 * ax_nu * y * (pipe_width - y) * 0.5 * (1 - cos(pi * t / start_time))
    else
        return 0.5 * ax_nu * y * (pipe_width - y)
    end
end
@inline function getPressureFromDensity(rho::Float64)::Float64
    return p_0 + c_f * c_f * (rho - rho_f)
end

const FLUID_TAG = 1
const WALL_TAG = 2
const MOVABLE_SOLID_TAG = 3
const FIXED_SOLID_TAG = 4

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
    vorticity_::Float64 = 0.0
    # shared properties
    p_::Float64 = p_0
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
    # calculate step time
    step_time_::Float64 = Inf64
end

@inline function calculateStepTime!(p::Particle)::Nothing
    if p.type_ in (MOVABLE_SOLID_TAG, FIXED_SOLID_TAG)
        dv_dt = norm(p.dv_vec_)
        if dv_dt > 0
            p.step_time_ = 0.6 * min(sqrt(h_s / dv_dt), h_s / (norm(p.v_vec_) + c_s))
        end
        return nothing
    elseif p.type_ == FLUID_TAG
        p.step_time_ = min(0.6 * h_f / (norm(p.v_vec_) + c_f), 0.25 * min(h_f / u_max, h_f * h_f / nu_f))
        return nothing
    end
    return nothing
end

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    p.rho_ = rho_f
    p.c_ = c_f
    p.mu_ = mu_f
    p.v_vec_[1] = getUFromY(p.x_vec_[2])
    p.p_ = p_0
    p.gap_ = gap_f
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.rho_ = rho_f
    p.c_ = c_f
    p.mu_ = mu_f * 100000
    p.v_vec_ = Vector2D(0.0, 0.0)
    p.p_ = p_0
    p.gap_ = gap_f
    return nothing
end

@inline function modifySolid!(p::Particle)::Nothing
    x = p.x_vec_[1] - cylinder_with_plate.x_
    y = p.x_vec_[2] - cylinder_with_plate.y_
    r = sqrt(x^2 + y^2)
    if r <= cylinder_with_plate.r_ && x <= cylinder_with_plate.plate_start_x_ + cylinder_with_plate.l_
        p.type_ = FIXED_SOLID_TAG
    else
        p.type_ = MOVABLE_SOLID_TAG
    end
    p.rho_ = rho_s
    p.c_ = c_s
    p.mu_ *= 100000
    p.v_vec_ = Vector2D(0.0, 0.0)
    p.p_ = p_0
    p.gap_ = gap_s
    return nothing
end

@inline function modifyParticle!(p::Particle)::Nothing
    x = p.x_vec_[1]
    y = p.x_vec_[2]
    if y < y_0 || y > y_0 + pipe_width
        modifyWall!(p)
    elseif isInsideShape(p.x_vec_, cylinder_with_plate)
        modifySolid!(p)
    else
        modifyFluid!(p)
    end
    return nothing
end

const particle_column =
    Rectangle(Point2D(x_0, y_0 - wall_width), Point2D(x_0 + pipe_length, y_0 + pipe_width + wall_width))
particles = createParticles(Particle, gap_f, particle_column; modify! = modifyParticle!)

solid_ids = findall(p -> p.type_ in (MOVABLE_SOLID_TAG, FIXED_SOLID_TAG), particles)
deleteat!(particles, solid_ids)

const solid_cylinder = Circle(Point2D(cylinder_with_plate.x_, cylinder_with_plate.y_), cylinder_with_plate.r_)
solid_particles = createParticles(Particle, gap_s, solid_cylinder; modify! = modifySolid!)
move_ids = findall(p -> p.type_ == MOVABLE_SOLID_TAG, solid_particles)
deleteat!(solid_particles, move_ids)
const solid_plate = Rectangle(
    Point2D(cylinder_with_plate.plate_start_x_, cylinder_with_plate.plate_start_y_),
    Point2D(
        cylinder_with_plate.plate_start_x_ + cylinder_with_plate.l_,
        cylinder_with_plate.plate_start_y_ + cylinder_with_plate.w_,
    ),
)
append!(solid_particles, createParticles(Particle, gap_s, solid_plate; modify! = modifySolid!))

total_system = ParticleSystem(
    Particle,
    h_f,
    Vector2D(x_0, y_0 - wall_width),
    Vector2D(x_0 + pipe_length, y_0 + pipe_width + wall_width),
)
append!(total_system, vcat(particles, solid_particles))
setPeriodicBoundaryAlongX!(total_system)

solid_system = ParticleSystem(
    Particle,
    h_s,
    Vector2D(cylinder_x - cylinder_r, cylinder_y - cylinder_r),
    Vector2D(cylinder_x + cylinder_r + plate_length, cylinder_y + cylinder_r),
)
append!(solid_system, solid_particles)

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

setStepDigit!(vtp_writer, 4)
setFileName!(vtp_writer, "fsi_cylinder_with_plate_case_1_")
setOutputPath!(vtp_writer, "example/results/fsi_cylinder_with_plate/fsi_cylinder_with_plate_case_1")
@info "$vtp_writer"

@inline function continuityAndPressureExtrapolation!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        dw = DWF(r)
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, MOVABLE_SOLID_TAG, FIXED_SOLID_TAG)
        dw = DWF(r)
        EtherSPH.libBalancedContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ in (WALL_TAG, MOVABLE_SOLID_TAG, FIXED_SOLID_TAG) && q.type_ == FLUID_TAG
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
    elseif p.type_ in (WALL_TAG, MOVABLE_SOLID_TAG, FIXED_SOLID_TAG)
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
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, MOVABLE_SOLID_TAG, FIXED_SOLID_TAG)
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
    elseif p.type_ == MOVABLE_SOLID_TAG && q.type_ == FLUID_TAG
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
    if p.type_ == MOVABLE_SOLID_TAG
        EtherSPH.libLDCMContinuumMomentum!(particle_system, p)
        EtherSPH.libLDCMContinuumArtificialStress!(
            particle_system,
            p;
            smooth_kernel = ker_s,
            alpha = artificial_alpha,
            beta = artificial_beta,
            h = h_s,
        )
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle; dt::Float64 = 0.0, t::Float64 = 0.0)::Nothing
    calculateStepTime!(p)
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
        if p.x_vec_[1] > x_0 + pipe_length
            p.x_vec_[1] -= pipe_length
            p.v_vec_[1] = getUFromY(p.x_vec_[2]; t = t)
            p.v_vec_[2] = 0.0
        elseif p.x_vec_[1] < x_0
            p.x_vec_[1] += pipe_length
        elseif x_0 + pipe_length - buffer_length < p.x_vec_[1] < x_0 + pipe_length
            xi = (p.x_vec_[1] - (x_0 + pipe_length - buffer_length)) / buffer_length
            p.v_vec_ .*= 1.0 - xi
            p.v_vec_[1] += getUFromY(p.x_vec_[2]; t = t) * xi
            p.rho_ = p.rho_ * (1.0 - xi) + rho_f * xi
            p.p_ = getPressureFromDensity(p.rho_)
        end
        return nothing
    elseif p.type_ == MOVABLE_SOLID_TAG
        p.dv_vec_copy_ .= p.dv_vec_
        EtherSPH.libAccelerateAndMove!(p; dt = dt)
        return nothing
    elseif p.type_ == FIXED_SOLID_TAG
        p.dv_vec_ .= 0.0
        return nothing
    end
    return nothing
end

@inline function calculateVorticity!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG
        p.vorticity_ += DWF(r) * q.mass_ / (q.rho_ * r) * cross(rpq, p.v_vec_ .- q.v_vec_)
    end
    return nothing
end

@inline function clearVorticity!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        p.vorticity_ = 0.0
    end
    return nothing
end

@inline getVorticity(p::Particle)::Float64 = p.vorticity_
addScalar!(vtp_writer, "Vorticity", getVorticity)

@inline function timeStep(p::Particle)::Float64
    return p.step_time_
end

@inline function findMonitorPointIndex()
    index = 1
    x = 0.0
    y = 0.0
    for i in eachindex(total_system.particles_)
        p = total_system.particles_[i]
        if p.type_ in (FLUID_TAG, FIXED_SOLID_TAG)
            continue
        else
            x_, y_ = p.x_vec_
            if x < x_ && abs(y_ - 2 * D) <= dr_s
                x = x_
                y = y_
                index = i
            end
        end
    end
    return index
end

const monitor_point_index = findMonitorPointIndex()
@info "monitor point index = $monitor_point_index"
const csv_save_path = vtp_writer.output_path_

function main()::Nothing
    t = 0.0
    dt = 0.0
    step = 0
    ts = Float64[]
    xs = Float64[]
    ys = Float64[]
    assurePathExist(vtp_writer)
    createNeighbourIndexList!(solid_system)
    applyReflection!(solid_system, initializeSolid!)
    createNeighbourIndexList!(total_system)
    saveVTP(vtp_writer, total_system, 0, t)
    progress = ProgressThresh(0.0; desc = "Computing...", showspeed = true)
    while t < total_time
        push!(ts, t)
        push!(xs, total_system.particles_[monitor_point_index].x_vec_[1])
        push!(ys, total_system.particles_[monitor_point_index].x_vec_[2])
        applyInteractionWithNeighbours!(total_system, continuityAndPressureExtrapolation!)
        applySelfaction!(total_system, updateDensityAndPressure!; dt = dt)
        applyInteractionWithNeighbours!(total_system, momentumTotal!)
        applyReflection!(solid_system, calculateDeformationGradient!; dt = dt)
        applyReflection!(solid_system, momentumSolid!)
        applySelfaction!(total_system, accelerateAndMove!; dt = dt, t = t)
        createNeighbourIndexList!(total_system)
        if step % density_filter_interval == 0
            applyInteractionWithNeighbours!(total_system, densityFilterInteraction!)
            applySelfaction!(total_system, densityFilterSelfaction!)
        end
        if step % output_step_interval == 0
            applyInteractionWithNeighbours!(total_system, calculateVorticity!)
            saveVTP(vtp_writer, total_system, step, t)
            applySelfaction!(total_system, clearVorticity!)
        end
        dt = collectmin(timeStep, total_system)
        dt = min(dt, reference_dt)
        step += 1
        t += dt
        percent = t / total_time
        update!(progress, 1 - percent; showvalues = [("progress percentage %", percent * 100), ("dt", dt)])
    end
    df = DataFrame(Time = ts, X = xs, Y = ys)
    CSV.write(joinpath(csv_save_path, "monitor_point.csv"), df)
    return nothing
end

function post()::Nothing
    splitParticlesType(
        vtp_writer;
        type_name_dict = Dict(1 => "fluid", 2 => "wall", 3 => "solid_movable", 4 => "solid_fixed"),
    )
    PyCall.@pyinclude "example/fsi_cylinder_with_plate/fsi_cylinder_with_plate.py"
    FSICylinderWithPlatePostProcess = py"FSICylinderWithPlatePostProcess"
    post_process = FSICylinderWithPlatePostProcess(key_word = "case_1")
    post_process.viewPlot()
    return nothing
end
