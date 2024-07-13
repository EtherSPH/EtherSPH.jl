#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/13 17:14:52
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const dim = 2
const dr = 2.0 / 100
const gap = dr
const h = 3 * dr
const kernel = WendlandC2{dim}(h)

const box_width = 4.0
const box_height = 2.0
const water_width = 2.0
const water_height = 1.0
const wall_width = 3 * dr
const rigid_body_width = 0.2
const rigid_body_height = 0.2

const x_0 = 0.0
const y_0 = 0.0
const rigid_body_center_x = x_0 + 0.3 * box_width
const rigid_body_center_y = y_0 + water_height + h * 0.0 + 0.5 * rigid_body_height

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const rho_fluid = 1000.0
const rho_rigid_body = 500.0
const mu_fluid = 1e-3
const gravity = 9.8
const g = Vector2D(0.0, -gravity)
const c_fluid = 10.0 * sqrt(2 * water_height * gravity)

const dt = 0.05 * h / c_fluid
const total_time = 15.0
const output_dt = 500 * dt
const density_filter_dt = 30 * dt
const center_correct_dt = 20 * dt

const FLUID_TAG = 1
const WALL_TAG = 2
const RIGID_BODY_TAG = 3

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c_fluid * c_fluid * (rho - rho_fluid)
end

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector0(dim)
    rho_::Float64 = rho_fluid
    mass_::Float64 = rho_fluid * dr^dim
    type_::Int64 = FLUID_TAG
    # usually should contain
    gap_::Float64 = dr
    v_vec_::Vector2D = Vector0(dim)
    dv_vec_::Vector2D = Vector0(dim)
    drho_::Float64 = 0.0
    # fluid group fields
    p_::Float64 = 0.0
    c_::Float64 = c_fluid
    mu_::Float64 = mu_fluid
    # sum kernel ...
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
    sum_kernel_weighted_p_::Float64 = 0.0
    # normal vector
    normal_vec_::Vector2D = Vector0(dim)
end

@inline function continuityAndPressureExtrapolation!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        dw = DW(r)
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, RIGID_BODY_TAG)
        dw = DW(r)
        EtherSPH.libBalancedContinuity!(p, q, rpq, r; kernel_gradient = dw)
        return nothing
    elseif p.type_ in (WALL_TAG, RIGID_BODY_TAG) && q.type_ == FLUID_TAG
        w = W(r)
        kernel_weight = w * q.mass_ / q.rho_
        p.sum_kernel_weight_ += kernel_weight
        p.sum_kernel_weighted_p_ += kernel_weight * (max(q.p_, 0.0) + q.rho_ * max(dot(g, rpq), 0.0))
        return nothing
    end
    return nothing
end

@inline function updateDensityAndPressureWithPressureExtrapolation!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
        return nothing
    elseif p.type_ in (WALL_TAG, RIGID_BODY_TAG)
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

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libTraditionalPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, RIGID_BODY_TAG)
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libBalancedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        return nothing
    elseif p.type_ == RIGID_BODY_TAG && q.type_ == FLUID_TAG
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libBalancedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        return nothing
    elseif p.type_ == RIGID_BODY_TAG && q.type_ == WALL_TAG
        dw = DW(r)
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        EtherSPH.libCompulsiveForce!(p, q, rpq, r; h = 0.5 * h)
        return nothing
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
    end
    return nothing
end

@inline function densityFilter!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libKernelAverageDensityFilterInteraction!(p, q, rpq, r; kernel_value = W(r))
    end
    return nothing
end

@inline function densityFilter!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libKernelAverageDensityFilterSelfaction!(p; kernel_value = kernel.kernel_value_0_)
    end
    return nothing
end

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    p.p_ = p.rho_ * gravity * (water_height - (p.x_vec_[2] - y_0))
    p.rho_ = p.p_ / c_fluid / c_fluid + rho_fluid
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    if (water_height - (p.x_vec_[2] - y_0)) > 0.0
        p.p_ = p.rho_ * gravity * (water_height - (p.x_vec_[2] - y_0))
        p.rho_ = p.p_ / c_fluid / c_fluid + rho_fluid
    end
    return nothing
end

@inline function modifyBottomWall!(p::Particle)::Nothing
    modifyWall!(p)
    if p.x_vec_[1] < x_0
        p.normal_vec_ = Vector2D(1.0, 1.0) ./ sqrt(2)
    elseif p.x_vec_[1] > x_0 + box_width
        p.normal_vec_ = Vector2D(-1.0, 1.0) ./ sqrt(2)
    else
        p.normal_vec_ = Vector2D(0.0, 1.0)
    end
    return nothing
end

@inline function modifyLeftWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = Vector2D(1.0, 0.0)
    return nothing
end

@inline function modifyRightWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = Vector2D(-1.0, 0.0)
    return nothing
end

@inline function modifyRigidBody!(p::Particle)::Nothing
    p.type_ = RIGID_BODY_TAG
    p.rho_ = rho_rigid_body
    return nothing
end

const fluid_column = Rectangle(Point2D(x_0, y_0), Point2D(x_0 + water_width, y_0 + water_height))
const botton_wall_column =
    Rectangle(Point2D(x_0 - wall_width, y_0 - wall_width), Point2D(x_0 + box_width + wall_width, y_0))
const left_wall_column = Rectangle(Point2D(x_0 - wall_width, y_0), Point2D(x_0, y_0 + box_height))
const right_wall_column =
    Rectangle(Point2D(x_0 + box_width, y_0), Point2D(x_0 + box_width + wall_width, y_0 + box_height))
const rigid_body_column = Rectangle(
    Point2D(rigid_body_center_x - 0.5 * rigid_body_width, rigid_body_center_y - 0.5 * rigid_body_height),
    Point2D(rigid_body_center_x + 0.5 * rigid_body_width, rigid_body_center_y + 0.5 * rigid_body_height),
)

particles = Particle[]

fluid_particles = createParticles(Particle, gap, fluid_column; modify! = modifyFluid!)
bottom_wall_particles = createParticles(Particle, gap, botton_wall_column; modify! = modifyBottomWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyLeftWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; modify! = modifyRightWall!)
rigid_body_particles = createParticles(Particle, gap, rigid_body_column; modify! = modifyRigidBody!)
append!(
    particles,
    vcat(fluid_particles, bottom_wall_particles, left_wall_particles, right_wall_particles, rigid_body_particles),
)

# rigid body particles can be reached either by system or by rigid body

rigid_body = SingleRigidBody2D(particles_ = rigid_body_particles)
computeMass!(rigid_body)
computeCenter!(rigid_body)
computeInertia!(rigid_body)

const lower = Vector2D(x_0 - wall_width, y_0 - wall_width)
const upper = Vector2D(x_0 + box_width + wall_width, y_0 + box_height + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline function getPressure(p::Particle)::Float64
    if p.type_ == FLUID_TAG
        return p.p_
    else
        return NaN
    end
end
@inline function getVelocity(p::Particle)::Vector2D
    if p.type_ == WALL_TAG
        return Vector2D(0.0, 0.0)
    else
        return p.v_vec_
    end
end
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "floating_object_2d_density_500_"
vtp_writer.output_path_ = "example/results/floating_object/floating_object_2d_density_500"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:round(Int, total_time / dt))
        applyInteraction!(system, continuityAndPressureExtrapolation!)
        applySelfaction!(system, updateDensityAndPressureWithPressureExtrapolation!)
        applyInteraction!(system, momentum!)
        applySelfaction!(system, accelerateAndMove!)
        motionRigidBody!(rigid_body; dt = dt, body_force_vec = g)
        createCellLinkList!(system)
        if step % round(Int, density_filter_dt / dt) == 0
            applyInteraction!(system, densityFilter!) # density filter
            applySelfaction!(system, densityFilter!) # density filter
        end
        if step % round(Int, center_correct_dt / dt) == 0
            computeCenter!(rigid_body)
        end
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end

function post()::Nothing
    PyCall.@pyinclude "example/floating_object/floating_object_2d.py"
    FloatingObjectPostProcess = py"FloatingObjectPostProcess"
    post_process = FloatingObjectPostProcess(key_word = "2d_density_500")
    post_process.viewPlot()
    return nothing
end
