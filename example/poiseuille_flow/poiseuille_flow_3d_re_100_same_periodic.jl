#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/04 20:09:21
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const x0 = 0.0
const y0 = 0.0
const z0 = 0.0

const dim = 3
const dr = 1e-2 / 20
const gap = dr
const h = 2 * dr

const kernel = WendlandC2{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const reynolds_number = 100
const pipe_diameter = 1e-2
const pipe_radius = pipe_diameter / 2
const pipe_length = 3 * pipe_diameter
const wall_width = 2 * dr
const reference_length = pipe_diameter

const rho_0 = 1e3
const mass = rho_0 * dr^dim
const mu_0 = 1e-3
const u_mean = reynolds_number * mu_0 / (rho_0 * pipe_diameter)
const u_max = 2 * u_mean
const fx = 4 * mu_0 * u_max / (pipe_radius^2)
const ax = fx / rho_0
const g = Vector3D(ax, 0.0, 0.0)
const c_0 = 10 * u_max
const p_0 = rho_0 * c_0^2 * 0.01
@info "u max: $u_max, c_0: $c_0"

@inline function getPressureFromDensity(rho::Float64)::Float64
    return p_0 + c_0 * c_0 * (rho - rho_0)
end

const dt = 0.1 * h / c_0
const total_time = 50.0
const output_dt = 200 * dt
const density_filter_dt = 20 * dt
@info "total step: $(round(Int, total_time / dt))"

const FLUID_TAG = 1
const WALL_TAG = 2

@kwdef mutable struct Particle <: AbstractParticle3D
    # must have
    x_vec_::Vector3D = Vector3D(0.0, 0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # additonal properties
    p_::Float64 = p_0
    drho_::Float64 = 0.0
    v_vec_::Vector3D = Vector3D(0.0, 0.0, 0.0)
    dv_vec_::Vector3D = Vector3D(0.0, 0.0, 0.0)
    c_::Float64 = c_0
    mu_::Float64 = mu_0
    gap_::Float64 = gap
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_value_::Float64 = 0.0
end

@inline function continuity!(p::Particle, q::Particle, rpq::Vector3D, r::Float64)::Nothing
    if p.type_ == WALL_TAG && q.type_ == WALL_TAG
        return nothing
    end
    EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = DW(r))
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    EtherSPH.libUpdateDensity!(p; dt = dt)
    p.p_ = getPressureFromDensity(p.rho_)
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector3D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG
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
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = h / 2)
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
        # ! periodic boundary condition here !
        if p.x_vec_[1] > x0 + pipe_length
            p.x_vec_[1] -= pipe_length
        elseif p.x_vec_[1] < x0
            p.x_vec_[1] += pipe_length
        end
    end
    return nothing
end

@inline function densityFilterInteraction!(p::Particle, q::Particle, rpq::Vector3D, r::Float64)::Nothing
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

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1
    return nothing
end

@inline function modify!(p::Particle)::Nothing
    r = sqrt((p.x_vec_[2] - y0)^2 + (p.x_vec_[3] - z0)^2)
    if r < pipe_radius
        modifyFluid!(p)
    else
        modifyWall!(p)
    end
    return nothing
end

const cylinder = CylinderX(Point3D(x0, y0, z0), pipe_radius + wall_width, pipe_length)

particles = createParticles(Particle, gap, cylinder; modify! = modify!)

const lower = Vector3D(x0, y0 - pipe_radius - wall_width, z0 - pipe_radius - wall_width)
const upper = Vector3D(x0 + pipe_length, y0 + pipe_radius + wall_width, z0 + pipe_radius + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system.particles_, particles)

@info "total particles: $(length(system.particles_))"

# ! set the periodic boundary along x direction
setPeriodicBoundaryAlongX!(system)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline function getVelocity(p::Particle)::Vector3D
    if p.type_ == FLUID_TAG
        return p.v_vec_
    else
        return Vector3D(NaN, NaN, NaN)
    end
end
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "poiseuille_flow_3d_re_100_same_periodic_"
vtp_writer.output_path_ = "example/results/poiseuille_flow/poiseuille_flow_3d_re_100_same_periodic"

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
    PyCall.@pyinclude "example/poiseuille_flow/poiseuille_flow_3d.py"
    PoiseuilleFlow3DPostProcess = py"PoiseuilleFlow3DPostProcess"
    post_process = PoiseuilleFlow3DPostProcess(reference_gap = gap)
    post_process.viewPlot()
    return nothing
end
