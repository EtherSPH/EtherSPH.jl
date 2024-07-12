#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/29 17:08:11
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const dim = 3
const dr = 0.114 / 20
const gap = dr
const h = 2.0 * dr

const kernel = WendlandC2{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const water_x_len = 0.114
const water_y_len = 0.114
const water_z_len = 0.228
const box_x_len = 0.42
const box_y_len = 0.44
const box_z_len = 0.228
const wall_width = 2 * dr

const rho_0 = 1000.0
const mass = rho_0 * dr^dim
const gravity = 9.8
const c = 10 * sqrt(2 * gravity * water_y_len)
const g = Vector3D(0.0, -gravity, 0.0)
const mu = 1e-3
const mu_wall = mu * 1000
const nu = mu / rho_0

const FLUID_TAG = 1
const WALL_TAG = 2

@inline function getPressureFromDensity(rho::Float64)::Float64
    return c^2 * (rho - rho_0)
end

@inline function getDensityFromPressure(p::Float64)::Float64
    return p / c^2 + rho_0
end

const dt = 0.1 * h / c
const t_end = 2.0
const output_dt = 100 * dt
const density_filter_dt = 20 * dt

@kwdef mutable struct Particle <: AbstractParticle3D
    # must have
    x_vec_::Vector3D = Vector3D(0.0, 0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # additonal properties
    p_::Float64 = 0.0
    drho_::Float64 = 0.0
    v_vec_::Vector3D = Vector3D(0.0, 0.0, 0.0)
    dv_vec_::Vector3D = Vector3D(0.0, 0.0, 0.0)
    c_::Float64 = c
    mu_::Float64 = mu
    gap_::Float64 = dr
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
    # * wall normal vector
    normal_vec_::Vector3D = Vector3D(0.0, 0.0, 0.0)
end

@inline function continuity!(p::Particle, q::Particle, rpq::Vector3D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = DW(r))
    end
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
    end
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector3D, r::Float64)::Nothing
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
    elseif p.type_ == FLUID_TAG && q.type_ == WALL_TAG
        dw = DW(r)
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.5 * h)
        EtherSPH.libCompulsiveForce!(p, q, rpq, r; h = h)
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

const x0 = 0.0
const y0 = 0.0
const z0 = 0.0

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1e3
    return nothing
end

particles = Particle[]
const fluid_column = Cuboid(Point3D(x0, y0, z0), Point3D(x0 + water_x_len, y0 + water_y_len, z0 + water_z_len))
fluid_particles = createParticles(Particle, gap, fluid_column; modify! = modifyFluid!)
append!(particles, fluid_particles)

# * ========== XOY plane ========== * #
@inline function modifyXOYBottomWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = Vector3D(0.0, 0.0, 1.0)
    return nothing
end
@inline function modifyXOYTopWall!(p::Particle)::Nothing
    modifyWall!(p)
    p.normal_vec_ = Vector3D(0.0, 0.0, -1.0)
    return nothing
end
const xoy_bottom_wall_column = Cuboid(Point3D(x0, y0, z0 - wall_width), Point3D(x0 + box_x_len, y0 + box_y_len, z0))
const xoy_top_wall_column =
    Cuboid(Point3D(x0, y0, z0 + box_z_len), Point3D(x0 + box_x_len, y0 + box_y_len, z0 + box_z_len + wall_width))
xoy_bottom_wall_particles = createParticles(Particle, gap, xoy_bottom_wall_column; modify! = modifyXOYBottomWall!)
xoy_top_wall_particles = createParticles(Particle, gap, xoy_top_wall_column; modify! = modifyXOYTopWall!)
append!(particles, vcat(xoy_bottom_wall_particles, xoy_top_wall_particles))

# * ========== YOZ plane ========== * #
@inline function modifyYOZLeftWall!(p::Particle)::Nothing
    modifyWall!(p)
    _, _, z = p.x_vec_
    if z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(1.0, 0.0, 0.0)
    elseif z < z0
        p.normal_vec_ = Vector3D(1.0, 0.0, 1.0) / sqrt(2)
    elseif z > z0 + box_z_len
        p.normal_vec_ = Vector3D(1.0, 0.0, -1.0) / sqrt(2)
    end
    return nothing
end
@inline function modifyYOZRightWall!(p::Particle)::Nothing
    modifyWall!(p)
    _, _, z = p.x_vec_
    if z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(-1.0, 0.0, 0.0)
    elseif z < z0
        p.normal_vec_ = Vector3D(-1.0, 0.0, 1.0) / sqrt(2)
    elseif z > z0 + box_z_len
        p.normal_vec_ = Vector3D(-1.0, 0.0, -1.0) / sqrt(2)
    end
    return nothing
end
const yoz_left_wall_column =
    Cuboid(Point3D(x0 - wall_width, y0, z0 - wall_width), Point3D(x0, y0 + box_y_len, z0 + box_z_len + wall_width))
const yoz_right_wall_column = Cuboid(
    Point3D(x0 + box_x_len, y0, z0 - wall_width),
    Point3D(x0 + box_x_len + wall_width, y0 + box_y_len, z0 + box_z_len + wall_width),
)
yoz_left_wall_particles = createParticles(Particle, gap, yoz_left_wall_column; modify! = modifyYOZLeftWall!)
yoz_right_wall_particles = createParticles(Particle, gap, yoz_right_wall_column; modify! = modifyYOZRightWall!)
append!(particles, vcat(yoz_left_wall_particles, yoz_right_wall_particles))

# * ========== ZOX plane ========== * #
@inline function modifyZOXFrontWall!(p::Particle)::Nothing
    modifyWall!(p)
    x, _, z = p.x_vec_
    if x0 < x < x0 + box_x_len && z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(0.0, 1.0, 0.0)
    elseif x < x0 && z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(1.0, 1.0, 0.0) / sqrt(2)
    elseif x > x0 + box_x_len && z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(-1.0, 1.0, 0.0) / sqrt(2)
    elseif x0 < x < x0 + box_x_len && z < z0
        p.normal_vec_ = Vector3D(0.0, 1.0, 1.0) / sqrt(2)
    elseif x0 < x < x0 + box_x_len && z > z0 + box_z_len
        p.normal_vec_ = Vector3D(0.0, 1.0, -1.0) / sqrt(2)
    elseif x < x0 && z < z0
        p.normal_vec_ = Vector3D(1.0, 1.0, 1.0) / sqrt(3)
    elseif x < x0 && z > z0 + box_z_len
        p.normal_vec_ = Vector3D(1.0, 1.0, -1.0) / sqrt(3)
    elseif x > x0 + box_x_len && z < z0
        p.normal_vec_ = Vector3D(-1.0, 1.0, 1.0) / sqrt(3)
    elseif x > x0 + box_x_len && z > z0 + box_z_len
        p.normal_vec_ = Vector3D(-1.0, 1.0, -1.0) / sqrt(3)
    end
    return nothing
end
@inline function modifyZOXBackWall!(p::Particle)::Nothing
    modifyWall!(p)
    x, _, z = p.x_vec_
    if x0 < x < x0 + box_x_len && z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(0.0, -1.0, 0.0)
    elseif x < x0 && z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(1.0, -1.0, 0.0) / sqrt(2)
    elseif x > x0 + box_x_len && z0 < z < z0 + box_z_len
        p.normal_vec_ = Vector3D(-1.0, -1.0, 0.0) / sqrt(2)
    elseif x0 < x < x0 + box_x_len && z < z0
        p.normal_vec_ = Vector3D(0.0, -1.0, 1.0) / sqrt(2)
    elseif x0 < x < x0 + box_x_len && z > z0 + box_z_len
        p.normal_vec_ = Vector3D(0.0, -1.0, -1.0) / sqrt(2)
    elseif x < x0 && z < z0
        p.normal_vec_ = Vector3D(1.0, -1.0, 1.0) / sqrt(3)
    elseif x < x0 && z > z0 + box_z_len
        p.normal_vec_ = Vector3D(1.0, -1.0, -1.0) / sqrt(3)
    elseif x > x0 + box_x_len && z < z0
        p.normal_vec_ = Vector3D(-1.0, -1.0, 1.0) / sqrt(3)
    elseif x > x0 + box_x_len && z > z0 + box_z_len
        p.normal_vec_ = Vector3D(-1.0, -1.0, -1.0) / sqrt(3)
    end
    return nothing
end
const zox_front_wall_column = Cuboid(
    Point3D(x0 - wall_width, y0 - wall_width, z0 - wall_width),
    Point3D(x0 + box_x_len + wall_width, y0, z0 + box_z_len + wall_width),
)
const zox_back_wall_column = Cuboid(
    Point3D(x0 - wall_width, y0 + box_y_len, z0),
    Point3D(x0 + box_x_len + wall_width, y0 + box_y_len + wall_width, z0 + box_z_len + wall_width),
)
zox_front_wall_particles = createParticles(Particle, gap, zox_front_wall_column; modify! = modifyZOXFrontWall!)
zox_back_wall_particles = createParticles(Particle, gap, zox_back_wall_column; modify! = modifyZOXBackWall!)
append!(particles, vcat(zox_front_wall_particles, zox_back_wall_particles))

const lower = Point3D(x0 - wall_width, y0 - wall_width, z0 - wall_width)
const upper = Point3D(x0 + box_x_len + wall_width, y0 + box_y_len + wall_width, z0 + box_z_len + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system.particles_, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline function getVelocity(p::Particle)::Vector3D
    if p.type_ == FLUID_TAG
        return p.v_vec_
    else
        return Vector3D(NaN, NaN, NaN)
    end
end
@inline function getNormal(p::Particle)::Vector3D
    if p.type_ == WALL_TAG
        return p.normal_vec_
    else
        return Vector3D(NaN, NaN, NaN)
    end
end
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
addVector!(vtp_writer, "Normal", getNormal)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "cruchaga_3d_"
vtp_writer.output_path_ = "example/results/cruchaga/cruchaga_3d"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createCellLinkList!(system)
    for step in ProgressBar(1:round(Int, t_end / dt))
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
    PyCall.@pyinclude "example/cruchaga/cruchaga_3d.py"
    CruchagaPostProcess = py"CruchagaPostProcess"
    post_process = CruchagaPostProcess()
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
