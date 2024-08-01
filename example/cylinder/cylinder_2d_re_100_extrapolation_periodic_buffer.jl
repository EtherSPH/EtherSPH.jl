#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/30 21:14:40
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall
using DataFrames
using CSV

const reynolds_number = 100.0
const u_mean = 1.0
const u_max = 1.5 * u_mean

const dim = 2
const dr = 0.41 / 64
const gap = dr
const h = 3.0 * dr

const x0 = 0.0
const y0 = 0.0
const wall_width = 3 * dr
const buffer_length = 10 * dr
const pipe_length = 2.2 + buffer_length
const pipe_width = 0.41
const cylinder_x = 0.2
const cylinder_y = 0.2
const cylinder_r = 0.05
const reference_length = 2 * cylinder_r

const kernel = WendlandC2{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const rho_0 = 1.0
const mass = rho_0 * dr^dim
const mu_0 = rho_0 * u_mean * reference_length / reynolds_number
const nu_0 = mu_0 / rho_0
const fx = 12 * mu_0 * u_mean / pipe_width^2
const ax = fx / rho_0
const ax_nu = ax / nu_0
const g = Vector2D(ax, 0.0)
const c_0 = 30 * u_mean # u_max is not real global u_max, just the inlet flow's u_max
const p_0 = 0.005 * rho_0 * c_0^2

@inline function getUFromY(y::Float64)::Float64
    return 0.5 * ax_nu * y * (pipe_width - y)
end

@inline function getPressureFromDensity(rho::Float64)::Float64
    return p_0 + c_0 * c_0 * (rho - rho_0)
end

const dt = 0.1 * h / c_0
const total_time = 5.0
const output_dt = 100 * dt
const density_filter_dt = 20 * dt
@info "total steps: $(round(Int, total_time / dt))"

const FLUID_TAG = 1
const WALL_TAG = 2
const CYLINDER_TAG = 3

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # additonal properties
    p_::Float64 = p_0
    drho_::Float64 = 0.0
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    c_::Float64 = c_0
    mu_::Float64 = mu_0
    gap_::Float64 = gap
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
    sum_kernel_weighted_p_::Float64 = 0.0
    # neighbour info
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer()
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer()
    # post
    vorticity_::Float64 = 0.0
end

@inline function setTheoreticalVelocity!(p::Particle)::Nothing
    p.v_vec_[1] = getUFromY(p.x_vec_[2] - y0)
    return nothing
end

@inline function continuityAndPressureExtrapolation!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = DW(r))
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, CYLINDER_TAG)
        EtherSPH.libBalancedContinuity!(p, q, rpq, r; kernel_gradient = DW(r))
        return nothing
    elseif p.type_ in (WALL_TAG, CYLINDER_TAG) && q.type_ == FLUID_TAG
        EtherSPH.libPressureExtrapolationInteraction!(p, q, rpq, r; kernel_value = W(r), body_force_vec = g)
        return nothing
    end
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = getPressureFromDensity(p.rho_)
        return nothing
    elseif p.type_ in (WALL_TAG, CYLINDER_TAG)
        EtherSPH.libPressureExtrapolationSelfaction!(p)
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
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.0 * h)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ in (WALL_TAG, CYLINDER_TAG)
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libDensityWeightedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.0 * h)
        return nothing
    elseif p.type_ in CYLINDER_TAG && q.type_ == FLUID_TAG
        # ! just for calculating the force exerted on the cylinder
        w = W(r)
        dw = DW(r)
        rw = W(0.5 * (p.gap_ + q.gap_))
        EtherSPH.libDensityWeightedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = dw,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = dw, h = 0.0 * h)
        return nothing
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
        # ! periodic boundary condition here !
        if p.x_vec_[1] > x0 + pipe_length
            p.x_vec_[1] -= pipe_length
            p.v_vec_[1] = getUFromY(p.x_vec_[2] - y0)
        elseif p.x_vec_[1] < x0
            p.x_vec_[1] += pipe_length
        elseif x0 + pipe_length - buffer_length < p.x_vec_[1] < x0 + pipe_length
            xi = (p.x_vec_[1] - (x0 + pipe_length - buffer_length)) / buffer_length
            p.v_vec_ .*= 1.0 - xi
            p.v_vec_[1] += getUFromY(p.x_vec_[2] - y0) * xi
            p.rho_ = p.rho_ * (1.0 - xi) + rho_0 * xi
            p.p_ = getPressureFromDensity(p.rho_)
        end
    elseif p.type_ == CYLINDER_TAG
        p.dv_vec_ .= 0.0
    end
    return nothing
end

@inline function densityFilterInteraction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
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

@inline function calculateVorticity!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == FLUID_TAG
        p.vorticity_ += DW(r) * q.mass_ / (q.rho_ * r) * cross(rpq, p.v_vec_ .- q.v_vec_)
    end
    return nothing
end

@inline function clearVorticity!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        p.vorticity_ = 0.0
    end
    return nothing
end

@inline function modifyFluid!(p::Particle)::Nothing
    p.type_ = FLUID_TAG
    p.v_vec_[1] = getUFromY(p.x_vec_[2] - y0)
    return nothing
end

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 100000
    return nothing
end

@inline function modifyCylinder!(p::Particle)::Nothing
    p.type_ = CYLINDER_TAG
    p.mu_ *= 100000
    return nothing
end

const fluid_column = Rectangle(Vector2D(x0, y0), Vector2D(x0 + pipe_length, y0 + pipe_width))
const bottom_wall_column = Rectangle(Vector2D(x0, y0 - wall_width), Vector2D(x0 + pipe_length, y0))
const top_wall_column =
    Rectangle(Vector2D(x0, y0 + pipe_width), Vector2D(x0 + pipe_length, y0 + pipe_width + wall_width))
const cylinder = Ring(Point2D(cylinder_x, cylinder_y), cylinder_r - wall_width, cylinder_r)

particles = Particle[]
fluid_particles = createParticles(Particle, dr, fluid_column, modify! = modifyFluid!)
bottom_wall_particles = createParticles(Particle, dr, bottom_wall_column, modify! = modifyWall!)
top_wall_particles = createParticles(Particle, dr, top_wall_column, modify! = modifyWall!)
cylinder_particles = createParticles(Particle, dr, cylinder, modify! = modifyCylinder!)

remove_index = Int[]
for (i, p) in enumerate(fluid_particles)
    r = norm(p.x_vec_ .- cylinder.center_)
    if r < cylinder.outer_radius_
        push!(remove_index, i)
    end
end

deleteat!(fluid_particles, remove_index)

append!(particles, vcat(fluid_particles, bottom_wall_particles, top_wall_particles, cylinder_particles))

@inline function computeCylinderForce(cylinder_particles::Vector{Particle})::Vector2D
    force = Vector2D(0.0, 0.0)
    for p in cylinder_particles
        force .+= p.dv_vec_ .* p.mass_
    end
    return force
end

const lower = Vector2D(x0, y0 - wall_width)
const upper = Vector2D(x0 + pipe_length, y0 + pipe_width + wall_width)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

# ! set the periodic boundary along x direction
setPeriodicBoundaryAlongX!(system)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
@inline getVorticity(p::Particle)::Float64 = p.vorticity_
addScalar!(vtp_writer, "Vorticity", getVorticity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "cyliner_2d_re_100_extrapolation_periodic_buffer_"
vtp_writer.output_path_ = "example/results/cylinder/cylinder_2d_re_100_extrapolation_periodic_buffer"
csv_save_path = vtp_writer.output_path_

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createNeighbourIndexList!(system)
    n_steps = round(Int, total_time / dt)
    time_vector = zeros(n_steps)
    cl_vector = zeros(n_steps)
    cd_vector = zeros(n_steps)
    for step in ProgressBar(1:n_steps)
        applyInteractionWithNeighbours!(system, continuityAndPressureExtrapolation!)
        applySelfaction!(system, updateDensityAndPressure!)
        applyInteractionWithNeighbours!(system, momentum!)
        force = computeCylinderForce(cylinder_particles)
        cl = 2 * force[2] / (u_mean * u_mean * reference_length)
        cd = 2 * force[1] / (u_mean * u_mean * reference_length)
        cl_vector[step] = cl
        cd_vector[step] = cd
        time_vector[step] = t
        applySelfaction!(system, accelerateAndMove!)
        createNeighbourIndexList!(system)
        if step % round(Int, density_filter_dt / dt) == 0
            applyInteractionWithNeighbours!(system, densityFilterInteraction!)
            applySelfaction!(system, densityFilterSelfaction!)
        end
        if step % round(Int, output_dt / dt) == 0
            applyInteractionWithNeighbours!(system, calculateVorticity!)
            saveVTP(vtp_writer, system, step, t)
            applySelfaction!(system, clearVorticity!)
        end
        t += dt
    end
    df = DataFrame(time = time_vector, cl = cl_vector, cd = cd_vector)
    CSV.write(joinpath(csv_save_path, "cl_cd.csv"), df)
    return nothing
end

function post()::Nothing
    splitParticlesType(vtp_writer; type_name_dict = Dict(1 => "fluid", 2 => "wall", 3 => "cylinder"))
    PyCall.@pyinclude "example/cylinder/cylinder_2d.py"
    Cylinder2DPostProcess = py"Cylinder2DPostProcess"
    post_process = Cylinder2DPostProcess(key_word = "re_100_extrapolation_periodic_buffer")
    post_process.viewPlot()
    post_process.clcdPlot()
    return nothing
end
