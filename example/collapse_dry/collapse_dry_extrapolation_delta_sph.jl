#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com> | MaiZLnuaa <mai-zl@nuaa.edu.cn>
  @ date: 2024/09/13 17:13:35
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const dim = 2
const dr = 0.02
const h = 3 * dr
const gap = dr

const kernel = CubicSpline{dim}(h)

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end
@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const water_width = 1.0
const water_height = 2.0
const box_width = 4.0
const box_height = 3.0
const wall_width = h

const rho_0 = 1e3
const mass = rho_0 * dr^dim
const gravity = 9.8
const g = Vector2D(0.0, -gravity)
const c = 10 * sqrt(2 * water_height * gravity)
const mu = 1e-3
const delta_sph_coefficient = 0.1
const artificial_alpha = 0.1
const artificial_beta = 0.1

const dt = 0.1 * h / c
const t_end = 4.0
const output_dt = 100 * dt

const FLUID_TAG = 1
const WALL_TAG = 2

@inline function eos(rho::Float64)::Float64
    return c * c * (rho - rho_0)
end

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # additonal properties
    p_::Float64 = 0.0
    drho_::Float64 = 0.0
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    c_::Float64 = c
    mu_::Float64 = mu
    gap_::Float64 = dr
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
    sum_kernel_weighted_p_::Float64 = 0.0
    # * neighbour
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer()
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer()
    # ! for speed issue, neighbour kenel gradient is stored
    neighbour_kernel_gradient_list_::RealNumberContainer = RealNumberContainer()
    # * δ SPH
    corrected_mat_::Matrix2D = Matrix2D(0.0, 0.0, 0.0, 0.0) # Lᵢ
    corrected_density_gradient_vec_::Vector2D = Vector2D(0.0, 0.0) # ∇ρᵢ
end

# ! for speed issue, neighbour kenel gradient is stored
@inline function findNeighbourKernelGradientSelfaction!(p::Particle)::Nothing
    reset!(p.neighbour_kernel_gradient_list_)
    @simd for i_neighbour in eachindex(p.neighbour_index_list_)
        @inbounds kernel_gradient = DW(p.neighbour_distance_list_[i_neighbour])
        push!(p.neighbour_kernel_gradient_list_, kernel_gradient)
    end
    return nothing
end

@inline function continuity!(
    p::Particle,
    q::Particle,
    rpq::Vector2D,
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libTraditionalContinuity!(p, q, rpq, r; kernel_gradient = kernel_gradient)
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ == WALL_TAG
        EtherSPH.libBalancedContinuity!(p, q, rpq, r; kernel_gradient = kernel_gradient)
        return nothing
    end
    return nothing
end

@inline function deltaSPHStep1!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libDeltaSPHClear!(p)
        return nothing
    end
    return nothing
end

@inline function deltaSPHStep2!(
    p::Particle,
    q::Particle,
    rpq::Vector2D,
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libDeltaSPHGatherCorrectedMatrix!(p, q, rpq, r; kernel_gradient = kernel_gradient)
        return nothing
    end
    return nothing
end

@inline function deltaSPHStep3!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libDeltaSPHGenerateCorrectedMatrix!(p)
        return nothing
    end
    return nothing
end

@inline function deltaSPHStep4!(
    p::Particle,
    q::Particle,
    rpq::Vector2D,
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libDeltaSPHGatherDensityGradientVector!(p, q, rpq, r; kernel_gradient = kernel_gradient)
        return nothing
    end
    return nothing
end

@inline function deltaSPHStep5!(
    p::Particle,
    q::Particle,
    rpq::Vector2D,
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        EtherSPH.libDeltaSPHDiffusiveFilter!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = kernel_gradient,
            delta_sph_coefficient = delta_sph_coefficient,
            h = h,
            sound_speed = c,
        )
        return nothing
    end
    return nothing
end

@inline function updateDensityAndPressure!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libUpdateDensity!(p; dt = dt)
        p.p_ = eos(p.rho_)
        return nothing
    end
    return nothing
end

@inline function pressureExtrapolationInteraction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64;)::Nothing
    if p.type_ == WALL_TAG && q.type_ == FLUID_TAG
        w = W(r)
        EtherSPH.libPressureExtrapolationInteraction!(p, q, rpq, r; kernel_value = w, body_force_vec = g)
        return nothing
    end
    return nothing
end

@inline function pressureExtrapolationSelfaction!(p::Particle)::Nothing
    if p.type_ == WALL_TAG
        EtherSPH.libPressureExtrapolationSelfaction!(p;)
        return nothing
    end
    return nothing
end

@inline function momentum!(p::Particle, q::Particle, rpq::Vector2D, r::Float64; kernel_gradient::Float64 = 0.0)::Nothing
    if p.type_ == FLUID_TAG && q.type_ == FLUID_TAG
        w = W(r)
        rw = W(0.5 * p.gap_ + q.gap_)
        EtherSPH.libTraditionalPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = kernel_gradient,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        # EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = kernel_gradient, h = 0.5 * h)
        # * apply artificial viscosity only
        EtherSPH.libArtificialViscosityForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = kernel_gradient,
            h = 0.5 * h,
            alpha = artificial_alpha,
            beta = artificial_beta,
        )
        return nothing
    elseif p.type_ == FLUID_TAG && q.type_ == WALL_TAG
        w = W(r)
        rw = W(0.5 * p.gap_ + q.gap_)
        EtherSPH.libDensityWeightedPressureForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = kernel_gradient,
            kernel_value = w,
            reference_kernel_value = rw,
        )
        # EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = kernel_gradient, h = 0.5 * h)
        # * apply artificial viscosity only
        EtherSPH.libArtificialViscosityForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = kernel_gradient,
            h = 0.5 * h,
            alpha = artificial_alpha,
            beta = artificial_beta,
        )
        return nothing
    end
    return nothing
end

@inline function accelerateAndMove!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        EtherSPH.libAccelerateAndMove!(p; dt = dt, body_force_vec = g)
        return nothing
    end
    return nothing
end

const x0 = 0.0
const y0 = 0.0

@inline function modifyWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1000
    return nothing
end

const fluid_column = Rectangle(Vector2D(x0, y0), Vector2D(x0 + water_width, y0 + water_height))

const bottom_wall_column =
    Rectangle(Vector2D(x0 - wall_width, y0 - wall_width), Vector2D(x0 + box_width + wall_width, y0))

const left_wall_column = Rectangle(Vector2D(x0 - wall_width, y0), Vector2D(x0, y0 + box_height))

const right_wall_column =
    Rectangle(Vector2D(x0 + box_width, y0), Vector2D(x0 + box_width + wall_width, y0 + box_height))

particles = Particle[]
fluid_particles = createParticles(Particle, gap, fluid_column)
bottom_wall_particles = createParticles(Particle, gap, bottom_wall_column; modify! = modifyWall!)
left_wall_particles = createParticles(Particle, gap, left_wall_column; modify! = modifyWall!)
right_wall_particles = createParticles(Particle, gap, right_wall_column; modify! = modifyWall!)

append!(particles, fluid_particles)
append!(particles, bottom_wall_particles)
append!(particles, left_wall_particles)
append!(particles, right_wall_particles)

const lower = Vector2D(x0 - wall_width, y0 - wall_width)
const upper = Vector2D(x0 + box_width + wall_width, y0 + box_height)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "collapse_dry_extropolation_delta_sph_"
vtp_writer.output_path_ = "example/results/collapse_dry/collapse_dry_extrapolation_delta_sph"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createNeighbourIndexList!(system)
    applySelfaction!(system, findNeighbourKernelGradientSelfaction!)
    for step in ProgressBar(1:round(Int, t_end / dt))
        applyInteractionByKernelGradient!(system, continuity!)
        applySelfaction!(system, deltaSPHStep1!)
        applyInteractionByKernelGradient!(system, deltaSPHStep2!)
        applySelfaction!(system, deltaSPHStep3!)
        applyInteractionByKernelGradient!(system, deltaSPHStep4!)
        applyInteractionByKernelGradient!(system, deltaSPHStep5!)
        applySelfaction!(system, updateDensityAndPressure!) # update the density and pressure
        applyInteractionWithNeighbours!(system, pressureExtrapolationInteraction!)
        applySelfaction!(system, pressureExtrapolationSelfaction!)
        applyInteractionByKernelGradient!(system, momentum!) # calculate the momentum
        applySelfaction!(system, accelerateAndMove!) # accelerate and move
        createNeighbourIndexList!(system) # create neighbour index list
        applySelfaction!(system, findNeighbourKernelGradientSelfaction!)
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end

function post()::Nothing
    PyCall.@pyinclude "example/collapse_dry/collapse_dry.py"
    CollapseDryPostProcess = py"CollapseDryPostProcess"
    post_process = CollapseDryPostProcess(key_word = "extrapolation_delta_sph")
    post_process.viewPlot()
    post_process.referencePlot()
    return nothing
end
