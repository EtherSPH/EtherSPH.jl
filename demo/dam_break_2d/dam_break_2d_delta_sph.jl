#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com> | MaiZLnuaa <mai-zl@nuaa.edu.cn>
  @ date: 2024/09/12 16:31:00
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars

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
const c_0 = 10 * sqrt(2 * water_height * gravity)
const mu = 1e-3
const artificial_alpha = 0.1
const artificial_beta = 0.1

const dt = 0.1 * h / c_0
const t_end = 4.0
const output_dt = 100 * dt

const delta_plus_sph = 0.1

const FLUID_TAG = 1
const WALL_TAG = 2

@inline function eos(rho::Float64)::Float64
    return c_0 * c_0 * (rho - rho_0)
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
    c_::Float64 = c_0
    mu_::Float64 = mu
    gap_::Float64 = dr
    sum_kernel_weight_::Float64 = 0.0
    sum_kernel_weighted_rho_::Float64 = 0.0
    sum_kernel_weighted_p_::Float64 = 0.0
    # neighbour
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer()
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer()
    # ! for speed issue, neighbour kenel gradient is stored
    neighbour_kernel_gradient_list_::RealNumberContainer = RealNumberContainer()
    # δ+ SPH
    corrected_mat_::Matrix2D = Matrix2D(0.0, 0.0, 0.0, 0.0) # Lᵢ
    corrected_nabla_rho_vec_::Vector2D = Vector2D(0.0, 0.0) # ∇ρᵢ
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

# * ==================== δ SPH begin ====================

@inline function deltaSPHClearSelfaction!(p::Particle)::Nothing
    p.corrected_mat_ .= 0.0
    p.corrected_nabla_rho_vec_ .= 0.0
    return nothing
end

@inline function deltaSPHCorrectedMatrixInteraction!(
    p::Particle,
    q::Particle,
    rpq::Vector2D,
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing
    EtherSPH.libGatherCorrectedMatrix!(p, q, rpq, r; kernel_gradient = kernel_gradient)
    return nothing
end

@inline function deltaSPHCorrectedMatrixSelfaction!(p::Particle; real_eps::Float64 = eps())::Nothing
    EtherSPH.libGenerateCorrectedMatrix!(p; real_eps = real_eps)
    return nothing
end

@inline function deltaSPHCalculateDensityGradientInteraction!(
    p::Particle,
    q::Particle,
    rpq::Vector2D,
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing
    p.corrected_nabla_rho_vec_ .+=
        (q.rho_ - p.rho_) * kernel_gradient * q.mass_ / (q.rho_ * r) * dot(p.corrected_mat_, rpq)
    return nothing
end

@inline function deltaSPHCorrectionInteraction!(
    p::Particle,
    q::Particle,
    rpq::Vector2D,
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing
    p.drho_ -=
        delta_plus_sph * h * c_0 * q.mass_ * kernel_gradient / (q.rho_ * r) *
        (2 * (q.rho_ - p.rho_) + dot(p.corrected_nabla_rho_vec_ .+ q.corrected_nabla_rho_vec_, rpq))
    return nothing
end

# * ==================== δ SPH end   ====================

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
        deltaSPHClearSelfaction!(p)
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
        deltaSPHCorrectedMatrixInteraction!(p, q, rpq, r; kernel_gradient = kernel_gradient)
        return nothing
    end
    return nothing
end

@inline function deltaSPHStep3!(p::Particle)::Nothing
    if p.type_ == FLUID_TAG
        deltaSPHCorrectedMatrixSelfaction!(p)
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
        deltaSPHCalculateDensityGradientInteraction!(p, q, rpq, r; kernel_gradient = kernel_gradient)
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
        deltaSPHCorrectionInteraction!(p, q, rpq, r; kernel_gradient = kernel_gradient)
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
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = kernel_gradient, h = 0.5 * h)
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
        EtherSPH.libTraditionalViscosityForce!(p, q, rpq, r; kernel_gradient = kernel_gradient, h = 0.5 * h)
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

function createRectangleParticles(
    ParticleType::DataType,
    x0::Float64,
    y0::Float64,
    width::Float64,
    height::Float64,
    reference_dr::Float64;
    modifyOnParticle!::Function = (p) -> nothing,
)::Vector{ParticleType}
    particles = Vector{ParticleType}()                 # create particles, first an empty array
    n_along_x = Int64(width / reference_dr |> round)
    n_along_y = Int64(height / reference_dr |> round)
    dx = width / n_along_x
    dy = height / n_along_y
    for i in 1:n_along_x, j in 1:n_along_y
        particle = ParticleType()
        x = x0 + (i - 0.5) * dx # calculate the x position, in the center of each square
        y = y0 + (j - 0.5) * dy # calculate the y position, in the center of each square
        particle.x_vec_ = Vector2D(x, y)
        modifyOnParticle!(particle) # modify the particle as you wish
        particle.mass_ = particle.rho_ * dx * dy # calculate the mass of the particle
        push!(particles, particle) # add the particle to the array
    end
    return particles
end

const x0 = 0.0
const y0 = 0.0

particles = Particle[]

fluid_particles = createRectangleParticles(Particle, x0, y0, water_width, water_height, dr)
append!(particles, fluid_particles)

@inline function initialWall!(p::Particle)::Nothing
    p.type_ = WALL_TAG
    p.mu_ *= 1000
    return nothing
end

bottom_wall_particles = createRectangleParticles(
    Particle,
    x0 - wall_width,
    y0 - wall_width,
    box_width + 2 * wall_width,
    wall_width,
    dr;
    modifyOnParticle! = initialWall!,
)
append!(particles, bottom_wall_particles)

left_wall_particles = createRectangleParticles(
    Particle,
    x0 - wall_width,
    y0,
    wall_width,
    box_height,
    dr;
    modifyOnParticle! = initialWall!,
)
append!(particles, left_wall_particles)

right_wall_particles =
    createRectangleParticles(Particle, x0 + box_width, y0, wall_width, box_height, dr; modifyOnParticle! = initialWall!)
append!(particles, right_wall_particles)

lower = Vector2D(x0 - wall_width, y0 - wall_width)
upper = Vector2D(x0 + box_width + wall_width, y0 + box_height)
system = ParticleSystem(Particle, h, lower, upper)
append!(system, particles)

vtp_writer = VTPWriter()
@inline getPressure(p::Particle)::Float64 = p.p_
@inline getVelocity(p::Particle)::Vector2D = p.v_vec_
addScalar!(vtp_writer, "Pressure", getPressure)
addVector!(vtp_writer, "Velocity", getVelocity)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "dam_break_2d_delta_sph_"
vtp_writer.output_path_ = "demo/results/dam_break_2d_delta_sph"

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
