#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/22 19:41:54
  @ license: MIT
  @ description:
 =#

using EtherSPH
using ProgressBars
using PyCall

const c_vec = Vector2D(1.0, 1.0)
const n_t = 2

const dim = 2
const dr = 0.01
const gap = dr
const h = 2.0 * dr

const kernel = CubicSpline{dim}(h)
const w0 = kernel.kernel_value_0_

@inline function W(r::Float64)::Float64
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)::Float64
    return kernelGradient(r, kernel)
end

const edge_length = 1.0
const x_0 = 0.0
const y_0 = 0.0

const v_max = norm(c_vec)
const c = 100 * v_max
const dt = 0.1 * h / c
const total_time = 2.0
const output_dt = 500 * dt
@info "total steps: $(round(Int, total_time / dt))"

@inline function initialize(x::Float64, y::Float64)::Float64
    return sin(2 * pi * x / edge_length * n_t) * sin(2 * pi * y / edge_length * n_t)
end

@kwdef mutable struct Particle <: AbstractParticle2D
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    mass_::Float64 = 1.0
    rho_::Float64 = 1.0
    type_::Int64 = 1
    # user defined
    u_::Float64 = 0.0
    du_::Float64 = 0.0
    sum_kernel_weight_::Float64 = 0.0
    # neighbour info
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer()
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer()
    # corrected mat
    corrected_mat_::Matrix2D = Matrix2D(0.0, 0.0, 0.0, 0.0)
end

@inline function initialize!(p::Particle)::Nothing
    p.u_ = initialize(p.x_vec_.x, p.x_vec_.y)
    return nothing
end

@inline function interaction!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    dw = DW(r)
    p.du_ += q.mass_ / (q.rho_ * r) * (p.u_ - q.u_) * dw * dot(rpq, c_vec)
    return nothing
end

@inline function selfaction!(p::Particle)::Nothing
    p.u_ += p.du_ * dt
    p.du_ = 0.0
    return nothing
end

const square = Rectangle(Vector2D(x_0, y_0), Vector2D(x_0 + edge_length, y_0 + edge_length))

particles = createParticles(Particle, gap, square; modify! = initialize!)

system = ParticleSystem(Particle, h, Vector2D(x_0, y_0), Vector2D(x_0 + edge_length, y_0 + edge_length))
append!(system, particles)

setPeriodicBoundaryAlongXY!(system)

vtp_writer = VTPWriter()
@inline getU(p::Particle) = p.u_
addScalar!(vtp_writer, "U", getU)
vtp_writer.step_digit_ = 4
vtp_writer.file_name_ = "single_wave_2d_"
vtp_writer.output_path_ = "example/results/single_wave/single_wave_2d"

function main()::Nothing
    t = 0.0
    assurePathExist(vtp_writer)
    saveVTP(vtp_writer, system, 0, t)
    createNeighbourIndexList!(system)
    for step in ProgressBar(1:round(Int, total_time / dt))
        applyInteractionWithNeighbours!(system, interaction!)
        applySelfaction!(system, selfaction!)
        if step % round(Int, output_dt / dt) == 0
            saveVTP(vtp_writer, system, step, t)
        end
        t += dt
    end
    return nothing
end

function post()::Nothing
    @pyinclude "example/single_wave/single_wave_2d.py"
    return nothing
end
