#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/17 14:59:14
  @ license: MIT
  @ description:
 =#

using EtherSPH
using SparseArrays

abstract type DemoFunction{Dimension} end

const DemoFunction2D = DemoFunction{2}

@kwdef struct BinarySquareFunction <: DemoFunction2D
    x0_::Float64 = 0.0
    y0_::Float64 = 0.0
    x1_::Float64 = 1.0
    y1_::Float64 = 1.0
    k_::Float64 = 1.0
end

function Base.show(io::IO, f::BinarySquareFunction)
    return println(io, "BinarySquareFunction: $(f.k_)(x - $(f.x0_)) * (x - $(f.x1_)) * (y - $(f.y0_)) * (y - $(f.y1_))")
end

function (f::BinarySquareFunction)(x::Float64, y::Float64)::Float64
    return f.k_ * (x - f.x0_) * (x - f.x1_) * (y - f.y0_) * (y - f.y1_)
end

function (f::BinarySquareFunction)(v::Vector2D)::Float64
    return f(v[1], v[2])
end

function poisson(x::Float64, y::Float64, f::BinarySquareFunction)::Float64
    return -(2 * (y - f.y0_) * (y - f.y1_) + 2 * (x - f.x0_) * (x - f.x1_)) * f.k_
end

function poisson(v::Vector2D, f::BinarySquareFunction)::Float64
    return poisson(v[1], v[2], f)
end

const TO_SOLVE_TAG = 1
const BOUNDARY_TAG = 2

@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector0(2)
    mass_::Float64 = 0.0
    rho_::Float64 = 1.0
    type_::Int64 = TO_SOLVE_TAG
    # addition
    x_index_::Int64 = 0
    y_index_::Int64 = 0
    dict_A_::Dict{CartesianIndex{2}, Float64} = Dict{CartesianIndex{2}, Float64}()
    b_::Float64 = 0.0
    u_::Float64 = 0.0
    u_theory_::Float64 = 0.0
end

@inline function add!(dict::Dict{CartesianIndex{2}, Float64}, index::CartesianIndex{2}, value::Float64)
    if haskey(dict, index)
        dict[index] += value
    else
        dict[index] = value
    end
end

abstract type PoissonProblem{Dimension} end

@kwdef mutable struct PoissonProblem2D <: PoissonProblem{2}
    x0_::Float64 = 0.0
    y0_::Float64 = 0.0
    x1_::Float64 = 1.0
    y1_::Float64 = 1.0
    reference_gap_::Float64 = 0.01
    n_x_::Int64 = round(Int64, (x1_ - x0_) / reference_gap_)
    n_y_::Int64 = round(Int64, (y1_ - y0_) / reference_gap_)
    n_addition_::Int64 = 3
    edge_::Float64 = reference_gap_ * n_addition_
    f_::DemoFunction2D = BinarySquareFunction(x0_ = x0_, y0_ = y0_, x1_ = x1_, y1_ = y1_, k_ = 1.0)
end

const dim = 2
problem = PoissonProblem2D(x0_ = 1.0, y0_ = 0.0, x1_ = 5.0, y1_ = 4.0, reference_gap_ = 0.05, n_addition_ = 3)
const kernel = WendlandC4{dim}(problem.edge_)

@inline function W(r::Float64)
    return kernelValue(r, kernel)
end

@inline function DW(r::Float64)
    return kernelGradient(r, kernel)
end

lower = Vector2D(problem.x0_ - problem.edge_, problem.y0_ - problem.edge_)
upper = Vector2D(problem.x1_ + problem.edge_, problem.y1_ + problem.edge_)
system = ParticleSystem(Particle, problem.edge_, lower, upper)
problem_domain = Rectangle(lower, upper)

function modify!(p::Particle)::Nothing
    x = p.x_vec_[1]
    y = p.x_vec_[2]
    if x < problem.x0_ || x > problem.x1_ || y < problem.y0_ || y > problem.y1_
        p.type_ = BOUNDARY_TAG
    else
        p.type_ = TO_SOLVE_TAG
        x_index = ceil(Int, (x - problem.x0_) / problem.reference_gap_)
        y_index = ceil(Int, (y - problem.y0_) / problem.reference_gap_)
        p.x_index_ = x_index
        p.y_index_ = y_index
    end
    return nothing
end

particles = createParticles(Particle, problem.reference_gap_, problem_domain; modify! = modify!);
append!(system, particles);

@inline function globalIndex(p::Particle)::Int64
    return (p.y_index_ - 1) * problem.n_x_ + p.x_index_
end

@inline function assembleAandb!(p::Particle, q::Particle, rpq::Vector2D, r::Float64)::Nothing
    if p.type_ == TO_SOLVE_TAG && q.type_ == TO_SOLVE_TAG
        p_index = globalIndex(p)
        q_index = globalIndex(q)
        Aij = -2 * q.mass_ / q.rho_ * DW(r) / r
        Aii = -Aij
        bi = -q.mass_ / q.rho_ * W(r) * poisson(q.x_vec_, problem.f_)
        add!(p.dict_A_, CartesianIndex(p_index, p_index), Aii)
        add!(p.dict_A_, CartesianIndex(p_index, q_index), Aij)
        p.b_ += bi
    elseif p.type_ == TO_SOLVE_TAG && q.type_ == BOUNDARY_TAG
        p_index = globalIndex(p)
        q_index = globalIndex(q)
        Aij = -2 * q.mass_ / q.rho_ * DW(r) / r
        Aii = -Aij
        bi = -q.mass_ / q.rho_ * W(r) * poisson(q.x_vec_, problem.f_)
        bi -= Aij * problem.f_(q.x_vec_)
        add!(p.dict_A_, CartesianIndex(p_index, p_index), Aii)
        p.b_ += bi
    end
    return nothing
end

@inline function assembleb!(p::Particle)::Nothing
    if p.type_ == TO_SOLVE_TAG
        bi = -p.mass_ / p.rho_ * poisson(p.x_vec_, problem.f_) * W(0.0)
        p.b_ += bi
    end
    return nothing
end

createCellLinkList!(system)
applyInteraction!(system, assembleAandb!)
applySelfaction!(system, assembleb!)

A = spzeros(Float64, problem.n_x_ * problem.n_y_, problem.n_x_ * problem.n_y_)
b = zeros(Float64, problem.n_x_ * problem.n_y_)
for p in system.particles_
    if p.type_ == TO_SOLVE_TAG
        for key in keys(p.dict_A_)
            A[key] += p.dict_A_[key]
        end
        b[globalIndex(p)] = p.b_
    end
end

u = A \ b;

@inline function writeu(p::Particle)::Nothing
    if p.type_ == TO_SOLVE_TAG
        p.u_ = u[globalIndex(p)]
        p.u_theory_ = problem.f_(p.x_vec_)
    end
    return nothing
end

applySelfaction!(system, writeu)

remove_index_list = Int64[]
for index in eachindex(system)
    p = system.particles_[index]
    if p.type_ == BOUNDARY_TAG
        push!(remove_index_list, index)
    end
end
deleteat!(system.particles_, remove_index_list)

vtp_writer = VTPWriter()
@inline getSPH(p::Particle) = p.u_
@inline getTheory(p::Particle) = p.u_theory_
addScalar!(vtp_writer, "U", getSPH)
addScalar!(vtp_writer, "UTheory", getTheory)
vtp_writer.step_digit_ = 1
vtp_writer.file_name_ = "poisson_equation_2d"
vtp_writer.output_path_ = "example/results/poisson_equation/poisson_equation_2d"

assurePathExist(vtp_writer)
saveVTP(vtp_writer, system, 0, 0.0)
