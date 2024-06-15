#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/16 01:22:24
  @ license: MIT
  @ description:
 =#

using EtherSPH
using SparseArrays
using CSV
using DataFrames

abstract type DemoFunction{Dimension} end

const DemoFunction1D = DemoFunction{1}

@kwdef struct SquareFunction <: DemoFunction1D
    a_::Float64 = 1.0
    b_::Float64 = 0.0
    c_::Float64 = 0.0
end

function Base.show(io::IO, f::SquareFunction)
    return print(io, "SquareFunction: $(f.a_)x^2 + $(f.b_)x + $(f.c_)")
end

function (f::SquareFunction)(x::Float64)::Float64
    return f.a_ * x * x + f.b_ * x + f.c_
end

function poisson(x::Float64, f::SquareFunction)::Float64
    return -2 * f.a_
end

abstract type PoissonProblem{Dimension} end

mutable struct PoissonProblem1D <: PoissonProblem{1}
    n_::Int64
    n_kernel_::Int64
    x_min_::Float64
    x_max_::Float64
    dx_::Float64
    f_::DemoFunction1D
    x_::Vector{Float64}
    u_::Vector{Float64}
    u_theory_::Vector{Float64}
end

function PoissonProblem1D(
    n::Int64,
    n_kernel::Int64,
    x_min::Float64,
    x_max::Float64,
    f::DemoFunction1D,
)::PoissonProblem1D
    dx = (x_max - x_min) / n
    x = [(i - 0.5) * dx + x_min for i in 1:n]
    u = zeros(n)
    u_theory = [f(x[i]) for i in 1:n]
    return PoissonProblem1D(n, n_kernel, x_min, x_max, dx, f, x, u, u_theory)
end

function solve!(problem::PoissonProblem1D, kernel::SmoothKernel)::Nothing
    n = problem.n_
    n_kernel = problem.n_kernel_
    x_min = problem.x_min_
    x_max = problem.x_max_
    dx = problem.dx_
    f = problem.f_
    x = problem.x_
    u = problem.u_
    A = spzeros(n, n)
    b = zeros(n)
    for i in 1:n
        for k in (-n_kernel):n_kernel
            j = i + k
            xi = x[i]
            xj = x_min + (j - 0.5) * dx
            rij = abs(xi - xj)
            wij = kernelValue(rij, kernel)
            dwij = kernelGradient(rij, kernel)
            m_devide_rho = dx
            if j in 1:n
                b[i] += -m_devide_rho * poisson(xj, f) * wij
                if i != j
                    A[i, j] += -m_devide_rho * dwij / rij * 2
                    A[i, i] += m_devide_rho * dwij / rij * 2
                end
            else
                uj = f(xj)
                b[i] += -m_devide_rho * poisson(xj, f) * wij
                b[i] += m_devide_rho * dwij / rij * uj * 2
                A[i, i] += m_devide_rho * dwij / rij * 2
            end
        end
    end
    u .= A \ b
    problem.u_ = u
    return nothing
end

const n = 100
const n_kernel = 3
const x_min = 1.0
const x_max = 8.0
const f = SquareFunction(3.0, 10.0, 2.0)
const problem = PoissonProblem1D(n, n_kernel, x_min, x_max, f)
const kernel = WendlandC4{1}(n_kernel * problem.dx_)

solve!(problem, kernel)

const save_path = "example/results/poisson_equation/poisson_equation_1d"
const save_filename = "poisson_equation_1d.csv"

assurePathExist(save_path)

df = DataFrame(x = problem.x_, u = problem.u_, u_theory = problem.u_theory_)
CSV.write(joinpath(save_path, save_filename), df)
