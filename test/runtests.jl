#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 16:46:51
  @ license: MIT
  @ description:
 =#

using Test
using EtherSPH
using StaticArrays
using Einsum
using Parameters

@testset "EtherSPH" begin
    @testset "Math" begin
        include("MathTest.jl")
    end
    @testset "Geometry" begin
        include("GeometryTest.jl")
    end
    @testset "Function" begin
        include("FunctionTest.jl")
    end
    @testset "Core" begin
        include("CoreTest.jl")
    end
end
