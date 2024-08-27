#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 16:32:00
  @ license: MIT
  @ description:
 =#

module EtherSPH

const Dimension1 = 1
const Dimension2 = 2
const Dimension3 = 3

using StaticArrays
import LinearAlgebra # LinearAlgebra import, instead of using.
using FLoops
using Parameters
using Dates
using WriteVTK
using ProgressBars
using ProgressMeter
using PyCall
using ExportAll

include("Math/Math.jl")
include("Container/Container.jl")
include("Geometry/Geometry.jl")
include("Function/Function.jl")
include("Core/Core.jl")
include("Post/Post.jl")
include("SPH/SPH.jl")
include("RigidBody/RigidBody.jl")

@inline greet() = println("Welcome to EtherSPH.jl, let's break down to particles!")

ExportAll.@exportAll()

end # module EtherSPH
