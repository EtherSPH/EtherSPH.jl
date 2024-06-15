#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/14 18:04:03
  @ license: MIT
  @ description:
 =#

""" 
# all Library function shoule be written in following format:

```julia
@inline function someSelfaction!(p::ParticleType; parameters...)::Nothing
    # do something to p
    return nothing
end

@inline function someInteraction!(p::ParticleType, q::ParticleType, rpq::RealVector; parameters...)::Nothing
    # do something to p and q
    return nothing
end
```
"""

@inline function harmonicMean(a::Float64, b::Float64)::Float64
    return 2 * a * b / (a + b)
end

include("Continuity.jl")
include("PressureForce.jl")
include("ViscosityForce.jl")
include("AccelerateAndMove.jl")
include("UpdateDensity.jl")
include("DensityFilter.jl")
include("BoundaryForce.jl")
