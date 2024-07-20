#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/19 17:04:53
  @ license: MIT
  @ description:
 =#

@inline function libGatherCorrectedMatrix!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.corrected_mat_ .-= q.mass_ / (q.rho_ * r) * kernel_gradient * dyad(rpq, rpq)
    return nothing
end

@inline function libGenerateCorrectedMatrix!(
    p::ParticleType;
    real_eps::Float64 = eps(),
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    if abs(det(p.corrected_mat_)) > real_eps
        p.corrected_mat_ .= inv(p.corrected_mat_)
        return nothing
    else
        p.corrected_mat_ .= 0.0
        @inbounds @simd for i in 1:Dimension
            p.corrected_mat_[i, i] = 1.0
        end
        return nothing
    end
    return nothing
end
