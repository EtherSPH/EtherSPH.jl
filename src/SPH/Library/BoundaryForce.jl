#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/15 17:58:35
  @ license: MIT
  @ description:
 =#

"""
SPH MODELING OF TSUNAMI WAVES, http://www.worldscientific.com/doi/abs/10.1142/9789812790910_0003, Rogers & Dalrymple - 2008
"""
@inline function libCompulsiveForce!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    h::Float64,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    psi = abs(dot(rpq, q.normal_vec_))
    xi = sqrt(max(0.0, r^2 - psi^2))
    eta = psi / q.gap_
    if eta > 1.0 || xi > q.gap_
        return nothing
    end
    p_xi = abs(1.0 + cos(pi * xi / q.gap_)) * 0.5
    verticle_v = dot(p.v_vec_ .- q.v_vec_, q.normal_vec_)
    beta = verticle_v > 0 ? 0 : 1
    r_psi = (0.01 * p.c_ * p.c_ + beta * p.c_ * abs(verticle_v)) * abs(1 - eta) / (sqrt(eta) * h)
    p.dv_vec_ .+= r_psi * p_xi * q.normal_vec_
    return nothing
end
