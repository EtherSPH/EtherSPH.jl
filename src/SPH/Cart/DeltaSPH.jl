#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com> | MaiZLnuaa <mai-zl@nuaa.edu.cn>
  @ date: 2024/09/12 21:00:26
  @ license: MIT
  @ description:
 =#

# see `example/collapse_dry/collapse_dry_extrapolation_delta_sph.jl`

# - [Free-surface flows solved by means of SPH schemes with numerical diffusive terms](https://linkinghub.elsevier.com/retrieve/pii/S0010465509003506)
# - [Propagation of gravity waves through an SPH scheme with numerical diffusive terms](https://linkinghub.elsevier.com/retrieve/pii/S0010465510004911)
# - [Numerical diffusive terms in weakly-compressible SPH schemes](https://linkinghub.elsevier.com/retrieve/pii/S0010465512002342)
# - [δ-SPH model for simulating violent impact flows](https://linkinghub.elsevier.com/retrieve/pii/S0045782510003725)
# - [Enhancement of δ-SPH for ocean engineering applications through incorporation of a background mesh scheme](https://linkinghub.elsevier.com/retrieve/pii/S0141118720310671)

#=
    DeltaSPH - δ-SPH,
    such calculation is so expansive.
    I do not recommend to use this method in your simulation.
    But I still keep this method in this file.
=#

@inline function libDeltaSPHClear!(
    p::ParticleType;
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.corrected_mat_ .= 0.0
    p.corrected_density_gradient_vec_ .= 0.0
    return nothing
end

@inline function libDeltaSPHGatherCorrectedMatrix!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    libGatherCorrectedMatrix!(p, q, rpq, r; kernel_gradient = kernel_gradient)
    return nothing
end

@inline function libDeltaSPHGenerateCorrectedMatrix!(
    p::ParticleType;
    real_eps::Float64 = eps(),
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    libGenerateCorrectedMatrix!(p; real_eps = real_eps)
    return nothing
end

@inline function libDeltaSPHGatherDensityGradientVector!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.corrected_density_gradient_vec_ .+=
        (q.rho_ - p.rho_) * kernel_gradient * q.mass_ / (q.rho_ * r) * dot(p.corrected_mat_, rpq)
    return nothing
end

@inline function libDeltaSPHDiffusiveFilter!(
    p::ParticleType,
    q::ParticleType,
    rpq::RealVector{Dimension},
    r::Float64;
    kernel_gradient::Float64 = 0.0,
    delta_sph_coefficient::Float64 = 0.1, # δ
    h::Float64 = 0.0,
    sound_speed::Float64 = 0.0, # c₀
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.drho_ -=
        delta_sph_coefficient * h * sound_speed * q.mass_ * kernel_gradient / (q.rho_ * r) *
        (2 * (q.rho_ - p.rho_) + dot(p.corrected_density_gradient_vec_ .+ q.corrected_density_gradient_vec_, rpq))
    return nothing
end
