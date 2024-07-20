#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/19 16:52:54
  @ license: MIT
  @ description:
 =#

# see `example/hanging_beam/hanging_beam_continuum.jl`

# - [An accurate FSI-SPH modeling of challenging fluid-structure interaction problems in two and three dimensions](https://www.sciencedirect.com/science/article/abs/pii/S0029801820314608)
# - [An SPH-based FSI framework for phase-field modeling of brittle fracture under extreme hydrodynamic events](https://link.springer.com/article/10.1007/s00366-023-01857-0)
# - [An integrative SPH method for heat transfer problems involving fluid-structure interaction](https://link.springer.com/10.1007/s10409-022-22248-x)
# - [SPH modeling of fluid-structure interaction](http://link.springer.com/10.1007/s42241-018-0006-9)
# - [Sphinxsys-Tutorial: Schemes for solid dynamics](https://www.sphinxsys.org/html/theory.html#schemes-for-solid-dynamics)
# - [知乎-连续介质力学简介](https://zhuanlan.zhihu.com/p/595941534)

#=
    Large Deformation Constitutive Model - LDCM, 
    or say Large Deformation Continuum Mechanics - LDCM, 
    for this file containes a series functions that should be used in sequence,
    and have relations with each other.
    LDCM is added to function name
=#

@inline function libLDCMInitializeContinuum!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    p::ParticleType;
    smooth_kernel::SmoothKernel{Dimension},
    real_eps::Float64 = eps(),
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.rho_0_ = p.rho_
    p.volume_0_ = p.mass_ / p.rho_
    p.initial_neighbour_index_list_ = deepcopy(p.neighbour_index_list_)
    reset!(p.kernel_gradient_0_vec_list_)
    reset!(p.corrected_kernel_gradient_0_vec_list_)
    @inbounds @simd for i in eachindex(p.initial_neighbour_index_list_)
        q_index = p.initial_neighbour_index_list_[i]
        q = particle_system[q_index]
        kernel_gradient = kernelGradient(p.neighbour_distance_list_[i], smooth_kernel)
        libGatherCorrectedMatrix!(
            p,
            q,
            p.neighbour_position_list_[i],
            p.neighbour_distance_list_[i];
            kernel_gradient = kernel_gradient,
        )
        push!(
            p.kernel_gradient_0_vec_list_,
            kernel_gradient / p.neighbour_distance_list_[i] * p.neighbour_position_list_[i],
        )
    end
    libGenerateCorrectedMatrix!(p; real_eps = real_eps)
    @inbounds @simd for i in eachindex(p.initial_neighbour_index_list_)
        push!(p.corrected_kernel_gradient_0_vec_list_, dot(p.kernel_gradient_0_vec_list_[i], p.corrected_mat_))
    end
    return nothing
end

@inline function libLDCMDeformationGradientEvolution!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    p::ParticleType;
    dt::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @inbounds @simd for i in eachindex(p.initial_neighbour_index_list_)
        q = particle_system[p.initial_neighbour_index_list_[i]]
        p.deformation_gradient_mat_ .+=
            q.volume_0_ * dt * dyad((q.v_vec_ .- p.v_vec_), p.corrected_kernel_gradient_0_vec_list_[i])
    end
    p.rho_ = p.rho_0_ / det(p.deformation_gradient_mat_)
    return nothing
end

@inline function libLDCMGreenLagrangeStrain!(
    p::ParticleType;
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    #= 
    p.green_lagrange_strain_mat_ .= p.deformation_gradient_mat_' * p.deformation_gradient_mat_
    @inbounds @simd for i in 1:Dimension
        p.green_lagrange_strain_mat_[i, i] -= 1.0
    end
    p.green_lagrange_strain_mat_ .*= 0.5
    =#
    # * above code is slower than below
    p.green_lagrange_strain_mat_ .=
        0.5 * (p.deformation_gradient_mat_' * p.deformation_gradient_mat_ .- ConstMatI{Dimension}())
    return nothing
end

@inline function libLDCMLinearElasticityStress!(
    p::ParticleType;
    lambda::Float64 = 0.0,
    mu::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    # more complicated constitutive model can be implemented here, instead of linear elasticity
    #=
    lambda_strain_trace = lambda * trace(p.green_lagrange_strain_mat_)
    p.piola_kirchhoff_2nd_stress_mat_ .= 2 * mu * p.green_lagrange_strain_mat_
    @inbounds @simd for i in 1:Dimension
        p.piola_kirchhoff_2nd_stress_mat_[i, i] += lambda_strain_trace
    end
    =#
    # * above code is slower than below
    p.piola_kirchhoff_2nd_stress_mat_ .=
        lambda * trace(p.green_lagrange_strain_mat_) * ConstMatI{Dimension}() .+ 2 * mu * p.green_lagrange_strain_mat_
    return nothing
end

@inline function libLDCMPioalKirchhoffStress!(
    p::ParticleType,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    p.piola_kirchhoff_1st_stress_mat_ .= p.deformation_gradient_mat_ * p.piola_kirchhoff_2nd_stress_mat_
    p.corrected_piola_kirchhoff_1st_stress_mat_ .= p.piola_kirchhoff_1st_stress_mat_ * p.corrected_mat_
    return nothing
end

@inline function libLDCMContinuumMomentum!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    p::ParticleType;
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @inbounds @simd for i in eachindex(p.initial_neighbour_index_list_)
        q = particle_system[p.initial_neighbour_index_list_[i]]
        p.dv_vec_ .+=
            p.volume_0_ * q.volume_0_ / p.mass_ * dot(
                p.corrected_piola_kirchhoff_1st_stress_mat_ .+ q.corrected_piola_kirchhoff_1st_stress_mat_,
                p.kernel_gradient_0_vec_list_[i],
            )
    end
    return nothing
end

@inline function libLDCMContinuumArtificialStress!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    p::ParticleType;
    smooth_kernel::SmoothKernel{Dimension},
    alpha::Float64 = 0.1,
    beta::Float64 = 0.1,
    h::Float64 = 0.0,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    @inbounds @simd for i in eachindex(p.initial_neighbour_index_list_)
        q = particle_system[p.initial_neighbour_index_list_[i]]
        rpq = p.x_vec_ .- q.x_vec_
        r = norm(rpq)
        libArtificialViscosityForce!(
            p,
            q,
            rpq,
            r;
            kernel_gradient = kernelGradient(r, smooth_kernel),
            alpha = alpha,
            beta = beta,
            h = h,
        )
    end
    return nothing
end
