#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/11 22:19:25
  @ license: MIT
  @ description:
 =#

@inline function applySelfaction!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    selfactionFunction!::Function;
    parameters...,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    Threads.@threads for i in eachindex(particle_system)
        @inbounds selfactionFunction!(particle_system[i]; parameters...)
    end
    return nothing
end

@inline function applyInteractionFunction!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    interactionFunction!::Function,
    p::ParticleType,
    p_index::Int64;
    parameters...,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    cell_cartesian_index = getPositionCartesianIndexFromCellLinkList(p.x_vec_, particle_system.cell_link_list_)
    cell = particle_system.cell_link_list_.cells_[cell_cartesian_index]
    @inbounds for i_neighbour_cell in 1:length(cell.neighbour_cell_cartesian_grid_index_list_)
        neighbour_cell = particle_system.cell_link_list_.cells_[neighbourCartesianIndex(cell, i_neighbour_cell)]
        for j in 1:length(neighbour_cell)
            @inbounds q_index = neighbour_cell[j]
            if q_index == 0
                break
            elseif q_index == p_index
                continue
            end
            rpq = minusplus(
                p.x_vec_,
                particle_system[q_index].x_vec_,
                neighbourRelativePositionDisplacement(cell, i_neighbour_cell),
            )
            r = norm(rpq)
            if r > particle_system.cell_link_list_.reference_gap_
                continue
            end
            interactionFunction!(p, particle_system[q_index], rpq, r; parameters...)
        end
    end
    return nothing
end

@inline function applyInteraction!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    interactionFunction!::Function;
    parameters...,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    Threads.@threads for i in eachindex(particle_system)
        @inbounds applyInteractionFunction!(particle_system, interactionFunction!, particle_system[i], i; parameters...)
    end
    return nothing
end
