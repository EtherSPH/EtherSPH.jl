#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/14 17:33:49
  @ license: MIT
  @ description:
 =#

@inline function clearAllCells!(
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    # * 1. clear all the cell's contained particles id
    Threads.@threads for cell in particle_system.cell_link_list_.cells_
        reset!(cell)
    end
    # * 2. reset to_be_removed_cell_
    reset!(particle_system.cell_link_list_.to_be_removed_cell_)
    return nothing
end

@inline function removeParticles!(
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    # * 4. sort and remove
    length_of_to_be_removed_particles = length(particle_system.cell_link_list_.to_be_removed_cell_)
    if length_of_to_be_removed_particles > 0
        @inbounds to_be_removed_particle_index_list =
            particle_system.cell_link_list_.to_be_removed_cell_.contained_particle_index_list_[1:length_of_to_be_removed_particles]
        sort!(to_be_removed_particle_index_list)
        deleteat!(particle_system, to_be_removed_particle_index_list)
    end
    return nothing
end

@inline function addParticlesToCells!(
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    # * 5. add particles to cells
    Threads.@threads for i in eachindex(particle_system)
        @inbounds cartesian_index =
            getPositionCartesianIndexFromCellLinkList(particle_system[i].x_vec_, particle_system.cell_link_list_)
        @inbounds push!(particle_system.cell_link_list_.cells_[cartesian_index], i)
    end
    return nothing
end

# * ========== create cell link list ========== * #

@inline function createCellLinkList!(
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    clearAllCells!(particle_system)
    # * 3. add particles to to_be_removed_cell_
    Threads.@threads for i in eachindex(particle_system)
        @inbounds if !isInsideCellLinkListStrictCalculationDomain(
            particle_system[i].x_vec_,
            particle_system.cell_link_list_,
        )
            push!(particle_system.cell_link_list_.to_be_removed_cell_, i)
        end
    end
    removeParticles!(particle_system)
    addParticlesToCells!(particle_system)
    return nothing
end

# * ========== create neighbour index list ========== * #

@inline function createNeighbourIndexList!(
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    clearAllCells!(particle_system)
    # * 3. add particles to to_be_removed_cell_
    Threads.@threads for i in eachindex(particle_system)
        # ! additional operation begin
        reset!(particle_system[i].neighbour_index_list_)
        reset!(particle_system[i].neighbour_position_list_)
        reset!(particle_system[i].neighbour_distance_list_)
        # ! additional operation end
        @inbounds if !isInsideCellLinkListStrictCalculationDomain(
            particle_system[i].x_vec_,
            particle_system.cell_link_list_,
        )
            push!(particle_system.cell_link_list_.to_be_removed_cell_, i)
        end
    end
    removeParticles!(particle_system)
    addParticlesToCells!(particle_system)
    Threads.@threads for i in eachindex(particle_system)
        p = particle_system.particles_[i]
        cell_cartesian_index = getPositionCartesianIndexFromCellLinkList(p.x_vec_, particle_system.cell_link_list_)
        cell = particle_system.cell_link_list_.cells_[cell_cartesian_index]
        @inbounds for i_neighbour_cell in 1:length(cell.neighbour_cell_cartesian_grid_index_list_)
            neighbour_cell = particle_system.cell_link_list_.cells_[neighbourCartesianIndex(cell, i_neighbour_cell)]
            for j in 1:length(neighbour_cell)
                @inbounds q_index = neighbour_cell[j]
                if q_index == 0
                    break
                elseif q_index == i
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
                # ! additional operation begin
                push!(p.neighbour_index_list_, q_index)
                push!(p.neighbour_position_list_, rpq)
                push!(p.neighbour_distance_list_, r)
                # ! additional operation end
            end
        end
    end
    return nothing
end

# * ========== directly update neighbours ========== * #

@inline function libDirectlyUpdateNeighbours!(
    p::ParticleType;
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    # NOTE: this function could be used in `apply...`
    @simd for j in eachindex(p.neighbour_index_list_)
        @inbounds q_index = p.neighbour_index_list_[j]
        @inbounds p.neighbour_position_list_[j] .= p.x_vec_ .- particle_system[q_index]
        @inbounds p.neighbour_distance_list_[j] = norm(p.neighbour_position_list_[j])
    end
    return nothing
end

@inline function directlyUpdateNeighbours!(
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    Threads.@threads for i in eachindex(particle_system)
        p = particle_system.particles_[i]
        @simd for j in eachindex(p.neighbour_index_list_)
            @inbounds q_index = p.neighbour_index_list_[j]
            @inbounds p.neighbour_position_list_[j] .= p.x_vec_ .- particle_system[q_index]
            @inbounds p.neighbour_distance_list_[j] = norm(p.neighbour_position_list_[j])
        end
    end
    return nothing
end
