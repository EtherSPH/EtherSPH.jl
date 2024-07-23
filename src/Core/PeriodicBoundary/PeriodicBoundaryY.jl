#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/22 14:45:07
  @ license: MIT
  @ description:
 =#

@inline function setPeriodicBoundaryAlongY!(cell_link_list::CellLinkList2D)::Nothing
    nx, ny = size(cell_link_list.cells_)
    displacement = Vector2D(0.0, 0.0)
    displacement[2] = cell_link_list.calculation_domain_box_.range_[2]
    Threads.@threads for i in 1:nx
        for delta_i in -1:1
            i_ = i + delta_i
            if i_ < 1 || i_ > nx
                continue
            else
                addNeighbour!(
                    cell_link_list.cells_[i, 1],
                    CartesianIndex2D(i_, ny);
                    relative_position_displacement = displacement,
                )
                addNeighbour!(
                    cell_link_list.cells_[i, ny],
                    CartesianIndex2D(i_, 1);
                    relative_position_displacement = -displacement,
                )
            end
        end
    end
    return nothing
end

@inline function setPeriodicBoundaryAlongY!(
    particle_system::ParticleSystem2D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle2D}
    setPeriodicBoundaryAlongX!(particle_system.cell_link_list_)
    return nothing
end

@inline function setPeriodicBoundaryAlongY!(cell_link_list::CellLinkList3D)::Nothing
    nx, ny, nz = size(cell_link_list.cells_)
    displacement = Vector3D(0.0, 0.0, 0.0)
    displacement[2] = cell_link_list.calculation_domain_box_.range_[2]
    Threads.@threads for index in 1:(nx * nz)
        i = div(index - 1, nz) + 1
        k = mod(index - 1, nz) + 1
        for delta_i in -1:1, delta_k in -1:1
            i_ = i + delta_i
            k_ = k + delta_k
            if i_ < 1 || i_ > nx || k_ < 1 || k_ > nz
                continue
            else
                addNeighbour!(
                    cell_link_list.cells_[i, 1, k],
                    CartesianIndex3D(i_, ny, k_);
                    relative_position_displacement = displacement,
                )
                addNeighbour!(
                    cell_link_list.cells_[i, ny, k],
                    CartesianIndex3D(i_, 1, k_);
                    relative_position_displacement = -displacement,
                )
            end
        end
    end
    return nothing
end

@inline function setPeriodicBoundaryAlongY!(
    particle_system::ParticleSystem3D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle3D}
    setPeriodicBoundaryAlongY!(particle_system.cell_link_list_)
    return nothing
end
