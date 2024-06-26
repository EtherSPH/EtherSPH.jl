#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/26 21:15:53
  @ license: MIT
  @ description:
 =#

@inline function setPeriodicBoundaryAlongX!(cell_link_list::CellLinkList2D)::Nothing
    nx, ny = size(cell_link_list.cells_)
    displacement = Vector2D(0.0, 0.0)
    displacement[1] = cell_link_list.calculation_domain_box_.range_[1]
    Threads.@threads for j in 1:ny
        addNeighbour!(
            cell_link_list.cells_[1, j],
            CartesianIndex2D(nx, j);
            relative_position_displacement = displacement,
        )
        addNeighbour!(
            cell_link_list.cells_[nx, j],
            CartesianIndex2D(1, j);
            relative_position_displacement = -displacement,
        )
        if j == 1
            addNeighbour!(
                cell_link_list.cells_[1, 1],
                CartesianIndex2D(nx, 2);
                relative_position_displacement = displacement,
            )
            addNeighbour!(
                cell_link_list.cells_[nx, 1],
                CartesianIndex2D(1, 2);
                relative_position_displacement = -displacement,
            )
        elseif j == ny
            addNeighbour!(
                cell_link_list.cells_[1, ny],
                CartesianIndex2D(nx, ny - 1);
                relative_position_displacement = displacement,
            )
            addNeighbour!(
                cell_link_list.cells_[nx, ny],
                CartesianIndex2D(1, ny - 1);
                relative_position_displacement = -displacement,
            )
        else
            for delta_j in [-1, 1]
                addNeighbour!(
                    cell_link_list.cells_[1, j],
                    CartesianIndex2D(nx, j + delta_j);
                    relative_position_displacement = displacement,
                )
                addNeighbour!(
                    cell_link_list.cells_[nx, j],
                    CartesianIndex2D(1, j + delta_j);
                    relative_position_displacement = -displacement,
                )
            end
        end
    end
    return nothing
end

@inline function setPeriodicBoundaryAlongX!(
    particle_system::ParticleSystem2D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle2D}
    setPeriodicBoundaryAlongX!(particle_system.cell_link_list_)
    return nothing
end
