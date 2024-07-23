#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/22 17:30:02
  @ license: MIT
  @ description:
 =#

@inline function setPeriodicBoundaryAlongXY!(cell_link_list::CellLinkList2D)::Nothing
    setPeriodicBoundaryAlongX!(cell_link_list)
    setPeriodicBoundaryAlongY!(cell_link_list)
    # four corner needs to be repaired
    nx, ny = size(cell_link_list.cells_)
    displacement = cell_link_list.calculation_domain_box_.range_
    dx = displacement[1]
    dy = displacement[2]
    addNeighbour!(
        cell_link_list.cells_[1, 1],
        CartesianIndex2D(nx, ny);
        relative_position_displacement = Vector2D(dx, dy),
    )
    addNeighbour!(
        cell_link_list.cells_[nx, ny],
        CartesianIndex2D(1, 1);
        relative_position_displacement = Vector2D(-dx, -dy),
    )
    addNeighbour!(
        cell_link_list.cells_[nx, 1],
        CartesianIndex2D(1, ny);
        relative_position_displacement = Vector2D(-dx, dy),
    )
    addNeighbour!(
        cell_link_list.cells_[1, ny],
        CartesianIndex2D(nx, 1);
        relative_position_displacement = Vector2D(dx, -dy),
    )
    return nothing
end

@inline function setPeriodicBoundaryAlongXY!(
    particle_system::ParticleSystem2D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle2D}
    setPeriodicBoundaryAlongXY!(particle_system.cell_link_list_)
    return nothing
end
