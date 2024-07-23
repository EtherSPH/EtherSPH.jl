#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/22 16:42:13
  @ license: MIT
  @ description:
 =#

@inline function setPeriodicBoundaryAlongZ!(cell_link_list::CellLinkList3D)::Nothing
    nx, ny, nz = size(cell_link_list.cells_)
    displacement = Vector3D(0.0, 0.0, 0.0)
    displacement[1] = cell_link_list.calculation_domain_box_.range_[3]
    Threads.@threads for index in 1:(nx * ny)
        i = div(index - 1, ny) + 1
        j = mod(index - 1, ny) + 1
        for delta_i in -1:1, delta_j in -1:1
            i_ = i + delta_i
            j_ = j + delta_j
            if i_ < 1 || i_ > nx || j_ < 1 || j_ > ny
                continue
            else
                addNeighbour!(
                    cell_link_list.cells_[i, j, 1],
                    CartesianIndex3D(i_, j_, nz);
                    relative_position_displacement = displacement,
                )
                addNeighbour!(
                    cell_link_list.cells_[i, j, nz],
                    CartesianIndex3D(i_, j_, 1);
                    relative_position_displacement = -displacement,
                )
            end
        end
    end
    return nothing
end

@inline function setPeriodicBoundaryAlongZ!(
    particle_system::ParticleSystem3D{ParticleType},
)::Nothing where {ParticleType <: AbstractParticle3D}
    setPeriodicBoundaryAlongZ!(particle_system.cell_link_list_)
    return nothing
end
