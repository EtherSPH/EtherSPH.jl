#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 21:57:50
  @ license: MIT
  @ description:
 =#

mutable struct Cell{Dimension}
    contained_particle_index_list_::IndexContainer
    thread_lock_::Base.Threads.ReentrantLock
    neighbour_cell_cartesian_grid_index_list_::Vector{CartesianIndex{Dimension}}
    # for periodic boundary condition, the relative position should be fixed by a displacement vector
    # x - y + v
    neighbour_cell_relative_position_displacement_list_::Vector{RealVector{Dimension}}
end

@inline function Cell(dimension::Int64)::Cell{dimension}
    return Cell(
        IndexContainer(),
        Base.Threads.ReentrantLock(),
        Vector{CartesianIndex{dimension}}(),
        Vector{RealVector{dimension}}(),
    )
end

const Cell2D = Cell{Dimension2}
const Cell3D = Cell{Dimension3}

@inline dimension(::Cell{Dimension}) where {Dimension} = Dimension

@inline Base.length(cell::Cell{Dimension}) where {Dimension} = length(cell.contained_particle_index_list_)
@inline Base.size(cell::Cell{Dimension}) where {Dimension} = (length(cell),)
@inline Base.getindex(cell::Cell{Dimension}, index::Int64) where {Dimension} =
    getindex(cell.contained_particle_index_list_, index)
@inline Base.setindex!(cell::Cell{Dimension}, value::Int64, index::Int64) where {Dimension} =
    setindex!(cell.contained_particle_index_list_, value, index)
@inline function Base.show(io::IO, cell::Cell{Dimension}) where {Dimension}
    print(io, "Cell dimension: $Dimension\n")
    print(io, "    contained particle index list: ")
    for i in 1:length(cell)
        print(io, cell[i])
        if i < length(cell)
            print(io, ", ")
        end
    end
    print(io, "\n    neighbour cell index list: ")
    for i in 1:length(cell.neighbour_cell_cartesian_grid_index_list_)
        print(io, cell.neighbour_cell_cartesian_grid_index_list_[i])
        if i < length(cell.neighbour_cell_cartesian_grid_index_list_)
            print(io, ", ")
        end
    end
    print(io, "\n")
    return nothing
end

@inline function reset!(cell::Cell{Dimension}) where {Dimension}
    reset!(cell.contained_particle_index_list_)
    return nothing
end

@inline function Base.push!(cell::Cell{Dimension}, value::Int64) where {Dimension}
    Base.Threads.lock(cell.thread_lock_)
    try
        orderedPush!(cell.contained_particle_index_list_, value)
    finally
        Base.Threads.unlock(cell.thread_lock_)
    end
    return nothing
end

@inline function neighbourCartesianIndex(
    cell::Cell{Dimension},
    index::Int64,
)::CartesianIndex{Dimension} where {Dimension}
    @inbounds return cell.neighbour_cell_cartesian_grid_index_list_[index]
end

@inline function neighbourRelativePositionDisplacement(
    cell::Cell{Dimension},
    index::Int64,
)::RealVector{Dimension} where {Dimension}
    @inbounds return cell.neighbour_cell_relative_position_displacement_list_[index]
end

@inline function addNeighbour!(
    cell::Cell{Dimension},
    cartesian_index::CartesianIndex{Dimension};
    relative_position_displacement::RealVector{Dimension} = Vector0(Dimension),
)::Nothing where {Dimension}
    push!(cell.neighbour_cell_cartesian_grid_index_list_, cartesian_index)
    push!(cell.neighbour_cell_relative_position_displacement_list_, relative_position_displacement)
    return nothing
end
