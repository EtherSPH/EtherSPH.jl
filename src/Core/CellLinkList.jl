#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 21:58:11
  @ license: MIT
  @ description:
 =#

mutable struct CellLinkList{Dimension}
    reference_gap_::Float64
    gaps_::RealVector{Dimension} # dx, dy, dz
    calculation_domain_box_::Box{Dimension}
    stirct_calcilation_domain_::Shape
    cartesian_range_::CartesianRange{Dimension}
    cells_::Array{Cell{Dimension}, Dimension}
    to_be_removed_cell_::Cell{Dimension}
end

const CellLinkList2D = CellLinkList{Dimension2}
const CellLinkList3D = CellLinkList{Dimension3}

Base.show(io::IO, cell_link_list::CellLinkList{Dimension}) where {Dimension} = print(
    io,
    "CellLinkList{",
    Dimension,
    "}:\n",
    "    reference gap: $(cell_link_list.reference_gap_)\n",
    "    gaps: $(cell_link_list.gaps_)\n",
    "    calculation domain box: $(cell_link_list.calculation_domain_box_)\n",
    "    strict calculation domain: $(cell_link_list.stirct_calcilation_domain_)\n",
    "    cartesian range: $(cell_link_list.cartesian_range_)\n",
    "    cells: $(typeof(cell_link_list.cells_)), $(size(cell_link_list.cells_))\n",
    "    to be removed cell: $(typeof(cell_link_list.to_be_removed_cell_))\n",
)

@inline function createCellNeighbours!(cells::Array{Cell2D, Dimension2}, cartesian_range::CartesianRange2D)::Nothing
    nx, ny = size(cells)
    Threads.@threads for index in 1:(nx * ny)
        i = div(index - 1, ny) + 1
        j = mod(index - 1, ny) + 1
        cells[i, j] = Cell(Dimension2)
        for delta_i in -1:1, delta_j in -1:1
            ii = i + delta_i
            jj = j + delta_j
            cartesian_index = CartesianIndex2D(ii, jj)
            if isInsideCartesianRange(cartesian_index, cartesian_range)
                addNeighbour!(cells[i, j], cartesian_index)
            else
                continue
            end
        end
    end
    return nothing
end

@inline function createCellNeighbours!(cells::Array{Cell3D, Dimension3}, cartesian_range::CartesianRange3D)::Nothing
    # ! warning: neighbour contains itself
    nx, ny, nz = size(cells)
    Threads.@threads for index in 1:(nx * ny * nz)
        i = div(index - 1, ny * nz) + 1
        j = div(mod(index - 1, ny * nz), nz) + 1
        k = mod(index - 1, nz) + 1
        cells[i, j, k] = Cell(Dimension3)
        for delta_i in -1:1, delta_j in -1:1, delta_k in -1:1
            ii = i + delta_i
            jj = j + delta_j
            kk = k + delta_k
            cartesian_index = CartesianIndex3D(ii, jj, kk)
            if isInsideCartesianRange(cartesian_index, cartesian_range)
                addNeighbour!(cells[i, j, k], cartesian_index)
            else
                continue
            end
        end
    end
    return nothing
end

@inline function CellLinkList(
    reference_gap::Float64,
    calculation_domain_box::Box{Dimension},
    stirct_calcilation_domain::Shape,
)::CellLinkList{Dimension} where {Dimension}
    dimension = Dimension
    gaps = Vector0(dimension)
    cartesian_index_list = zeros(Int64, dimension)
    @simd for dim in 1:dimension
        @inbounds calculation_range = calculation_domain_box.range_[dim]
        @inbounds cartesian_index_list[dim] = Int64(floor(calculation_range / reference_gap))
        @inbounds gaps[dim] = calculation_range / cartesian_index_list[dim]
    end
    cartesian_index = CartesianIndex(cartesian_index_list...)
    cartesian_range = CartesianRange(cartesian_index)
    cells = Array{Cell{Dimension}, Dimension}(undef, cartesian_index_list...)
    createCellNeighbours!(cells, cartesian_range)
    to_be_removed_cell = Cell(dimension)
    return CellLinkList(
        reference_gap,
        gaps,
        calculation_domain_box,
        stirct_calcilation_domain,
        cartesian_range,
        cells,
        to_be_removed_cell,
    )
end

@inline function CellLinkList(
    reference_gap::Float64,
    calculation_domain_box::Box{Dimension},
)::CellLinkList{Dimension} where {Dimension}
    return CellLinkList(reference_gap, calculation_domain_box, calculation_domain_box)
end

@inline function CellLinkList(
    reference_gap::Float64,
    lower::Point{Dimension},
    upper::Point{Dimension},
)::CellLinkList{Dimension} where {Dimension}
    calculation_domain_box = Box(lower, upper)
    return CellLinkList(reference_gap, calculation_domain_box)
end

@inline function getPositionInGap(x::Float64, gap::Float64)::Int64
    return max(1, ceil(Int64, x / gap))
end

@inline function getPositionCartesianIndexFromCellLinkList(
    position::Point2D,
    cell_link_list::CellLinkList2D,
)::CartesianIndex2D
    @inbounds i =
        getPositionInGap(position[1] - cell_link_list.calculation_domain_box_.lower_[1], cell_link_list.gaps_[1])
    @inbounds j =
        getPositionInGap(position[2] - cell_link_list.calculation_domain_box_.lower_[2], cell_link_list.gaps_[2])
    return CartesianIndex2D(i, j)
end

@inline function getPositionCartesianIndexFromCellLinkList(
    position::Point3D,
    cell_link_list::CellLinkList3D,
)::CartesianIndex3D
    @inbounds i =
        getPositionInGap(position[1] - cell_link_list.calculation_domain_box_.lower_[1], cell_link_list.gaps_[1])
    @inbounds j =
        getPositionInGap(position[2] - cell_link_list.calculation_domain_box_.lower_[2], cell_link_list.gaps_[2])
    @inbounds k =
        getPositionInGap(position[3] - cell_link_list.calculation_domain_box_.lower_[3], cell_link_list.gaps_[3])
    return CartesianIndex3D(i, j, k)
end

@inline function isInsideCellLinkListStrictCalculationDomain(
    position::Point{Dimension},
    cell_link_list::CellLinkList{Dimension},
)::Bool where {Dimension}
    return isInsideShape(position, cell_link_list.stirct_calcilation_domain_)
end
