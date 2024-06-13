#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/09 00:50:14
  @ license: MIT
  @ description:
 =#

@testset "Cell" begin
    cell = Cell(2)
    for j in 1:10
        push!(cell, j)
    end

    @test length(cell) == 10
    @test cell[7] == 7

    reset!(cell)

    @test length(cell) == 0

    addNeighbour!(cell, CartesianIndex2D(1, 2))
    @test neighbourCartesianIndex(cell, 1) == CartesianIndex(1, 2)
    @test neighbourRelativePositionDisplacement(cell, 1) == Vector2D(0.0, 0.0)
end

@testset "CellLinkList" begin
    reference_gap = 0.13
    cell_link_list_2d = CellLinkList(reference_gap, Vector2D(-0.12, -0.23), Vector2D(1.22, 0.77))

    @test cell_link_list_2d.gaps_[1] ≈ 0.13399999999999998
    @test cell_link_list_2d.gaps_[2] ≈ 0.14285714285714285

    @test size(cell_link_list_2d.cells_) == (10, 7)

    @test getPositionCartesianIndexFromCellLinkList(Vector2D(0.0, 0.0), cell_link_list_2d) == CartesianIndex2D(1, 2)

    @test isInsideCellLinkListStrictCalculationDomain(Vector2D(0.0, 0.0), cell_link_list_2d) == true
    @test isInsideCellLinkListStrictCalculationDomain(Vector2D(-0.121, 0.0), cell_link_list_2d) == false

    cell_link_list_3d = CellLinkList(reference_gap, Vector3D(-0.12, -0.23, -0.27), Vector3D(1.22, 0.77, 0.88))

    @test cell_link_list_3d.gaps_[3] ≈ 0.14375

    @test size(cell_link_list_3d.cells_, 3) == 8

    @test getPositionCartesianIndexFromCellLinkList(Vector3D(0.0, 0.0, 0.0), cell_link_list_3d) ==
          CartesianIndex3D(1, 2, 2)

    @test isInsideCellLinkListStrictCalculationDomain(Vector3D(0.0, 0.0, 0.0), cell_link_list_3d) == true
    @test isInsideCellLinkListStrictCalculationDomain(Vector3D(0.0, 0.0, 0.88001), cell_link_list_3d) == false
end

@testset "ParticleSystem and Apply" begin
    @kwdef mutable struct P2 <: AbstractParticle2D
        x_vec_::Vector2D = Vector2D(0.0, 0.0)
        v_vec_::Vector2D = Vector2D(0.0, 0.0)
        mass_::Float64 = 0.0
        rho_::Float64 = 0.0
        type_::Float64 = 0.0
        t_::Float64 = 1.0
    end
    reference_gap = 0.13
    ps = ParticleSystem(P2, reference_gap, Vector2D(-0.12, -0.23), Vector2D(1.22, 0.77))

    @test length(ps) == 0
    @test size(ps)[1] == 0

    append!(ps, [P2(t_ = Float64(i)) for i in 1:5])
    push!(ps, P2(t_ = 6.0))
    @test length(ps) == 6
    pop!(ps)
    @test length(ps) == 5
    resize!(ps, 10)
    @test length(ps) == 10
    deleteat!(ps, [1, 6, 8, 10])
    @test length(ps) == 6
    resize!(ps, 3)

    append!(ps, [P2(t_ = Float64(i)) for i in 10:20])

    @test eltype(ps) == P2
    @test typeof(ps) <: ParticleSystem2D{P2}

    for i in eachindex(ps)
        j = div(i, 3) + 1
        k = i - 3 * (j - 1) + 1
        ps[i].x_vec_[1] = (0.12 + 1.22) / 7 * k - 0.5
        ps[i].x_vec_[2] = (0.23 + 0.77) / 4 * j - 0.7
    end

    createCellLinkList!(ps)

    @test length(ps) == 8
    @test ps.cell_link_list_.cells_[1, 1].contained_particle_index_list_[1] == 1
    @test ps.cell_link_list_.cells_[1, 2].contained_particle_index_list_[1] == 3
    @test ps.cell_link_list_.cells_[1, 4].contained_particle_index_list_[1] == 5
    @test ps.cell_link_list_.cells_[1, 6].contained_particle_index_list_[1] == 7
    @test ps.cell_link_list_.cells_[2, 1].contained_particle_index_list_[1] == 2
    @test ps.cell_link_list_.cells_[2, 2].contained_particle_index_list_[1] == 4
    @test ps.cell_link_list_.cells_[2, 4].contained_particle_index_list_[1] == 6
    @test ps.cell_link_list_.cells_[2, 6].contained_particle_index_list_[1] == 8

    @inline function move!(p::P2)::Nothing
        p.x_vec_ .+= p.v_vec_
        return nothing
    end

    @inline function force!(p::P2, q::P2, rpq::Vector2D, r::Float64)::Nothing
        p.v_vec_ .+= (q.x_vec_ - p.x_vec_) / r + rpq
        return nothing
    end

    applyInteraction!(ps, force!)
    applySelfaction!(ps, move!)

    @test ps[1].x_vec_ ≈ Vector2D(-0.11714285714285716, -0.19999999999999996)
end
