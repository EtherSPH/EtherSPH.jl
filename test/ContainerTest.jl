#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/14 20:24:02
  @ license: MIT
  @ description:
 =#

@testset "IndexContainer" begin
    idc = IndexContainer()
    for i in 1:10
        push!(idc, i)
    end
    @test length(idc) == 10
    @test idc[1] == 1
    @test idc[10] == 10
    @test idc[5] == 5
    reset!(idc)
    @test length(idc) == 0
end

@testset "CartesianIndex" begin
    @test CartesianIndex2D(3, 4) == CartesianIndex(3, 4)
    @test CartesianIndex3D(3, 4, 5) == CartesianIndex(3, 4, 5)

    cr2d = CartesianRange(CartesianIndex(3, 4))
    @test isInsideCartesianRange(CartesianIndex(2, 2), cr2d) == true

    cr3d = CartesianRange(CartesianIndex(3, 4, 5))
    @test isInsideCartesianRange(CartesianIndex(2, 3, 4), cr3d) == true
end

@testset "CartesianGridIndex" begin
    end_2d = CartesianIndex2D(3, 4)
    cri = CartesianRange(end_2d)
    @test begin
        passed = true
        for i in 1:3, j in 1:4
            passed = passed && isInsideCartesianRange(CartesianIndex2D(i, j), cri)
        end
        passed
    end

    end_3d = CartesianIndex3D(3, 4, 5)
    cri = CartesianRange(end_3d)
    @test begin
        passed = true
        for i in 1:3, j in 1:4, k in 1:5
            passed = passed && isInsideCartesianRange(CartesianIndex3D(i, j, k), cri)
        end
        passed
    end
end
