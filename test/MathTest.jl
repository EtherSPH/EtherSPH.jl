#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 16:47:50
  @ license: MIT
  @ description:
 =#

@testset "Vector" begin
    x2 = @MArray [1.0, 2.0]
    y2 = @MArray [3.0, 4.0]
    @test x2 + y2 ≈ Vector2D(4.0, 6.0)

    x3 = @MArray [1.0, 2.0, 3.0]
    y3 = @MArray [4.0, 5.0, 6.0]
    @test x3 + y3 ≈ Vector3D(5.0, 7.0, 9.0)

    @test Vector0(2) ≈ Vector2D(0.0, 0.0)
    @test Vector0(3) ≈ Vector3D(0.0, 0.0, 0.0)

    @test VectorX(2) ≈ Vector2D(1.0, 0.0)
    @test VectorY(2) ≈ Vector2D(0.0, 1.0)
    @test VectorX(3) ≈ Vector3D(1.0, 0.0, 0.0)
    @test VectorY(3) ≈ Vector3D(0.0, 1.0, 0.0)
    @test VectorZ(3) ≈ Vector3D(0.0, 0.0, 1.0)

    @test dot(x2, y2) ≈ 11.0
    @test dot(x3, y3) ≈ 32.0

    @test norm(x2) ≈ sqrt(5.0)
    @test norm(x3) ≈ sqrt(14.0)

    @test normalize(x2) ≈ Vector2D(1.0 / sqrt(5.0), 2.0 / sqrt(5.0))
    @test normalize(x3) ≈ Vector3D(1.0 / sqrt(14.0), 2.0 / sqrt(14.0), 3.0 / sqrt(14.0))

    @test cross(VectorX(3), VectorY(3)) ≈ VectorZ(3)
end

@testset "Matrix" begin
    A = @MArray [1.0 2.0; 3.0 4.0]
    B = @MArray [5.0 6.0; 7.0 8.0]
    @test A + B ≈ Matrix2D([6.0 8.0; 10.0 12.0])

    @test A * B ≈ Matrix2D([19.0 22.0; 43.0 50.0])
    @test A * 2.0 ≈ Matrix2D([2.0 4.0; 6.0 8.0])

    @test dot(A, B) ≈ 70.0

    @einsum C[i, j] := A[i, k] * B[k, j]
    @test C ≈ A * B
end

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
