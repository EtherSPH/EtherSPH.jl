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

    C = Matrix0(2)
    for i in 1:2, j in 1:2
        for k in 1:2
            C[i, j] += A[i, k] * B[k, j]
        end
    end
    @test C ≈ A * B

    @test trace(A) ≈ 5.0
    @test det(A) ≈ -2.0

    D = @MArray [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    E = @MArray [9.0 8.0 7.0; 6.0 5.0 4.0; 3.0 2.0 1.0]
    @test D + E ≈ Matrix3D([10.0 10.0 10.0; 10.0 10.0 10.0; 10.0 10.0 10.0])

    @test D * E ≈ Matrix3D([30.0 24.0 18.0; 84.0 69.0 54.0; 138.0 114.0 90.0])
    @test D * 2 ≈ Matrix3D([2.0 4.0 6.0; 8.0 10.0 12.0; 14.0 16.0 18.0])

    @test dot(D, E) ≈ 165.0

    @test trace(D) ≈ 15.0
    @test det(D) ≈ 0.0
end

@testset "TinyLinearAlgebra" begin
    x = @MArray [1.0, 2.0]
    y = @MArray [3.0, 4.0]
    @test dyad(x, y) ≈ x * y'

    a = @MArray [1.0, 2.0, 3.0]
    b = @MArray [4.0, 5.0, 6.0]
    @test dyad(a, b) ≈ a * b'
end
