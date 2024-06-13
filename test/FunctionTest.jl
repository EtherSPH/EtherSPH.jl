#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 21:23:11
  @ license: MIT
  @ description:
 =#

@testset "SmoothKernel" begin
    for smooth_type in [CubicSpline, Gaussian, WendlandC2, WendlandC4]
        kernel = smooth_type{2}(1.0)
        W, DW = asCallable(kernel)
        @test W(0.0) ≈ kernel.kernel_value_0_
        @test DW(0.0) ≈ 0.0
        @test W(1.0) ≈ 0.0
    end
end
