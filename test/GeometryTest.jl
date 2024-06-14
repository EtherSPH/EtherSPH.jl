#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 22:34:06
  @ license: MIT
  @ description:
 =#

x_2 = @Vec [1.0, 2.0]
y_2 = @Vec [3.0, 4.0]
box_2d = Box(x_2, y_2)
@test box_2d.range_ ≈ Vector2D(2.0, 2.0)
@test isInsideShape(Vector2D(1.5, 3.5), box_2d) == true
@test isInsideShape(Vector2D(0.5, 3.5), box_2d) == false
@test isInsideShape(Vector2D(1.5, 4.5), box_2d) == false
@test isInsideShape(Vector2D(0.5, 4.5), box_2d) == false

x_3 = @Vec [1.0, 2.0, 3.0]
y_3 = @Vec [4.0, 5.0, 6.0]
box_3d = Box(x_3, y_3)
@test box_3d.range_ ≈ Vector3D(3.0, 3.0, 3.0)
@test isInsideShape(Vector3D(2.5, 4.5, 5.5), box_3d) == true
@test isInsideShape(Vector3D(0.5, 4.5, 5.5), box_3d) == false
@test isInsideShape(Vector3D(2.5, 6.5, 5.5), box_3d) == false
@test isInsideShape(Vector3D(2.5, 4.5, 7.5), box_3d) == false
