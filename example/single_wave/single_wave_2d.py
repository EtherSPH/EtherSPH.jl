'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/07/23 21:02:44
 # @ license: MIT
 # @ description:
 '''

import os
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
pv.set_jupyter_backend("static")
pv.start_xvfb()

working_directory: str = "example/single_wave/"
data_path: str = os.path.join(working_directory, "../results/single_wave/single_wave_2d")
file_list: str = os.listdir(data_path)
file_number: int = len(file_list)

def getStepPolyData(step: int) -> pv.PolyData:
    file_name: str = os.path.join(data_path, file_list[step])
    return pv.read(file_name)
    pass

def getInterpolatedPolyData() -> pv.PolyData:
    begin: float = 0.02
    end: float = 1 - 0.02
    n: int = 100
    interpolated_x: np.ndarray = np.linspace(begin, end, n)
    interpolated_y: np.ndarray = interpolated_x.copy()
    interpolated_points = np.column_stack((interpolated_x, interpolated_y, np.zeros(n)))
    interpolated_poly_data: pv.PolyData = pv.PolyData(interpolated_points)
    return interpolated_poly_data
    pass

def cmapPlot(step: int) -> None:
    poly_data: pv.PolyData = getStepPolyData(step)
    plotter: pv.Plotter = pv.Plotter()
    plotter.add_mesh(poly_data, scalars="U", cmap="jet", point_size=10)
    tmstep: int = int(poly_data.field_data["TMSTEP"][0])
    time_value: float = poly_data.field_data["TimeValue"][0]
    plotter.add_text(f"step: {tmstep}  time: {time_value}")
    plotter.camera_position = "xy"
    plotter.show()
    pass

def plot(step: int, i: int) -> None:
    poly_data: pv.PolyData = getStepPolyData(step)
    x, y, _ = poly_data.points.T
    u = poly_data.point_data["U"]
    plt.figure(figsize=(9, 4), facecolor="white")
    ax1 = plt.subplot(1, 2, 1)
    ax1.set_aspect("equal")
    ax1.scatter(x, y, c=u, cmap="jet", s=10)
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$y$")
    tmstep: int = int(poly_data.field_data["TMSTEP"][0])
    time_value: float = poly_data.field_data["TimeValue"][0]
    ax1.set_title(f"step: {tmstep}  time: {time_value:.2f}")
    interpolated_poly_data: pv.PolyData = getInterpolatedPolyData()
    interpolated_poly_data = interpolated_poly_data.interpolate(poly_data, radius=0.02, sharpness=2)
    i_x, _, _ = interpolated_poly_data.points.T
    i_u = interpolated_poly_data.point_data["U"]
    ax2 = plt.subplot(1, 2, 2)
    ax2.plot(i_x, i_u, label="$u$ along $x=y$")
    ax2.set_xlabel("$x$")
    ax2.set_ylabel("$u$")
    ax2.grid(True)
    ax2.legend(loc="lower right")
    ax2.set_title(f"step: {tmstep}  time: {time_value:.2f}")
    plt.savefig(os.path.join(working_directory, f"image/single_wave_2d_{i}.png"), bbox_inches="tight", dpi = 300)
    pass

plot_list: list = [0, 2, 4, 6, 8, 10, 12]
for i in range(len(plot_list)):
    step: int = int(np.floor(plot_list[i] / 100 * (file_number-1)))
    plot(step, i)
    pass

def initial(x: float, y: float) -> float:
    return np.sin(4 * np.pi * x) * np.sin(4 * np.pi * y)
    pass

in_x: np.ndarray = np.linspace(0, 1, 101)
in_y: np.ndarray = in_x.copy()
in_X, in_Y = np.meshgrid(in_x, in_y)
in_U = initial(in_X, in_Y)
plt.figure(figsize=(6, 6), facecolor="white")
ax1 = plt.subplot(1, 1, 1, projection="3d")
ax1.plot_surface(in_X, in_Y, in_U, cmap="jet")
ax1.set_xlabel("$x$")
ax1.set_ylabel("$y$")
ax1.set_zlabel("$u$")
ax1.set_title(r"Initial Condition $u=\sin(4\pi x)\sin(4\pi y)$")
plt.savefig(os.path.join(working_directory, "image/single_wave_2d_initial.png"), bbox_inches="tight", dpi=300)