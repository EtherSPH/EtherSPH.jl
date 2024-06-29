'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/27 22:14:58
 # @ license: MIT
 # @ description:
 '''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
pv.set_jupyter_backend("static")
pv.start_xvfb()

POST_NAME: str = "poiseuille_flow"

class PoiseuilleFlow2DPostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    
    def __init__(
        self,
        working_directory: str = f"example/{POST_NAME}",
        reference_gap: float = 0.01,
        key_word: str = "2d_re_100_same_periodic",
        reynolds_number: float = 100,
        pipe_width: float = 0.01,
        pipe_length: float = 0.04,
        nu: float = 1e-6,
        ax: float = 1,
        analytical_order: int = 10,
    ) -> None:
        self.working_directory: str = working_directory
        self.data_path: str = os.path.join(self.working_directory, f"../results/{POST_NAME}/{POST_NAME}_{key_word}")
        self.reference_gap: float = reference_gap
        self.h: float = 3 * self.reference_gap
        self.key_word: str = key_word
        self.reynolds_number: float = reynolds_number
        self.pipe_width: float = pipe_width
        self.pipe_length: float = pipe_length
        self.nu: float = nu
        self.ax: float = ax
        self.ax_nu: float = ax / nu
        self.analytical_order: int = analytical_order
        self.readData()
        self.makeInterpolatedGrid()
        pass
    
    def readData(self) -> None:
        self.file_list: list = os.listdir(self.data_path)
        self.view_step: int = len(self.file_list) - 1
        pass
    
    def readStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.file_list[step])
        return pv.read(file)
        pass
    
    def getTime(self, poly_data: pv.PolyData) -> float:
        return poly_data.field_data["TimeValue"][0]
        pass
    
    def getStep(self, poly_data: pv.PolyData) -> int:
        return poly_data.field_data["TMSTEP"][0]
        pass
    
    def fastplot(self, step: int, scalars: str = "Velocity", point_size: str = 5) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self, scalars: str = "Velocity", point_size: str = 5) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(self.view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_cmap.png"))
        pass
    
    def makeInterpolatedGrid(self) -> None:
        n: int = int(np.round(self.pipe_width / self.reference_gap)) + 1
        x: np.ndarray = np.zeros(n) + self.pipe_length / 4 * 3
        y: np.ndarray = np.linspace(0, self.pipe_width, n)
        z: np.ndarray = np.zeros(n)
        points: np.ndarray = np.column_stack((x, y, z))
        grid: pv.PolyData = pv.PolyData(points)
        self.interpolated_grid: pv.PolyData = grid
        pass
    
    def interpolateStep(self, step: int) -> None:
        poly_data: pv.PolyData = self.readStep(step)
        self.interpolated_grid.clear_point_data()
        self.interpolated_grid = self.interpolated_grid.interpolate(poly_data, radius=self.h, sharpness=2)
        t: float = self.getTime(poly_data)
        ux_analytical: np.ndarray = np.zeros(self.interpolated_grid.n_points)
        for i in range(self.interpolated_grid.n_points):
            y: float = self.interpolated_grid.points[i, 1]
            ux_analytical[i] = self.analyticalSolution(y, t)
            pass
        self.interpolated_grid.point_data["UxAnalytical"] = ux_analytical
        pass
    
    def analyticalSolution(self, y: float, t: float) -> float:
        l: float = self.pipe_width
        ux: float = self.ax_nu / 2 * y * (l - y)
        for i in range(self.analytical_order):
            ux -= 4 * self.ax_nu * l**2 / np.pi**3 / (2 * i + 1)**3 * np.sin(np.pi * y / l * (2 * i + 1)) * np.exp(-(2*i+1)**2 * np.pi**2 * self.nu * t / l**2)
            pass
        return ux
        pass
    
    def referencePlot(self) -> None:
        plt.figure(figsize=(10, 5))
        color: list = ["r", "b", "g", "purple", "orange", "brown", "black"]
        percentage: list = [5, 10, 20, 30, 40, 50, 100]
        for i in range(len(percentage)):
            step: int = int(np.round(percentage[i]/100 * self.view_step))
            poly_data: pv.PolyData = self.readStep(step)
            time = self.getTime(poly_data)
            self.interpolateStep(step)
            plt.plot(
                self.interpolated_grid.point_data["Velocity"][:, 0],
                self.interpolated_grid.points[:, 1],
                color=color[i],
                label=f"SPH - {time:.4f}s",
            )
            plt.scatter(
                self.interpolated_grid.point_data["UxAnalytical"],
                self.interpolated_grid.points[:, 1],
                color=color[i],
                label=f"Analytical - {time:.2f}s",
            )
            pass
        plt.xlabel("$U_x$ (m/s)")
        plt.ylabel("$Y$ (m)")
        # legend outside the plot
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.title(f"Re = {self.reynolds_number}")
        plt.grid(True)
        plt.savefig(
            os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_reference.png"),
            bbox_inches="tight",
            dpi=300,
        )
        pass
    
    pass