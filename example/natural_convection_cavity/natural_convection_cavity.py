'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/18 21:53:52
 # @ license: MIT
 # @ description:
 '''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
import scipy as scp
pv.set_jupyter_backend("static")
pv.start_xvfb()

POST_NAME: str = "natural_convection_cavity"

class NaturalConvectionCavityPostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    THERMOSTATIC_WALL_TAG: int = 3
    cavity_length: float = 1.0
    
    def __init__(self, working_directory: str = f"example/{POST_NAME}", reference_gap: float = 0.01, key_word: str = "compulsive", rayleigh_number: str = "1e3"):
        self.working_directory: str = working_directory
        self.reference_gap: float = reference_gap
        self.h: float = reference_gap * 3
        self.key_word: str = key_word
        self.rayleigh_number_str: str = rayleigh_number
        self.rayleigh_number: float = float(rayleigh_number)
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_ra_{self.rayleigh_number_str}_{self.key_word}")
        self.readData()
        self.makeInterpolationGrid()
        pass
    
    def readData(self) -> None:
        self.file_list: list = os.listdir(self.data_path)
        self.view_step: int = len(self.file_list) - 1
        pass
    
    def getTime(self, poly_data: pv.PolyData) -> float:
        return poly_data.field_data["TimeValue"][0]
        pass
    
    def getStep(self, poly_data: pv.PolyData) -> int:
        return poly_data.field_data["TMSTEP"][0]
        pass
    
    def readStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.file_list[step])
        return pv.read(file)
        pass
    
    def fastplot(self, step: int, scalars: str = "Temperature", point_size: str = 8) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="coolwarm", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self, scalars: str = "Temperature", point_size: str = 8) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(self.view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="coolwarm", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_ra_{self.rayleigh_number_str}_{self.key_word}_cmap.png"))
        pass
    
    def makeInterpolationGrid(self) -> None:
        n: int = int(np.floor(self.cavity_length / self.reference_gap) + 1)
        x: np.ndarray = np.linspace(0, self.cavity_length, n)
        y: np.ndarray = np.linspace(0, self.cavity_length, n)
        z: np.ndarray = np.array([0.0])
        x, y, z = np.meshgrid(x, y, z, indexing="ij")
        self.interpolation_grid: pv.StructuredGrid = pv.StructuredGrid(x, y, z)
        pass
    
    def removeNanTemperaturePoint(self, poly_data: pv.PolyData) -> pv.PolyData:
        mask: np.ndarray = np.isnan(poly_data.point_data["Temperature"])
        keep_list: np.ndarray = np.arange(poly_data.n_points)[~mask]
        poly_data = poly_data.extract_points(keep_list)
        return poly_data
        pass
    
    def getInterpolatedStep(self, step: int) -> pv.PolyData:
        self.interpolation_grid.clear_point_data()
        poly_data: pv.PolyData = self.readStep(step)
        poly_data = self.removeNanTemperaturePoint(poly_data)
        self.interpolation_grid = self.interpolation_grid.interpolate(poly_data, radius=self.h, sharpness=2)
        return self.interpolation_grid
        pass
    
    def referencePlot(self) -> None:
        poly_data: pv.PolyData = self.getInterpolatedStep(self.view_step)
        n: int = int(np.floor(self.cavity_length / self.reference_gap) + 1)
        x = poly_data.points[:, 0].reshape(-1, n)
        y = poly_data.points[:, 1].reshape(-1, n)
        u = poly_data.point_data["Velocity"][:, 0].reshape(-1, n)
        v = poly_data.point_data["Velocity"][:, 1].reshape(-1, n)
        t = poly_data.point_data["Temperature"].reshape(-1, n)
        plt.figure(figsize=(12, 6), facecolor="white")
        ax1 = plt.subplot(121)
        ax1.set_aspect("equal")
        cf1 = ax1.contourf(x, y, t, cmap="coolwarm", levels=100)
        plt.colorbar(cf1, ticks=np.linspace(0, 1, 11))
        ax1.streamplot(x, y, u, v, color="black", density=1, linewidth=0.3)
        ax1.set_xlim(0, self.cavity_length)
        ax1.set_ylim(0, self.cavity_length)
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$y$")
        ax1.set_title("Temperature Field & Streamlines")
        ax2 = plt.subplot(122)
        ax2.set_aspect("equal")
        cf2 = ax2.contourf(x, y, t, cmap="coolwarm", levels=100)
        c2 = ax2.contour(x, y, t, levels=11, linewidths=0.3, colors="black")
        ax2.clabel(c2, inline=True, fontsize=10, fmt="%.1f")
        plt.colorbar(cf2, ticks=np.linspace(0, 1, 11))
        ax2.set_xlim(0, self.cavity_length)
        ax2.set_ylim(0, self.cavity_length)
        ax2.set_xlabel("$x$")
        ax2.set_ylabel("$y$")
        ax2.set_title("Temperature Field & Contours")
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_ra_{self.rayleigh_number_str}_{self.key_word}_reference.png"), bbox_inches="tight", dpi=300)
        pass
    
    def calculateNusseltNumber(self) -> float:
        n: int = int(np.floor(self.cavity_length / self.reference_gap) + 1)
        poly_data: pv.PolyData = self.getInterpolatedStep(self.view_step)
        poly_data = poly_data.compute_derivative(scalars="Temperature", gradient="TemperatureGradient")
        y: np.ndarray = np.array(poly_data.points[0:-1:n][:, 1])
        dtdx: np.array = np.array(poly_data.point_data["TemperatureGradient"][0:-1:n][:, 0])
        nu: float = -np.trapz(dtdx, y)
        return nu
        pass
    
    pass