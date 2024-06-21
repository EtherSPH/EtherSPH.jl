'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/21 14:24:35
 # @ license: MIT
 # @ description:
 '''

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
import scipy as scp
pv.set_jupyter_backend("static")
pv.start_xvfb()

POST_NAME: str = "natural_convection_ring"

class NaturalConvectionRingPostProcess:
    
    FLUID_TAG: int = 1
    THERMOSTATIC_WALL_TAG: int = 2
    radius_outer: float = 1.0
    
    def __init__(self, working_directory: str = f"example/{POST_NAME}", reference_gap: float = 0.01, key_word: str = "ratio_2.6_ra_4.7e4", rayleigh_number: str = "4.7e4"):
        self.working_directory: str = working_directory
        self.reference_gap: float = reference_gap
        self.h: float = 3 * self.reference_gap
        self.key_word: str = key_word
        self.rayleigh_number_str: str = rayleigh_number
        self.rayleigh_number: float = float(rayleigh_number)
        self.radius_ratio: float = float(re.findall(r"ratio_(\d+\.\d+)", key_word)[0])
        self.radius_inner: float = self.radius_outer / self.radius_ratio
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_{self.key_word}")
        self.readData()
        self.makeInterpolationGrids()
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
    
    def fastplot(self, step: int, scalars: str = "Temperature", point_size: str = 4) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="nipy_spectral", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self, scalars: str = "Temperature", point_size: str = 5) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(self.view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="nipy_spectral", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_{self.key_word}_cmap.png"))
        pass
    
    def makeInterpolationGrids(self) -> None:
        n: int = int(np.floor((self.radius_outer - self.radius_inner) / self.reference_gap)) + 1
        theta_array: np.ndarray = np.linspace(0, np.pi, 7) + np.pi / 2
        grid_list: list = []
        r: np.ndarray = np.linspace(self.radius_inner, self.radius_outer, n)
        for theta in theta_array:
            x: np.ndarray = r * np.cos(theta)
            y: np.ndarray = r * np.sin(theta)
            z: np.ndarray = np.zeros_like(x)
            points: np.ndarray = np.column_stack((x, y, z))
            grid: pv.PolyData = pv.PolyData(points)
            grid_list.append(grid)
            pass
        self.interpolated_r: np.ndarray = r
        self.interpolated_theta: np.ndarray = theta_array
        self.interpolated_benchmark_r: np.ndarray = (r - self.radius_inner) / (self.radius_outer - self.radius_inner)
        self.interpolated_grid_list: list[pv.PolyData] = grid_list
        pass
    
    def readBenchmarFile(self, index: int) -> np.ndarray:
        file: str = self.working_directory + f"/data/x{index}.dat"
        return np.loadtxt(file, delimiter=",").T
        pass
    
    def getInterpolatedStep(self, step: int) -> list[pv.PolyData]:
        poly_data: pv.PolyData = self.readStep(step)
        for i in range(len(self.interpolated_grid_list)):
            self.interpolated_grid_list[i].clear_point_data()
            self.interpolated_grid_list[i] = self.interpolated_grid_list[i].interpolate(poly_data, radius=self.h, sharpness=2)
            pass
        return self.interpolated_grid_list
        pass
    
    def referencePlot(self) -> None:
        self.getInterpolatedStep(self.view_step)
        plt.figure(figsize=(10, 6))
        theta_deg: np.ndarray = np.linspace(90, 270, 7)
        for i in range(len(self.interpolated_grid_list)):
            plt.plot(self.interpolated_benchmark_r, self.interpolated_grid_list[i].point_data["Temperature"], label=f"SPH - {theta_deg[i]} deg")
            bench_data: np.ndarray = self.readBenchmarFile(i)
            plt.scatter(bench_data[0], bench_data[1], label=f"Benchmark - {theta_deg[i]} deg", zorder=10)
            pass
        plt.xlabel("$R*$")
        plt.ylabel("$T^*$")
        plt.grid(True)
        # legnd outside
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_{self.key_word}_reference.png"), bbox_inches="tight", dpi=300)
        plt.title(f"Natural Convection Ring Benchmark - Ra={self.rayleigh_number_str}")
        pass
    
    pass