'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/30 02:01:11
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

POST_NAME: str = "cruchaga"

class CruchagaPostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    
    def __init__(self, working_directory: str = f"example/{POST_NAME}", referece_gap: float = 0.002, key_word: str = "2d") -> None:
        self.working_directory: str = working_directory
        self.referece_gap: float = referece_gap
        self.h: float = 3 * self.referece_gap
        self.key_word: str = key_word
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_" + self.key_word)
        self.readData()
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
    
    def fastplot(self, step: int, scalars: str = "Velocity", point_size: str = 4) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="coolwarm", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self) -> None:
        item = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        for i in range(len(item)):
            step = int(np.floor((item[i] / 100 * self.view_step)))
            self.viewPlotStep(i, step)
            pass
        pass
    
    def viewPlotStep(self, i: int, view_step: int, scalars: str = "Velocity", point_size: str = 5) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="coolwarm", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + f"_step_{i}_cmap.png"))
        pass
    
    def getReferencePlotData(self) -> tuple:
        t: list = []
        x: list = []
        step: int = 0
        current_t: float = 0.0
        while current_t < 0.3:
            poly_data: pv.PolyData = self.readStep(step)
            current_t = self.getTime(poly_data)
            # points with tag 1
            mask: np.ndarray = poly_data.point_data["Type"] == self.FLUID_TAG
            new_poly_data: pv.PolyData = poly_data.extract_points(mask)
            xmax: float = new_poly_data.points[:, 0].max()
            x.append(xmax)
            t.append(current_t)
            step += 1
            pass
        return np.array(t), np.array(x) / 0.114
        pass
    
    def getReferenceData(self) -> tuple:
        t: np.ndarray = np.array([0.0, 0.1, 0.2, 0.3])
        x: np.ndarray = np.array([1.0, 1.55, 2.6, 3.7])
        return t, x
        pass
    
    def referencePlot(self) -> None:
        t, x = self.getReferencePlotData()
        t_r, x_r = self.getReferenceData()
        plt.figure(figsize=(10, 5), facecolor="white")
        plt.plot(t, x, label="SPH - Cruchaga 2D", color="blue")
        plt.scatter(t_r, x_r, color="r", marker="x", s=90, zorder=10, label="Experimental - Cruchaga 3D")
        plt.xlabel("Time [s]")
        plt.ylabel("Dimensionless Horizontal Displacement")
        plt.title("Cruchaga 2D - SPH vs Experimental 3D")
        plt.grid()
        plt.legend()
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_{self.key_word}_reference.png"), dpi=300, bbox_inches="tight")
        pass
        
    pass