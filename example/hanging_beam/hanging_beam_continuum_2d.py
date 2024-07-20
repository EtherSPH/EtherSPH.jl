'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/07/20 17:50:16
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

POST_NAME = "hanging_beam"

class HangingBeam2DPostProcess:
    
    MATERIAL_MOVABLE_TAG = 1
    MATERIAL_FIXED_TAG = 2
    
    def __init__(self, working_directory: str = f"example/{POST_NAME}", reference_gap: float = 0.01, key_word: str = "continuum_2d") -> None:
        self.working_directory: str = working_directory
        self.reference_gap: float = reference_gap
        self.h: float = 2.6 * reference_gap
        self.key_word: str = key_word
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_{key_word}")
        self.readData()
        self.monitor_point_index: int = self.getMonitorPointIndex()
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
    
    def fastplot(self, step: int, scalars: str = "VonMisesStress", point_size: str = 8) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self) -> None:
        for plot_step in [0, 20, 40, 60]:
            self.viewPlotStep(plot_step)
            pass
        pass
    
    def viewPlotStep(self, plot_step: int = 0, scalars: str = "VonMisesStress", point_size: str = 8) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(plot_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.add_mesh(poly_data.extract_points([self.monitor_point_index]), color="red", point_size=10)
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_{self.key_word}_step_{plot_step}" + "_cmap.png"))
        pass
    
    def getMonitorPointIndex(self) -> int:
        poly_data: pv.PolyData = self.readStep(0)
        length: float = 0.2
        width: float = 0.02
        for i in range(poly_data.n_points):
            x, y, _ = poly_data.points[i]
            if x < length and x > length - self.reference_gap:
                if abs(y - width / 2) <= self.reference_gap:
                    return i
                    pass
                else:
                    continue
                    pass
                pass
            else:
                continue
                pass
            pass
        return 
        pass
    
    def referencePlot(self) -> None:
        ts: list = []
        ys: list = []
        for i in range(self.view_step + 1):
            poly_data: pv.PolyData = self.readStep(i)
            ts.append(self.getTime(poly_data))
            ys.append(poly_data.points[self.monitor_point_index][1])
            pass
        ts: np.ndarray = np.array(ts)
        ys: np.ndarray = np.array(ys) - 0.02 / 2
        T: float = 0.25406
        ts_ref: np.ndarray = np.linspace(0, 1.0, 81)
        ys_ref: np.ndarray = np.sin(2 * np.pi * ts_ref / T)
        ys_ref = (ys.max() - ys.min()) / 2 * ys_ref + (ys.max() + ys.min()) / 2
        plt.figure(figsize=(8, 6))
        plt.plot(ts, ys, color="blue", linewidth=2, label="Monitor Point - SPH")
        plt.scatter(ts_ref, ys_ref, color="red", label="Monitor Point - Analytical", zorder = 10, s=30)
        plt.xlabel("Time (s)")
        plt.ylabel("Displacement (m)")
        plt.title("Monitor Point Displacement")
        plt.grid()
        plt.legend(loc="upper right")
        plt.ylim(-0.025, 0.030)
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_{self.key_word}_comparison.png"), bbox_inches="tight", dpi=300)
        pass
    
    pass