'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/08/27 17:22:55
 # @ license: MIT
 # @ description:
 '''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.fft import fft, fftfreq
pv.set_jupyter_backend("static")
pv.start_xvfb()

POST_NAME: str = "fsi_cylinder_with_plate"

class FSICylinderWithPlatePostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    MOVABLE_SOLID_TAG: int = 3
    FIXED_SOLID_TAG: int = 4
    
    def __init__(
        self,
        working_directory: str = f"example/fsi_cylinder_with_plate",
        reference_gap: float = 0.1 / 20,
        key_word: str = "case_1",
    ) -> None:
        self.working_directory: str = working_directory
        self.reference_gap: float = reference_gap
        self.h: float = 3 * reference_gap
        self.key_word: str = key_word
        self.data_path: str = os.path.join(working_directory, f"../results/fsi_cylinder_with_plate/{POST_NAME}_{key_word}")
        self.readData()
        self.plot_times: int = 0
        pass
    
    def readData(self) -> None:
        self.file_list: list = os.listdir(self.data_path)
        # remove end with .csv
        self.file_list = [file for file in self.file_list if not file.endswith(".csv")]
        self.fluid_file_list: list = [file for file in self.file_list if file.startswith("fluid")]
        self.wall_file_list: list = [file for file in self.file_list if file.startswith("wall")]
        self.solid_movable_file_list: list = [file for file in self.file_list if file.startswith("solid_movable")]
        self.solid_fixed_file_list: list = [file for file in self.file_list if file.startswith("solid_fixed")]
        self.file_list = [file for file in self.file_list if file.endswith(".vtp")]
        self.view_step: int = len(self.file_list) - 1
        pass
    
    def readStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.file_list[step])
        return pv.read(file)
        pass
    
    def readFluidStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.fluid_file_list[step])
        return pv.read(file)
        pass
    
    def readWallStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.wall_file_list[step])
        return pv.read(file)
        pass
    
    def readSolidMovableStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.solid_movable_file_list[step])
        return pv.read(file)
        pass
    
    def readSolidFixedStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.solid_fixed_file_list[step])
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
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlotEachStep(self, step: int = 0, point_size_1: float = 4, point_size_2: float = 2) -> None:
        self.plot_times += 1
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        fluid_poly_data: pv.PolyData = self.readStep(step)
        wall_poly_data: pv.PolyData = self.readWallStep(step)
        solid_movable_poly_data: pv.PolyData = self.readSolidMovableStep(step)
        solid_fixed_poly_data: pv.PolyData = self.readSolidFixedStep(step)
        plotter.add_mesh(fluid_poly_data, scalars="Velocity", cmap="jet", point_size=point_size_1)
        plotter.add_mesh(wall_poly_data, color="gray", point_size=point_size_1)
        plotter.add_mesh(solid_movable_poly_data, scalars="VonMisesStress", cmap="coolwarm", point_size=point_size_2)
        plotter.add_mesh(solid_fixed_poly_data, color="gray", point_size=point_size_2)
        t, step = self.getTime(fluid_poly_data), self.getStep(fluid_poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + f"_cmap_{self.plot_times}.png"))
        pass
    
    def viewPlot(self, point_size_1: float = 4, point_size_2: float = 2) -> None:
        for i in range(6):
            step = 200 + i
            self.viewPlotEachStep(step, point_size_1, point_size_2)
            pass
        pass
    
    def curvePlot(self) -> None:
        df: pd.DataFrame = pd.read_csv(os.path.join(self.working_directory, f"../results/fsi_cylinder_with_plate/{POST_NAME}_{self.key_word}/monitor_point.csv"))
        ts: np.ndarray = df["Time"].values
        xs: np.ndarray = df["X"].values
        ys: np.ndarray = df["Y"].values
        D: float = 0.1
        U: float = 1.0
        ts *= U / D
        xs /= D
        ys /= D
        # * figure 1
        plt.figure(figsize=(12, 4), facecolor="white")
        plt.plot(ts, xs)
        plt.xlabel(r"Time: $t\bar{U}/D$")
        plt.ylabel(r"Displacement: $x/D$")
        plt.title("Displacement in x direction of the monitor point")
        plt.grid()
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_displacement_x.png"), bbox_inches="tight", dpi=300)
        # * figure 2
        y_max: float = np.max(ys)
        y_min: float = np.min(ys)
        amp: float = (y_max - y_min) / 2
        plt.figure(figsize=(12, 4), facecolor="white")
        plt.plot(ts, ys, label="Amplitude: {:.4f}".format(amp))
        plt.xlabel(r"Time: $t\bar{U}/D$")
        plt.ylabel(r"Displacement: $y/D$")
        plt.title("Displacement in y direction of the monitor point")
        plt.grid()
        plt.legend()
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_displacement_y.png"), bbox_inches="tight", dpi=300)
        # * figure 3
        plt.figure(figsize=(8, 8), facecolor="white")
        plt.gca().set_aspect("equal")
        plt.plot(xs, ys, lw=0.3)
        plt.xlabel(r"Displacement: $x/D$")
        plt.ylabel(r"Displacement: $y/D$")
        plt.title("Displacement in x-y direction of the monitor point")
        plt.grid()
        plt.xticks(np.arange(5, 7, 0.5))
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_displacement_xy.png"), bbox_inches="tight", dpi=300)
        # * figure 4
        sample_rate: float = 1 / (ts[1] - ts[0])
        n: int = len(ts)
        freqs: np.ndarray = fftfreq(n, 1 / sample_rate)
        yf: np.ndarray = fft(ys)
        freqs = freqs[:n // 2]
        yf = np.abs(yf[:n // 2])
        index = freqs > 0.1
        freqs = freqs[index]
        yf = yf[index]
        max_freq_index: int = np.argmax(yf)
        max_freq: float = freqs[max_freq_index]
        plt.figure(figsize=(12, 4), facecolor="white")
        plt.plot(freqs, yf, label="Max Frequency: {:.4f} Hz".format(max_freq))
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Amplitude")
        plt.title("Frequency Spectrum of the monitor point")
        plt.xlim(0, 1)
        plt.grid()
        plt.legend()
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_frequency_spectrum.png"), bbox_inches="tight", dpi=300)
        pass
    
    pass