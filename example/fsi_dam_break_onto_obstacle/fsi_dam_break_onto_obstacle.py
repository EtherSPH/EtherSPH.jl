'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/08/28 21:25:26
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

POST_NAME: str = "fsi_dam_break_onto_obstacle"

class FSIDamBreakOntoObstaclePostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    BEAM_MOVABLE_TAG: int = 3
    BEAM_FIXED_TAG: int = 4
    
    def __init__(
        self,
        working_directory: str = f"example/{POST_NAME}",
        reference_gap: float = 0.002,
        key_word: str = "original"
    ) -> None:
        self.working_directory: str = working_directory
        self.reference_gap: float = reference_gap
        self.h: float = 3 * reference_gap
        self.key_word: str = key_word
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_{key_word}")
        self.readData()
        self.plot_times: int = 0
        pass
    
    def readData(self) -> None:
        self.file_list: list = os.listdir(self.data_path)
        # remove end with .csv
        self.file_list = [file for file in self.file_list if not file.endswith(".csv")]
        self.fluid_file_list: list = [file for file in self.file_list if file.startswith("fluid")]
        self.wall_file_list: list = [file for file in self.file_list if file.startswith("wall")]
        self.beam_movable_file_list: list = [file for file in self.file_list if file.startswith("beam_movable")]
        self.beam_fixed_file_list: list = [file for file in self.file_list if file.startswith("beam_fixed")]
        self.file_list = [file for file in self.file_list if file.endswith(".vtp")]
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
    
    def readFluidStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.fluid_file_list[step])
        return pv.read(file)
        pass
    
    def readWallStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.wall_file_list[step])
        return pv.read(file)
        pass
    
    def readBeamMovableStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.beam_movable_file_list[step])
        return pv.read(file)
        pass
    
    def readBeamFixedStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.beam_fixed_file_list[step])
        return pv.read(file)
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
    
    def viewPlotEachStep(self, step: int = 0, point_size_1: float = 3, point_size_2: float = 1) -> None:
        self.plot_times += 1
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        fluid_poly_data: pv.PolyData = self.readStep(step)
        wall_poly_data: pv.PolyData = self.readWallStep(step)
        beam_movable_poly_data: pv.PolyData = self.readBeamMovableStep(step)
        beam_fixed_poly_data: pv.PolyData = self.readBeamFixedStep(step)
        plotter.add_mesh(fluid_poly_data, scalars="Velocity", cmap="coolwarm", point_size=point_size_1)
        plotter.add_mesh(wall_poly_data, color="gray")
        plotter.add_mesh(beam_movable_poly_data, scalars="VonMisesStress", cmap="jet", point_size=point_size_2)
        plotter.add_mesh(beam_fixed_poly_data, scalars="VonMisesStress", cmap="jet", point_size=point_size_2)
        t, step = self.getTime(fluid_poly_data), self.getStep(fluid_poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + f"_cmap_{self.plot_times}.png"))
        pass
    
    def viewPlot(self, point_size_1: float = 3, point_size_2: float = 1) -> None:
        self.plot_times = 0
        for step in [25, 30, 35, 40, 45, 50, 120, 130]:
            self.viewPlotEachStep(step, point_size_1, point_size_2)
            pass
        pass
    
    def getReference(self, name: str) -> tuple:
        file_name: str = os.path.join(self.working_directory, f"data/{name}.csv")
        df: pd.DataFrame = pd.read_csv(file_name)
        ts: np.ndarray = df["Time"].values
        xs: np.ndarray = df["X"].values
        return ts, xs / 1e3
        pass
    
    def referencePlot(self) -> None:
        df: pd.DataFrame = pd.read_csv(os.path.join(self.working_directory, f"../results/{POST_NAME}/{POST_NAME}_{self.key_word}/monitor_point.csv"))
        ts: np.ndarray = df["Time"].values
        xs: np.ndarray = df["X"].values
        ys: np.ndarray = df["Y"].values
        plt.figure(figsize=(8, 4))
        plt.plot(ts, xs - xs[0], label="My Code: X", color="k", zorder=10, lw=3)
        plt.xlabel("Time [s]")
        plt.ylabel("Displacement [m]")
        for name in ["li", "marti", "rafiee", "rahimi", "rahimi_sph"]:
            ts_ref, xs_ref = self.getReference(name)
            plt.plot(ts_ref, xs_ref - xs_ref[0], label=f"{name.capitalize()}: X", marker="o", markersize=3)
            pass
        plt.legend()
        plt.grid()
        plt.xlim(0, 0.5)
        plt.ylim(-0.001, 0.05)
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_reference.png"), bbox_inches="tight", dpi=300)
        pass
    
    pass