'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/18 02:27:37
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

class LidDrivenCavityPostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    lid_length: float = 1.0
    reference_middle_u: pd.DataFrame = pd.read_csv("data/middle_u.csv")
    reference_middle_v: pd.DataFrame = pd.read_csv("data/middle_v.csv")
    
    def __init__(self, reference_gap: float = 0.01, key_word: str = "compulsive", reynolds_number: int = 100) -> None:
        self.reference_gap: float = reference_gap
        self.h: float = 3 * reference_gap
        self.key_word: str = key_word
        self.reynolds_number: int = reynolds_number
        self.data_path: str = "../results/lid_driven_cavity/lid_driven_cavity_re_%d_%s" % (reynolds_number, key_word)
        self.readData()
        self.makeMiddleU()
        self.makeMiddleV()
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
    
    def fastplot(self, step: int, scalars: str = "Velocity", point_size: str = 8) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self, scalars: str = "Velocity", point_size: str = 8) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(self.view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(f"image/lid_driven_cavity_re_{self.reynolds_number}_" + self.key_word + "_cmap.png")
        pass
    
    def getPosition(self, poly_data: pv.PolyData) -> tuple[np.ndarray, np.ndarray]:
        x: np.ndarray = poly_data.points[:, 0]
        y: np.ndarray = poly_data.points[:, 1]
        return x, y
        pass
    
    def getProperty(self, poly_data: pv.PolyData, key: str="Velocity") -> np.ndarray:
        return poly_data.point_data[key]
        pass
    
    def makeMiddleU(self) -> None:
        n: int = int(np.floor(self.lid_length / self.reference_gap))
        y: np.ndarray = np.linspace(0.0, self.lid_length, n+1)
        x: np.ndarray = np.ones_like(y) * 0.5 * self.lid_length
        points: np.ndarray = np.column_stack((x, y, np.zeros_like(x)))
        poly_data: pv.PolyData = pv.PolyData(points)
        self.middle_u_poly_data: pv.PolyData = poly_data
        pass
    
    def makeMiddleV(self) -> None:
        n: int = int(np.floor(self.lid_length / self.reference_gap))
        x: np.ndarray = np.linspace(0.0, self.lid_length, n+1)
        y: np.ndarray = np.ones_like(x) * 0.5 * self.lid_length
        points: np.ndarray = np.column_stack((x, y, np.zeros_like(x)))
        poly_data: pv.PolyData = pv.PolyData(points)
        self.middle_v_poly_data: pv.PolyData = poly_data
        pass
    
    def interpolateMiddle(self, step: int) -> None:
        poly_data: pv.PolyData = self.readStep(step)
        self.middle_u_poly_data.clear_point_data()
        self.middle_v_poly_data.clear_point_data()
        self.middle_u_poly_data = self.middle_u_poly_data.interpolate(poly_data, radius=self.h, sharpness=0.2)
        self.middle_v_poly_data = self.middle_v_poly_data.interpolate(poly_data, radius=self.h, sharpness=0.2)
        pass
    
    def getReferenceMiddleYAndU(self) -> tuple[np.ndarray, np.ndarray]:
        y: np.ndarray = self.reference_middle_u["y"].values
        u: np.ndarray = self.reference_middle_u["re%d" % self.reynolds_number].values
        return y, u
        pass
    
    def getReferenceMiddleXAndV(self) -> tuple[np.ndarray, np.ndarray]:
        x: np.ndarray = self.reference_middle_v["x"].values
        v: np.ndarray = self.reference_middle_v["re%d" % self.reynolds_number].values
        return x, v
        pass
    
    def referencePlot(self) -> None:
        self.interpolateMiddle(len(self.file_list)-1)
        plt.figure(figsize=(12, 5), facecolor="white")
        ax1 = plt.subplot(1, 2, 1)
        ax1.plot(self.middle_u_poly_data.points[:, 1], self.middle_u_poly_data.point_data["Velocity"][:, 0], label="$u$ at $x=0.5$ - SPH", color="b")
        ref_y, ref_u = self.getReferenceMiddleYAndU()
        ax1.scatter(ref_y, ref_u, label="$u$ at $x=0.5$ - Reference", color="r", zorder=10)
        ax1.grid(True)
        ax1.legend()
        ax1.set_title("Middle Velocity $U$ Re=%d" % self.reynolds_number)
        ax1.set_xlabel("$y$")
        ax1.set_ylabel("$u$")
        ax2 = plt.subplot(1, 2, 2)
        ax2.plot(self.middle_v_poly_data.points[:, 0], self.middle_v_poly_data.point_data["Velocity"][:, 1], label="$v$ at $y=0.5$ - SPH", color="b")
        ref_x, ref_v = self.getReferenceMiddleXAndV()
        ax2.scatter(ref_x, ref_v, label="$v$ at $y=0.5$ - Reference", color="r", zorder=10)
        ax2.grid(True)
        ax2.legend()
        ax2.set_title("Middle Velocity $V$ Re=%d" % self.reynolds_number)
        ax2.set_xlabel("$x$")
        ax2.set_ylabel("$v$")
        plt.savefig("image/lid_driven_cavity_re_%d_" % self.reynolds_number + self.key_word + "_reference.png", bbox_inches="tight", dpi=300)
        pass
    
    pass

re100_compulsive: LidDrivenCavityPostProcess = LidDrivenCavityPostProcess(reynolds_number=100, reference_gap=0.01)
re100_compulsive.viewPlot()
re100_compulsive.referencePlot()

re400_compulsive: LidDrivenCavityPostProcess = LidDrivenCavityPostProcess(reynolds_number=400, reference_gap=0.01)
re400_compulsive.viewPlot()
re400_compulsive.referencePlot()