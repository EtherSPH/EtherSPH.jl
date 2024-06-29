'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/15 16:47:57
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

POST_NAME: str = "collapse_dry"

class CollapseDryPostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    water_width: float = 1.0
    water_height: float = 2.0
    gravity: float = 9.81
    view_step: int = 200
    
    def __init__(self, working_directory: str = f"example/{POST_NAME}", referece_gap: float = 0.02, key_word: str = "same") -> None:
        self.working_directory: str = working_directory
        self.referece_gap: float = referece_gap
        self.h: float = 3 * self.referece_gap
        self.key_word: str = key_word
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_" + self.key_word)
        self.readData()
        self.readReference()
        pass
    
    def readReference(self) -> None:
        self.koshizuka_and_oka_1996_width: pd.DataFrame = pd.read_csv(os.path.join(self.working_directory, "data/koshizuka_and_oka_1996_width.csv"))
        self.koshizuka_and_oka_1996_height: pd.DataFrame = pd.read_csv(os.path.join(self.working_directory, "data/koshizuka_and_oka_1996_height.csv"))
        self.violeau_and_issa_2007_width: pd.DataFrame = pd.read_csv(os.path.join(self.working_directory, "data/violeau_and_issa_2007_width.csv"))
        self.violeau_and_issa_2007_height: pd.DataFrame = pd.read_csv(os.path.join(self.working_directory, "data/violeau_and_issa_2007_height.csv"))
        pass
    
    def readData(self) -> None:
        self.file_list: list = os.listdir(self.data_path)
        pass
    
    def readStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.file_list[step])
        return pv.read(file)
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
    
    def getTime(self, poly_data: pv.PolyData) -> float:
        return poly_data.field_data["TimeValue"][0]
        pass
    
    def getStep(self, poly_data: pv.PolyData) -> int:
        return poly_data.field_data["TMSTEP"][0]
        pass
    
    def getPosition(self, poly_data: pv.PolyData) -> tuple[np.ndarray, np.ndarray]:
        x: np.ndarray = poly_data.points[:, 0]
        y: np.ndarray = poly_data.points[:, 1]
        return x, y
        pass
    
    def getProperty(self, poly_data: pv.PolyData, key: str="Velocity") -> np.ndarray:
        return poly_data.point_data[key]
        pass
    
    def getWidthAndHeight(self, poly_data: pv.PolyData) -> tuple[float, float]:
        x, y = self.getPosition(poly_data)
        types = self.getProperty(poly_data, "Type")
        width: float = 0.0
        height: float = 0.0
        for i in range(len(x)):
            if types[i] == self.WALL_TAG:
                continue
            else:
                width = max(width, x[i])
                if self.h < x[i] < self.water_width - self.h:
                    height = max(height, y[i])
                    pass 
                pass
            pass
        return width, height
        pass
    
    def getCurveForReference(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ts: list = []
        ws: list = []
        hs: list = []
        for step in range(len(self.file_list)):
            poly_data: pv.PolyData = self.readStep(step)
            t: float = self.getTime(poly_data)
            if t > 1.0:
                break
                pass
            w, h = self.getWidthAndHeight(poly_data)
            ts.append(t)
            ws.append(w)
            hs.append(h)
            pass
        ts: np.ndarray = np.array(ts) * np.sqrt(2 * self.gravity / self.water_width)
        ws: np.ndarray = np.array(ws) / self.water_width
        hs: np.ndarray = np.array(hs) / self.water_height
        return ts, ws, hs
        pass
    
    def referencePlot(self) -> None:
        ts, ws, hs = self.getCurveForReference()
        marker_size: int = 50
        plt.figure(figsize=(10, 10), facecolor="white")
        plt.subplot(2, 1, 1)
        plt.plot(ts, ws, label="Collapse Dry Width", color="blue")
        plt.scatter(
            self.koshizuka_and_oka_1996_width["time"],
            self.koshizuka_and_oka_1996_width["width"],
            label="Koshizuka and Oka 1996 Width",
            color="red",
            marker="o",
            zorder=10,
            s=marker_size
        )
        plt.scatter(
            self.violeau_and_issa_2007_width["time"],
            self.violeau_and_issa_2007_width["width"],
            label="Violeau and Issa 2007 Width",
            color="green",
            marker="x",
            zorder=10,
            s=marker_size
        )
        plt.xlabel("t*")
        plt.ylabel("W/W0")
        plt.legend()
        plt.grid(True)
        plt.subplot(2, 1, 2)
        plt.plot(ts, hs, label="Collapse Dry Height", color="blue")
        plt.scatter(
            self.koshizuka_and_oka_1996_height["time"],
            self.koshizuka_and_oka_1996_height["height"],
            label="Koshizuka and Oka 1996 Height",
            color="red",
            marker="o",
            zorder=10,
            s=marker_size
        )
        plt.scatter(
            self.violeau_and_issa_2007_height["time"],
            self.violeau_and_issa_2007_height["height"],
            label="Violeau and Issa 2007 Height",
            color="green",
            marker="x",
            zorder=10,
            s=marker_size
        )
        plt.xlabel("t*")
        plt.ylabel("H/H0")
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_reference.png"), bbox_inches="tight", dpi=300)
        pass
    
    pass