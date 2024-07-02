'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/07/02 16:40:39
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

POST_NAME: str = "cylinder_2d"

class Cylinder2DPostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    CYLINDER_TAG: int = 3
    
    def __init__(self, working_directory: str = f"example/cylinder",  reference_gap: float = 0.41 / 82, key_word: str = "re_20_same_periodic") -> None:
        self.working_directory: str = working_directory
        self.reference_gap: float = reference_gap
        self.h: float = 3 * reference_gap
        self.key_word: str = key_word
        self.data_path: str = os.path.join(working_directory, f"../results/cylinder/{POST_NAME}_{key_word}")
        self.readData()
        pass
    
    def readData(self) -> None:
        self.file_list: list = os.listdir(self.data_path)
        # remove end with .csv
        self.file_list = [file for file in self.file_list if not file.endswith(".csv")]
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
    
    def fastplot(self, step: int, scalars: str = "Velocity", point_size: str = 3) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self, scalars: str = "Velocity", point_size: str = 3) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(self.view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_cmap.png"))
        pass
    
    pass