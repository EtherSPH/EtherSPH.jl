'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/07/04 22:22:13
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

class PoiseuilleFlow3DPostProcess:
    
    FLUID_TAG: int = 1
    WALL_TAG: int = 2
    
    def __init__(
        self,
        working_directory: str = f"example/{POST_NAME}",
        reference_gap: float = 0.01,
        key_word: str = "3d_re_100_same_periodic",
        reynolds_number: float = 100,
    ) -> None:
        self.working_directory: str = working_directory
        self.reference_gap: float = reference_gap
        self.h: float = 2 * reference_gap
        self.key_word: str = key_word
        self.reynolds_number: float = reynolds_number
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_{key_word}")
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
    
    def fastplot(self, step: int, scalars: str = "Velocity", point_size: str = 10) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size, nan_opacity=0.1, clim=[0, 0.02])
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(jupyter_backend="static")
        return plotter
        pass
    
    def viewPlot(self, scalars: str = "Velocity", point_size: int = 10) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(self.view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="jet", point_size=point_size, nan_opacity=0.1, clim=[0, 0.02])
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + "_cmap.png"))
        pass
    
    pass