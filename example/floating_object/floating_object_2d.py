'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/07/13 19:47:32
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

POST_NAME: str = "floating_object"

class FloatingObjectPostProcess:
    
    def __init__(
        self, 
        working_directory: str = f"example/{POST_NAME}", 
        key_word: str = "2d_density_500"
    ) -> None:
        self.working_directory: str = working_directory
        self.key_word: str = key_word
        self.data_path: str = os.path.join(working_directory, f"../results/{POST_NAME}/{POST_NAME}_" + self.key_word)
        self.readData()
        self.view_plot_count: int = 0
        pass
    
    def readData(self) -> None:
        self.file_list: list = os.listdir(self.data_path)
        self.view_step: int = len(self.file_list) - 1
        pass
    
    def readStep(self, step: int = 0) -> pv.PolyData:
        file: str = os.path.join(self.data_path, self.file_list[step])
        return pv.read(file)
        pass
    
    def fastplot(self, step: int, scalars: str = "Type", point_size: str = 8) -> pv.Plotter:
        plotter: pv.Plotter = pv.Plotter()
        poly_data: pv.PolyData = self.readStep(step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="coolwarm", point_size=point_size, render_points_as_spheres=True)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.show(cpos="xy", jupyter_backend="static")
        return plotter
        pass
    
    def getTime(self, poly_data: pv.PolyData) -> float:
        return poly_data.field_data["TimeValue"][0]
        pass
    
    def getStep(self, poly_data: pv.PolyData) -> int:
        return poly_data.field_data["TMSTEP"][0]
        pass
    
    def viewPlotStep(self, view_step: int = 0, scalars: str = "Type", point_size: str = 8) -> None:
        plotter: pv.Plotter = pv.Plotter(off_screen=True)
        poly_data: pv.PolyData = self.readStep(view_step)
        plotter.add_mesh(poly_data, scalars=scalars, cmap="coolwarm", point_size=point_size, render_points_as_spheres=True)
        t, step = self.getTime(poly_data), self.getStep(poly_data)
        plotter.add_text(f"Time: {t:.4f} s  Step: {step}", font_size=15, position="upper_edge")
        plotter.camera_position = "xy"
        plotter.screenshot(os.path.join(self.working_directory, f"image/{POST_NAME}_" + self.key_word + f"_cmap{self.view_plot_count}.png"))
        self.view_plot_count += 1
        pass
    
    def viewPlot(self) -> None:
        for i in [0, 10, 20, 30, 40, 50]:
            step: int = int(np.floor(i / 100 * self.view_step))
            self.viewPlotStep(step)
            pass
        pass
    
    pass