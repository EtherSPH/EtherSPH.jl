'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/16 01:47:19
 # @ license: MIT
 # @ description:
 '''

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

working_directory: str = "example/poisson_equation/"

df_1d: pd.DataFrame = pd.read_csv(os.path.join(working_directory, "../results/poisson_equation/poisson_equation_1d/poisson_equation_1d.csv"))

plt.figure(figsize=(6, 4), facecolor="white")
plt.plot(df_1d["x"], df_1d["u"], label="SPH Numerical Solution", color="blue")
plt.scatter(df_1d["x"][0:-1:4], df_1d["u_theory"][0:-1:4], label="Analytical Solution", color="red", s=20, zorder=10)
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.grid(True)
plt.title("1D Poisson Equation $3x^2+10x+2$")
plt.savefig(os.path.join("example/poisson_equation/", "image/poisson_equation_1d.png"), bbox_inches="tight", dpi=300)