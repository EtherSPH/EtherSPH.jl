'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/06/17 15:29:22
 # @ license: MIT
 # @ description:
 '''

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
pv.set_jupyter_backend("static")
pv.start_xvfb()

file: str = "../results/poisson_equation/poisson_equation_2d/poisson_equation_2d0.vtp"
poly_data: pv.PolyData = pv.read(file)

x: np.ndarray = np.array(poly_data.points[:, 0])
y: np.ndarray = np.array(poly_data.points[:, 1])
u: np.ndarray = np.array(poly_data.point_data["U"])
u_theory: np.ndarray = np.array(poly_data.point_data["UTheory"])

plt.figure(figsize=(10, 6), facecolor="white")
ax1 = plt.subplot(1, 2, 1)
# aspect="equal"
ax1.set_aspect("equal")
cbar1 = plt.colorbar(ax1.scatter(x, y, c=u, cmap="jet", s=6), ax=ax1, orientation="horizontal")
cbar1.set_ticks(np.linspace(0, 15, 7))
ax1.set_title("SPH Numerical Solution")
ax1.set_xlabel("$x$")
ax1.set_ylabel("$y$")
ax2 = plt.subplot(1, 2, 2)
ax2.set_aspect("equal")
cbar2 = plt.colorbar(ax2.scatter(x, y, c=u_theory, cmap="jet", s=6), ax=ax2, orientation="horizontal")
cbar2.set_ticks(np.linspace(0, 15, 7))
ax2.set_title("Analytical Solution")
ax2.set_xlabel("$x$")
ax2.set_ylabel("$y$")
plt.savefig("image/poisson_equation_2d.png", bbox_inches="tight", dpi=300)