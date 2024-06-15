# Poisson Equation

$$
\begin{equation}
    \nabla^2 u + f = 0 \to
    \sum_j\frac{2m_j}{\rho_j}\frac{W^\prime_{ij}}{r_{ij}}(u_i-u_j) 
    = -\sum_j \frac{m_j}{\rho_j} f_j W_{ij}
\end{equation}
$$

# 1D problem

Take a function as example:

$$
\begin{equation}
    u(x) = 3 x^2 + 10 x + 2
\end{equation}
$$

Which satisfies the Poisson equation:

$$
\begin{equation}
    \frac{\mathrm{d}^2 u}{\mathrm{d}x^2} - 6 = 0
\end{equation}
$$

Use SPH method to discretize the equation:

$$
\begin{equation}
    \sum_j\frac{2}{\Delta x}\frac{W^\prime(r_{ij})}{r_{ij}}(u_i-u_j) 
    = \sum_j \frac{6}{\Delta x} W(r_{ij})
\end{equation}
$$

Choose `WendlandC4` kernel function, and the smoothing length is $h=3\Delta x$, we will have:

<center>
<img src="image/poisson_equation_1d.png" width=70%>
</br>
fig. 1D Poisson equation demo
</center>

This case indicates that the SPH method can solve the Poisson equation in 1D space. In addition, a regular distribution of particles leads to a better result, as shown in the figure above.