[toc]

# Floating Object

Particle method is amazing that it will automatically handle the floating object as long as the density of fluid and rigid body are correctly set. The floating object will be pushed up by the pressure of fluid, and the gravity will pull it down. The balance of these two forces will determine the floating height of the object.

As rigid body particles get close to the fluid particles, a balanced continuity equation will increas the density of fluid particles:

$$
\begin{equation}
    \begin{aligned}
        \left(\frac{\mathrm{d}\rho}{\mathrm{d}t}\right)_f=
        \sum_s \rho_f \frac{m_s}{\rho_s}\vec{v}_{fs}\cdot \nabla W_{fs}
    \end{aligned}
\end{equation}
$$

Pressure should be extrapolated to get rigid body particles' pressure:

$$
\begin{equation}
    \begin{aligned}
    p_s=\sum_f \frac{\frac{m_f}{\rho_f}W_{pf}(p_f+\rho_f\vec{g}\cdot\vec{r}_{sf})}{\frac{m_f}{\rho_f}W_{pf}}
    \end{aligned}
\end{equation}
$$

With pressure extrapolation from fluid particles to rigid body particles, pressure force will prevent fluid particles from penetrating the rigid body:

$$
\begin{equation}
    \begin{aligned}
        m_p\vec{f}_{f}^p = -\sum_s \frac{m_pm_s}{\rho_s\rho_f}
        (p_f+p_s)\nabla W_{pf}
    \end{aligned}
\end{equation}
$$

It's easy to find that the force applied on fluid and rigid body particles are equal in magnitude but opposite in direction, which obeys Newton's third law.

# Different density of floating object in 2D case

From middle school, we know that the floating object will sink if its density is larger than the fluid. The following figures show the floating object with different density ratios. The density ratio is defined as the density of floating object divided by the density of fluid.

Assuming $r=\frac{\rho_{\text{rigid body}}}{\rho_f}$, the little cube will sink $r\times 100\%$ volume into the fluid. When $r > 1$, the cube will sink totally into the fluid.

## Density ratio = 0.2

<center>
<img src="image/floating_object_2d_density_200_cmap0.png" width=45%>
<img src="image/floating_object_2d_density_200_cmap1.png" width=45%>
</br>
<img src="image/floating_object_2d_density_200_cmap2.png" width=45%>
<img src="image/floating_object_2d_density_200_cmap3.png" width=45%>
</br>
<img src="image/floating_object_2d_density_200_cmap4.png" width=45%>
<img src="image/floating_object_2d_density_200_cmap5.png" width=45%>
</br>
fig. floating object density = 0.2 fluid density
</center>

## Density ratio = 0.5

<center>
<img src="image/floating_object_2d_density_500_cmap0.png" width=45%>
<img src="image/floating_object_2d_density_500_cmap1.png" width=45%>
</br>
<img src="image/floating_object_2d_density_500_cmap2.png" width=45%>
<img src="image/floating_object_2d_density_500_cmap3.png" width=45%>
</br>
<img src="image/floating_object_2d_density_500_cmap4.png" width=45%>
<img src="image/floating_object_2d_density_500_cmap5.png" width=45%>
</br>
fig. floating object density = 0.5 fluid density
</center>

## Density ratio = 0.8

<center>
<img src="image/floating_object_2d_density_800_cmap0.png" width=45%>
<img src="image/floating_object_2d_density_800_cmap1.png" width=45%>
</br>
<img src="image/floating_object_2d_density_800_cmap2.png" width=45%>
<img src="image/floating_object_2d_density_800_cmap3.png" width=45%>
</br>
<img src="image/floating_object_2d_density_800_cmap4.png" width=45%>
<img src="image/floating_object_2d_density_800_cmap5.png" width=45%>
</br>
fig. floating object density = 0.8 fluid density
</center>

## Density ratio = 2.0

<center>
<img src="image/floating_object_2d_density_2000_cmap0.png" width=45%>
<img src="image/floating_object_2d_density_2000_cmap1.png" width=45%>
</br>
<img src="image/floating_object_2d_density_2000_cmap2.png" width=45%>
<img src="image/floating_object_2d_density_2000_cmap3.png" width=45%>
</br>
<img src="image/floating_object_2d_density_2000_cmap4.png" width=45%>
<img src="image/floating_object_2d_density_2000_cmap5.png" width=45%>
</br>
fig. floating object density = 2.0 fluid density
</center>

## Conclusion on 2D case

As can bee seen from the figures, the floating object will sink into the fluid when its density is larger than the fluid. The floating object will float on the fluid when its density is smaller than the fluid. The floating height of the object is determined by the balance of pressure and gravity. Density ratio determines the depth of the floating object in the fluid.