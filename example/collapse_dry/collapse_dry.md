[toc]

# Collapse Dry Problem

This case is from [Fluid mechanics and the SPH method: theory and applications](https://academic.oup.com/book/8904), on page 484, by Damien Violeau. In other papers, the case can also be called `Dam break problem`.

<center>
<image src="image/collapse_dry.png">
</br>
<b>fig. collapse onto a dry bottom</b>
</center>

It's also included in [SmoothedParticles.jl](https://github.com/OndrejKincl/SmoothedParticles.jl).

# Different Equation Model Comparison

## Same: treat wall particles as fixed water particles

This method treats wall particles as fixed water particles. They have continuity & momentum interaction with normal fluid particles.

<center>
<image src="image/collapse_dry_same_reference.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - treat wall as the same</b>
</center>

<center>
<image src="image/collapse_dry_same_cmap.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - treat wall as the same</b>
</center>

> note 1: Actualluy, i find this method will meet numerical unstability when `dr` gets smaller. Setting `dr` as 0.01 will more easily meet negative pressure problem. -- 2024.06.15

> note 2: Wall particles provide a pressure foce to avoid fluid particles penetrating them. This will create a blank gap between fluid particles and wall particles, which is not the case in real life. -- 2024.06.15

> note 3: a few water particles will hang on the wall. -- 2024.06.15

## Compulsive: a compulsive force prevents fluid particles from penetrating

This method is from [SPH MODELING OF TSUNAMI WAVES, Rogers & Dalrymple - 2008](http://www.worldscientific.com/doi/abs/10.1142/9789812790910_0003). A compulsive force is applied on fluid particles by wall particles.

<center>
<image src="image/collapse_dry_compulsive_reference.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - compulsive wall</b>
</center>

<center>
<image src="image/collapse_dry_compulsive_cmap.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - compulsive wall</b>
</center>

> note 1: WendlandC2 kernel performs better stability than CubicSpline. -- 2024.06.15

> note 2: Still there will be a 'blank gap'. A thin layer of water particles close to wall boundary performs the pressure force on the other water particles, introducing a 'blank layer' between wall and fluid. -- 2024.06.15

## Pressure Extrapolation: extrapolate pressure from fluid particles to wall particles

The basic idea of pressure extrapolation is to calculate the pressure of wall particles by extrapolating the pressure of fluid particles:

$$
\begin{equation}
    \begin{aligned}
    p_w=\sum_f \frac{\frac{m_f}{\rho_f}W_{pf}(p_f+\rho_f\vec{g}\cdot\vec{r}_{wf})}{\frac{m_f}{\rho_f}W_{pf}}
    \end{aligned}
\end{equation}
$$

However, when negative pressure or wall particles locating upon fluid particles, equation above contributes negative pressure on wall pressure, causing abnormal absorption of fluid particles. A cutoff is added to prevent negative pressure:

$$
\begin{equation}
    \begin{aligned}
        p_f^* = \max(p_f,0.0)
        \quad
        \rho_f\vec{g}\cdot\vec{r}_{wf}^* = \max(\rho_f\vec{g}\cdot\vec{r}_{wf},0.0)
    \end{aligned}
\end{equation}
$$

At the same time, a balanced continuity and pressure force form should be applied:

$$
\begin{equation}
    \begin{aligned}
        \left(\frac{\mathrm{d}\rho}{\mathrm{d}t}\right)_i &= 
        \sum_j\rho_i\frac{m_j}{\rho_j}\vec{v}_{ij}\cdot \nabla W_{ij}\\
        \vec{f}_{i}^p &=
        \sum_j -m_j \frac{p_i+p_j}{\rho_i \rho_j}\nabla W_{ij}
    \end{aligned}
\end{equation}
$$

<center>
<image src="image/collapse_dry_extrapolation_reference.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - pressure extrapolation</b>
</center>

<center>
<image src="image/collapse_dry_extrapolation_cmap.png" width=70%>
</br>
<b>fig. collapse onto a dry bottom - pressure extrapolation</b>
</center>

> *what's good?*: No blank gap between wall and fluid particles. Fit for all kind of wall material (any value of density is allowed). No need to generate normal vector.
> *what's bad?*: a little bit unstable than the other 2 methods mentioned above. -- 2024.07.11 