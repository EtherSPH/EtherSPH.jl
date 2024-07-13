[toc]

# Rigid Body

Rigid body is an ideal model of a solid body in which deformation is neglected. In other words, the distance between any two given points on a rigid body remains constant in time regardless of external forces exerted on it. Sound speed in a rigid body is $c=\infty$. The motion of a rigid body can be described by the motion of its center of mass and the rotation of the body about its center of mass.

By the law of center of mass, the mass center of a rigid body can be calculated as:

$$
\begin{equation}
    \begin{aligned}
        \vec{r}_{c} = \frac{\int\vec{r}\mathrm{d}m}{\int\mathrm{d}m}
        =\frac{\sum_j m_j\vec{r}_j}{\sum_j m_j}
    \end{aligned}
\end{equation}
$$

As the body discretizes into $n$ particles, by Newton's second law, the acceleration of the mass center is:

$$
\begin{equation}
    \begin{aligned}
        m_c\vec{a}_c &=\sum_j m_j\vec{a}_j\\
        I_c\cdot\vec{\alpha} &=\sum_j m_j(\vec{r}_j-\vec{r}_c)\times\vec{a}_j
    \end{aligned}
\end{equation}
$$

$I_c$ is inertia tensor of the rigid body, $\vec{\alpha}$ is angular acceleration, $\vec{a}_j$ is acceleration of particle $j$. In Einstein notation, $I_c$ can be written as (in cartesian coordinates instead of curved coordinates):

$$
\begin{equation}
    \begin{aligned}
        I_c^{ik} = \int_{R.B.}(
            \delta^{ik}r^m r^m - r^i r^k
        )\mathrm{d}m
    \end{aligned}
\end{equation}
$$

Currently (2024.07.13) I just finish a simple rigid body simulation in 2D, without constraints or multiple-body dynamics:

```julia
@kwdef mutable struct SingleRigidBody2D{ParticleType <: AbstractParticle2D}
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    v_vec_::Vector2D = Vector2D(0.0, 0.0)
    dv_vec_::Vector2D = Vector2D(0.0, 0.0)
    omega_::Float64 = 0.0
    alpha_::Float64 = 0.0
    force_vec_::Vector2D = Vector2D(0.0, 0.0)
    torque_::Float64 = 0.0
    mass_::Float64 = 0.0
    inertia_::Float64 = 0.0
    particles_::Vector{ParticleType} = ParticleType[]
end
```

I will add more features in the future maybe with aid of [RigidBodyDynamics.jl](https://github.com/JuliaRobotics/RigidBodyDynamics.jl)

# Examples of rigid body dynamics

**TODO**