# EtherSPH.jl

Maybe the last version in julia, largely thanks to [SmoothedParticles.jl](https://github.com/OndrejKincl/SmoothedParticles.jl).

`EtherSPH.jl` is a Smoothed Particle Hydrodynamics (SPH) code written in Julia. It is designed to be a general-purpose SPH code that can be used to simulate a wide range of problems in fluid dynamics, solid mechanics, and fluid-structure interaction. The code is designed to be easy to use and easy to extend, with a focus on readability and simplicity. As ether was once considered to be a particle-like medium that permeated all of space which was proved to be non-existent, the name `EtherSPH` is chosen to reflect the particle-based nature of the simulation method.

# Setup

In julia repl, type:

```julia
julia> ']' # enter package mode
pkg> activate .
pkg> instantiate # or resolve maybe needed
pkg> test
```

Also, to plot results and do post-processing, a python environment with packages in [requirements.txt](requirements.txt) is needed. Either `pip` or `conda` can be used to install the packages:

```bash
pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple # use tsinghua mirror if needed
```

# Examples

1. classical pde problems
    - [poisson equation](example/poisson_equation/poisson_equation.md)
    - [single wave](example/single_wave/single_wave.md)
    - ...
2. free surface problems
    - [collapse dry](example/collapse_dry/collapse_dry.md)
    - [cruchaga](example/cruchaga/cruchaga.md)
    - [hydrostatic pressure](example/hydrostatic_pressure/hydrostatic_pressure.md)
    - ...
3. classical fluid dynamics problems
    - [lid-driven cavity](example/lid_driven_cavity/lid_driven_cavity.md)
    - [poiseuille flow](example/poiseuille_flow/poiseuille_flow.md)
    - [cylinder](example/cylinder/cylinder.md)
    - ...
4. heat transfer problems
    - [natural convection cavity](example/natural_convection_cavity/natural_convection_cavity.md)
    - [natural convection ring](example/natural_convection_ring/natural_convection_ring.md)
    - ...
5. solid mechanics problems
    - [hanging beam](example/hanging_beam/hanging_beam.md)
    - ...
6. rigid body problems
    - [rigid body](example/rigid_body/rigid_body.md)
    - ...
7. fsi (fluid-structure interaction) problems
    - [floating object](example/floating_object/floating_object.md)
    - [fsi: cylinder with plate](example/fsi_cylinder_with_plate/fsi_cylinder_with_plate.md)
    - [fsi: dam break onto obstacle](example/fsi_dam_break_onto_obstacle/fsi_dam_break_onto_obstacle.md)
    - ...