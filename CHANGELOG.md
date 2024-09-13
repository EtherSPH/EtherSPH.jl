[toc]

---

# Changelog

We start maintain [Changelog](CHANGELOG.md) from **2024-07-14** for important changes in the project. [TODO](TODO.md) has already in maintenance while it's not capable to record the changes in time sequence.

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [0.1.0] - 2024.09.13

### Brief

#### FSI: Dambreak onto Obstacle

A fluid-structure interaction problem is solved. See [fsi dambreak onto obstacle](example/fsi_dam_break_onto_obstacle/fsi_dam_break_onto_obstacle.md) for details.

#### `VTPWriter` for `binary` choice

As `vtk` file can be either stored in `ascii` or `binary` format, `VTPWriter` is updated to support `use_binary` choice when calling `splitParticlesType`.

#### New Author

Welcome [`MaiZLnuaa`](mai-zl@nuaa.edu.cn) to join the project!

#### $\delta$-SPH

$\delta$-SPH is introduced to filter the density in SPH method. See `SPH/Cart/DeltaSPH.jl` for details. And a demo in [collapse dry](example/collapse_dry/collapse_dry_extrapolation_delta_sph.jl) is included.

#### `applyInteractionByKernelGradient!`: An important function in `src/Core/Apply.jl`

`applyInteractionByKernelGradient!` is introduced to apply interaction by kernel gradient. As the kernel gradient $\nabla W_{ij}$ can be called several times in the interaction, it's recommended to calculate the kernel gradient first and store them. Thus, `applyInteractionByKernelGradient!` is introduced to parse `kernel_gradient` into the interaction function. You may see the change in file [Apply.jl](src/Core/Apply.jl).

### Added

- `Post/VTPWriter.jl` adds `use_binary` choice in `splitParticlesType` function.
- `SPH/Cart/DeltaSPH.jl` adds $\delta$-SPH to filter the density in SPH method.
- `example/fsi_dam_break_onto_obstacle` adds fluid-structure interaction problem: dambreak onto obstacle.
- `example/collapse_dry` adds `collapse_dry_extrapolation_delta_sph.jl` to show how to use $\delta$-SPH in SPH method.
- `demo/dam_break_2d/dam_break_2d_delta_sph.jl` is added to show how to use $\delta$-SPH in SPH method.
- `Core/Apply.jl` adds `applyInteractionByKernelGradient!` function to apply interaction by kernel gradient.

### Changed

- Nothing

### Fixed

- Nothing

### Removed

- Nothing

---

## [0.1.0] - 2024.08.27

### Brief

#### FSI: Cylinder with Plate

A fluid-structure interaction problem is solved. The cylinder is fixed and the plate is free to move. It's a classic problem in fluid-structure interaction. The fluid is incompressible and the flow is laminar. The Reynolds number varies in $20, 100, 200$, marked as fsi1, fsi2 and fsi3.

See [fsi_cylinder_with_plate](example/fsi_cylinder_with_plate/fsi_cylinder_with_plate.md) for details.

#### `Project.toml`

Add `FLoops` and `ProgressMeter` package for adaptive time step in simulations.

`FLoop.jl` package is used to collect the time step information in the simulation. `ProgressMeter.jl` package is used to show the progress with unknown total steps.

#### `collectmin` and `collectmax`

`collectmin` and `collectmax` are added in `src/Core`. They are used to collect the minimum and maximum value particle system.

#### Some functions for `VTPWriter`

`VTPWriter` is used to write the vtk file. Some functions are added to write the field data in the vtk file.

### Added

- `FLoops` and `ProgressMeter` in `Project.toml`.
- `Core/Collect.jl` adds `collectmin` and `collectmax` functions to collect the minimum and maximum value particle system.
- `example/fsi_cylinder_with_plate` adds fluid-structure interaction problem: cylinder with plate.
- `Post/VTPWriter.jl` adds `setStepDigit!`, `setFileName!` and `setOutputPath!` functions to modify the output file name and path.

### Changed

- Nothing

### Fixed

- Nothing

### Removed

- Nothing

---

## [0.1.0] - 2024.08.01

### Brief

#### Example: Flow Over Cylinder

Use pressure extrapolation method to handle wall boundary, and add a new boundary: **buffer + periodic** which fits the reference well.

#### Pressure Force

In `SPH/Library/PressureForce.jl`, a density weighted pressure force is added which implies: $p^* = (p_i\rho_i + p_j\rho_j)/(\rho_i + \rho_j)$. Also, pressure extrapolation for wall boundary is added. 

### Added

- `SPH/Library/PressureForce.jl` adds `libDensityWeightedPressureForce!` function to calculate density weighted pressure force.
- `SPH/Library/PressureForce.jl` adds `libPressureExtrapolationInteraction!` function to apply pressure extrapolation for wall boundary.
- `SPH/Library/PressureForce.jl` adds `libPressureExtrapolationSelfaction!` function to apply pressure extrapolation for wall boundary.
- `example/cylinder` adds some case and updates `cylinder.md`.

### Changed

- Nothing

### Fixed

- Nothing

### Removed

- Nothing

---

## [0.1.0] - 2024.07.26

### Brief

#### Split Particles by Type

Sometimes post process requires to split particles by type. Using `PyCall`, I import [`pyvsita`](https://docs.pyvista.org/) to read `.vtp` and split particles by type into different vtk files. The code is shown as below:

```julia
splitParticlesType(
    vtp_writer;
    type_name_dict = Dict(
        1 => "Water", 
        2 => "Wall", 
        3 => "Beam_Move", 
        4 => "Beam_Fixed"),
)
```

When `type_name_dict` is not given, the default type name is `Particle_1_xxx.vtk`, `Particle_2_xxx.vtk`, etc. However, such function requires user to add `PyCall` in their python environment which is correctly configured with julia.

#### Marching towards FSI

FSI problem is interesting and can not be solved easily by traditional numerical methods. I start to march towards FSI problem. A simple FSI problem from [An SPH‑based FSI framework for phase‑feld modeling of brittle fracture under extreme hydrodynamic events](https://link.springer.com/10.1007/s00366-023-01857-0) is included in demo and waits for further validation. Although I should not publish this version since it's not been validated, the `splitParticlesType` function matters more for other users, which makes me decide to push this version into the main branch.

### Added

- `Post/VTPWriter.jl/splitParticlesType` function to split particles by type into different vtk files.
- `demo/fsi_water_onto_obstacle` is added.

### Changed

- Nothing

### Fixed

- Nothing

### Removed

- Nothing

---

## [0.1.0] - 2024-07-23

### Brief

#### More periodic boundary conditions

In the previous version, we introduced periodic boundary conditions in one `x` direction for 2D and 3D cases. In this version, we add periodic boundary conditions in both `x` and `y` directions for 2D cases. And single periodic direction can be added for `y` or `z` direction in both 2D and 3D cases.

#### Model pde: single wave

My tutor once required me to solve the propogation of sound using SPH method. Although it's of none sense for SPH focues on Lagrangian view, I still try to solve the single wave propogation in 2D case. See `example/single_wave/single_wave_2d.jl`. The results shows that particle discretization is accurate in some aspect.

### Added

- Periodic boundary condition in one direction at each dimension. See `Core/PeriodicBoundary/PeriodicBoundary.jl` for details.
- Periodic boundary condition in both `x` and `y` direction in 2D case. See `Core/PeriodicBoundary/PeriodicBoundaryXY.jl` for details.
- `example/single_wave` folder for single wave propogation in 2D case.

### Changed

- move `Core/PeriodicBoundary.jl` to `Core/PeriodicBoundary/PeriodicBoundary.jl`

### Fixed

- Nothing

### Removed

- Nothing

---

## [0.1.0] - 2024-07-20

### Brief

#### Call `ConstVec` or `ConstMat`

In continuum mechanics, `ConstVec` or `ConstMat` is recommended to use for constant vectors or matrices. We need to use a vector or matrix with different dimensions in calculation.

How it works? Let's see a simple example:

```julia
const MatrixI2D = MatrixI(2)
const MatrixI3D = MatrixI(3)

struct ConstMatI{Dimension} end
@inline function ConstMatI{Dimension2}()::Matrix2D
    return MatrixI2D
end
@inline function ConstMatI{Dimension3}()::Matrix3D
    return MatrixI3D
end
```

The struct `ConstMatI` is defined to call `MatrixI2D` or `MatrixI3D` according to the dimension option. We can use a constant matrix $\mathrm{I}$ like:

```julia
A = Matrix2D(0., 0., 0., 0.)
A .= ConstMatI{2}() # without any cost
```

Supported types:

- `ConstVec0{Dimension}`
- `ConstVecX{Dimension}`
- `ConstVecY{Dimension}`
- `ConstVecZ{Dimension}`
- `ConstMatI0Dimension}`
- `ConstMatI{Dimension}`

#### `TinyLinearAlgebra` for GC problems

In our tests, `MArray` from `StaticArrays` suffers from GC problems. A small example is that:

```julia
x = Vector2D(1., 2.)
A = Matrix2D(1., 2., 3., 4.)
c = Vector2D(0., 0.)
B = Matrix2D(0., 0., 0., 0.)

c .= A * x # suffers from GC problems in for-loop
c .= x' * A # suffers from GC problems in for-loop
c .= dot(A, x) # recommended
c .= dot(x, A) # recommended

B .= x * x' # suffers from GC problems in for-loop
B .= dyad(x, x) # recommended

B .= A * A # recommended, the same as below ↓
B .= A * A' # recommended
B .= dot(A, A) # recommended
B .= dot(A, A') # recommended
```
It seems that gc problems are caused by converting `Vector` to `Matrix` or `Matrix` to `Vector`. Always carefully check the GC problems in your code. If you find abnormal speed down in your code, it's recommended use `profview` to check the GC problems. Functions in `Math` are recommended to use `TinyLinearAlgebra` to avoid GC problems.

#### `VTPWriter` for field data

Sometimes, additional field data in one vtp file is needed. We introduce an optional argument `field_dict::Dict{String, Any} = Dict{String, Any}()` in `saveVTP` function to save field data in VTP file. Construct a `Dict` and parse it into `saveVTP` function to save field data in VTP file.

#### `NeighbourSearch` and `Apply` for continuum mechanics

In continuum mechanics, time-varied neighbour search is not needed. We should record all neighbour index first (ensure the index not changes in your particle system). Then, some function should be applied on particles in each time step.

For the order of neighbour particle indexes contained by one particle matters, the particle action should be applied with the sequence given. `directlyUpdateNeighbours!` is introduced to update neighbour particles in continuum mechanics. `applyReflection!` is introduced to apply reflection particle action if related particle ids are recorded. The code is simply like:

```julia
@inline function directlyUpdateNeighbours!(
    particle_system::ParticleSystem{Dimension, ParticleType},
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    Threads.@threads for i in eachindex(particle_system)
        p = particle_system.particles_[i]
        @simd for j in eachindex(p.neighbour_index_list_)
            @inbounds q_index = p.neighbour_index_list_[j]
            @inbounds p.neighbour_position_list_[j] .= p.x_vec_ .- particle_system[q_index]
            @inbounds p.neighbour_distance_list_[j] = norm(p.neighbour_position_list_[j])
        end
    end
    return nothing
end

@inline function applyReflection!(
    particle_system::ParticleSystem{Dimension, ParticleType},
    reflectionFunction!::Function;
    parameters...,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    Threads.@threads for i in eachindex(particle_system)
        @inbounds reflectionFunction!(particle_system, particle_system[i]; parameters...)
    end
    return nothing
end
```

#### Continuum Mechanics Model

In `SPH/Cart/LargeDeformationContinuumMechanics.jl`, a large deformation continuum mechanics model is introduced. Any elastic material can now be calculated with the model. The model is based on the continuum mechanics theory and the SPH discretization. The calculation steps is briefly shown as below:

- initialize the neighbour particles and kernel gradient.
- calculate the deformation gradient and green-lagrange strain.
- use constitutive relation (either elastic or plastic) to calculate stress.
- calculate the acceleration of particles with its stress (or artificial stress should be added)

Examples can be found [here: haning beam continuum](example/hanging_beam/hanging_beam_continuum_2d.jl). A [demo](demo/deformation_2d/deformation_2d.jl) when I was developing such model is also included.

To achieve the model, `SPH/Library` folder adds:

- `SPH/Library/ViscosityForce.jl` adds `libArtificialViscosityForce!` function to apply artificial stress to avoid particles being too close to each other.
- `SPH/Library/Correction.jl` is added with corrected matrix calculated. In paper, such matrix is usually denoted as $\mathbf{B}$ or $\mathbf{L}^{-1}$.

### Added

- `Math/Vector.jl` adds `ConstVec` xxx.
- `Math/Matrix.jl` adds `ConstMat` xxx.
- `Math/TinyLinearAlgebra.jl` adds `dot` for `Matrix` and `Vector`.
- `Post/VTPWriter.jl/saveVTP` function adds an optional argument .`field_dict::Dict{String, Any} = Dict{String, Any}()` to save field data in VTP file.
- `Core/NeighbourSearch.jl` adds `directlyUpdateNeighbours!` for continuum mechanics.
- `Core/Apply.jl` adds `applyReflection!` for users to apply reflection
particle action if related particle ids are recorded.
- `SPH/Cart/Cart.jl` is added which is just like a cart in shopping mall. Feel free to add model described in papers here!
- `SPH/Cart/LargeDeformationContinuumMechanics.jl` adds a large deformation continuum mechanics model.
- `SPH/Library/ViscosityForce.jl` adds `libArtificialViscosityForce!` function to apply artificial stress.
- `SPH/Library/Correction.jl` adds corrected matrix calculation.
- `demo/deformation_2d/deformation_2d.jl` is added to show how I tested and developed the large deformation continuum mechanics model.
- `example/hanging_beam/hanging_beam_continuum_2d.jl` is added to show how to use the large deformation continuum mechanics model contained in `src/SPH/Cart/LargeDeformationContinuumMechanics.jl`.

### Changed

- `Math/Matrix.jl/dot` is changed to `Math/Matrix.jl/ddot` to avoid conflict.

### Fixed

- Nothing

### Removed

- Nothing

---

## [0.1.0] - 2024-07-14

### Brief

This is the first version of Changelog. We start maintain Changelog from this version.

For interaction between particles could be applied for over one time in each time step, original `applyInteraction!` will do neighbour search in cells over and over again, which is not efficient. 

Original `applyInteraction!` usage looks like:

```julia
for step in 1: total_step
    # ...
    applyInteraction!(system, some_interaction!)
    applySelfaction!(system, some_selfaction!)
    # ...
    # cell contains particle indexes ↓
    # instead of particles containing their neighbours' indexes ↓
    createCellLinkList!(system)
    # ...
end
```

As long as complicated or multiple-time interaction in one time step is applied, the efficiency of `applyInteraction!` will be decreased for one particle will first locate itself in one cell, and then search its neighbours in cell's neighbourhood.

From this version on, a more efficient way is introduced. You may see a [demo](demo/dam_break_2d/dam_break_2d_neighbour_list.jl). We require user define a particle type as follows:

```julia
@kwdef mutable struct Particle <: AbstractParticle2D
    # must have
    x_vec_::Vector2D = Vector2D(0.0, 0.0)
    rho_::Float64 = rho_0
    mass_::Float64 = mass
    type_::Int64 = FLUID_TAG
    # neighbour information field is now required
    neighbour_index_list_::IndexContainer = IndexContainer()
    neighbour_position_list_::Vector2DContainer = Vector2DContainer() # to store `rpq`
    neighbour_distance_list_::RealNumberContainer = RealNumberContainer() # to store `r`
    # the property you defined
    # ...
end
```

And the usage of `applyInteraction!` is changed to:

```julia
for step in 1: total_step
    # ...
    applyInteractionWithNeighbours!(system, some_interaction!)
    applySelfaction!(system, some_selfaction!)
    createNeighbourIndexList!(system)
    # ...
end
```

Original `applyInteraction` is still kept for compatibility. `applyInteractionWithNeighbours!` will act function with the neighbour particles recorded in `neighbour_index_list_`. At the same time, `createCellLinkList!` should be replaced by `createNeighbourIndexList!` to record the neighbour particles.

Users can still adopt `applyInteraction! + createCellLinkList!` for simple interaction cases. But for complicated or multiple-time interaction in one time step, `applyInteractionWithNeighbours! + createNeighbourIndexList!` is recommended.

### Added

- `VariableContainer{VariableType}`, a variable-length container for any type `VariableType` is introduced.
- `RealNumberContainer = VariableContainer{Float64}` is introduced to act as `neighbour_distance_list_`.
- `RealVectorContainer{Dimension} = VariableContainer{RealVector{Dimension}}` is introduced to act as `neighbour_position_list_`.
- `Math/TinyLinearAlgebra.jl` to perform `dyad` in continuum mechanics is introduced.
- `Matrix`'s `det` is introduced to calculate determinant of a matrix.
- `Core/NeighbourSearch.jl` is introduced to perform neighbour search in cells. And the particles' neighbour particles' indexes can be added into `neighbour_index_list_` in user-defined particle type through `createNeighbourIndexList!`.
- `createNeighbourIndexList!` is introduced to record the neighbour particles' indexes in `neighbour_index_list_`.
- `applyInteractionWithNeighbours!` is introduced to perform interaction.
- [`dam_break_2d_neighbour_list.jl`](demo/dam_break_2d/dam_break_2d_neighbour_list.jl), a demo to show how the noval neighbour search and interaction perform is used.
- CHANGELOG.md.

### Changed

- Move `IndexContainer.jl`, `CartesianIndex.jl` and `VariableContainer.jl` from `Math` folder to `Container` folder.
- Move `createCellLinkList!` function from `Core/ParticleSystem.jl` to `Core/NeighbourSearch.jl`.

### Fixed

- Nothing

### Removed

- Nothing

---