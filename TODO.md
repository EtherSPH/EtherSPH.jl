# TODO LIST

1. Shape
    - [ ] add more standard shape types to `Geometry`
        - [x] Rectangle
        - [x] Ring
        - [x] Circle
        - [x] RingColumn
        - [x] Cylinder
        - [ ] ...
    - [ ] generation of particles for standard shape
    - [ ] generation of particles for arbitary shape (via `point-cloud-utils` using `PyCall`)
    - [ ] ...
2. Boundary
    - [x] compulsive wall force
    - [ ] periodic boundary in each dimension
        - [x] 1 dimension
        - [x] 2 dimension
        - [ ] 3 dimension
    - [ ] wall boundary
    - [x] buffer particle zone boundary
    - [ ] ...
3. Model
    - [x] linear elastic equation model (see hanging beam example)
    - [x] large deformation constitutive relation (see `SPH/Cart/LargeDeformationContinuumMechanics.jl`)
    - [ ] arbitray constitutive relation
    - [ ] fsi (fluid-solid interaction)
        - [x] fluid-rigid body interaction
        - [x] fluid-elastic body interaction
        - [ ] ...
    - [ ] rigid body dynamics
        - [x] 2D rigid body dynamics
        - [ ] 3D rigid body dynamics
    - [x] thermal transfer model
    - [ ] ...
4. Trick
    - [x] kernel average density filter
    - [x] $\delta$ SPH
    - [ ] XSPH
    - [ ] ...
5. Example
    - [x] poisson equation
    - [x] single wave
    - [x] hydrostatic pressure
    - [x] collapse dry
    - [x] cruchaga
    - [x] lid-driven cavity
    - [x] poiseuille flow
    - [x] flow over cylinder
    - [x] natural convection cavity
    - [x] natural convection ring
    - [x] hanging beam
    - [x] floating object
    - [ ] ...
6. Others
    - [ ] `DEM` method
    - [ ] documentation github pages
    - [ ] neural network acceleration
    - [ ] SPH interpolation to grid
    - [ ] particle id reordering
    - [ ] ...