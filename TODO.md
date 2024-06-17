# TODO LIST

1. Shape
    - [ ] add more standard shape types to `Geometry`
    - [ ] generation of particles for standard shape
    - [ ] generation of particles for arbitary shape (via `point-cloud-utils` using `PyCall`)
    - [ ] ...
2. Boundary
    - [x] compulsive wall force
    - [ ] periodic boundary in each dimension
    - [ ] wall boundary
    - [ ] buffer particle zone boundary
    - [ ] ...
3. Model
    - [ ] linear elastic equation model
    - [ ] arbitray constitutive relation (via `Einsum.jl` fpor tensor operation)
    - [ ] fsi (fluid-solid interaction)
    - [ ] rigid body dynamics
    - [ ] thermal transfer model
    - [ ] ...
4. Trick
    - [x] kernel average density filter
    - [ ] XSPH
    - [ ] ...
5. Example
    - [ ] poisson equation
        - [x] 1D
        - [x] 2D
    - [x] collapse dry
    - [ ] lid-driven cavity
    - [ ] poiseuille flow
    - [ ] flow over cylinder
    - [ ] natural convection cavity
    - [ ] natural convection circular
    - [ ] ...
6. Others
    - [ ] `DEM` method
    - [ ] documentation github pages
    - [ ] neural network acceleration
    - [ ] ...