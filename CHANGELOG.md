[toc]

---

# Changelog

We start maintain [Changelog](CHANGELOG.md) from **2024-07-14** for important changes in the project. [TODO](TODO.md) has already in maintenance while it's not capable to record the changes in time sequence.

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [Unreleased] - 2024-07-14

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