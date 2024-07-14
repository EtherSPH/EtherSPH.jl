#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/07/14 20:03:54
  @ license: MIT
  @ description:
 =#

@kwdef mutable struct VariableContainer{VariableType}
    variable_list_::Vector{VariableType} = Vector{VariableType}()
    length_::Int64 = 0
end

@inline function Base.length(variable_container::VariableContainer{VariableType})::Int64 where {VariableType}
    return variable_container.length_
end

@inline function Base.getindex(
    variable_container::VariableContainer{VariableType},
    i::Int64,
)::VariableType where {VariableType}
    return variable_container.variable_list_[i]
end

@inline function Base.getindex(
    variable_container::VariableContainer{VariableType},
    range::UnitRange{Int64},
)::Vector{VariableType} where {VariableType}
    return variable_container.variable_list_[range]
end

@inline function Base.setindex!(
    variable_container::VariableContainer{VariableType},
    value::VariableType,
    i::Int64,
)::Nothing where {VariableType}
    variable_container.variable_list_[i] = value
    return nothing
end

@inline function Base.eachindex(
    variable_container::VariableContainer{VariableType},
)::UnitRange{Int64} where {VariableType}
    return Base.OneTo(length(variable_container))
end

@inline function Base.show(io::IO, variable_container::VariableContainer{VariableType})::Nothing where {VariableType}
    print(io, "VariableContainer{", VariableType, "}(")
    for i in 1:length(variable_container)
        print(io, variable_container[i])
        if i < length(variable_container)
            print(io, ", ")
        end
    end
    print(io, ")\n")
    return nothing
end

@inline function capacity(variable_container::VariableContainer{VariableType})::Int64 where {VariableType}
    return length(variable_container.variable_list_)
end

@inline function setlength!(
    variable_container::VariableContainer{VariableType},
    length::Int64,
)::Nothing where {VariableType}
    variable_container.length_ = length
    return nothing
end

@inline function reset!(variable_container::VariableContainer{VariableType})::Nothing where {VariableType}
    setlength!(variable_container, 0)
    return nothing
end

@inline function Base.resize!(
    variable_container::VariableContainer{VariableType},
    new_capacity::Int64,
)::Nothing where {VariableType}
    if new_capacity > capacity(variable_container)
        resize!(variable_container.variable_list_, new_capacity)
    end
    return nothing
end

@inline function Base.push!(
    variable_container::VariableContainer{VariableType},
    value::VariableType,
)::Nothing where {VariableType}
    if length(variable_container) == capacity(variable_container)
        resize!(variable_container, capacityExpandPolicy(capacity(variable_container)))
    end
    setlength!(variable_container, length(variable_container) + 1)
    setindex!(variable_container, value, length(variable_container))
    return nothing
end

const RealNumberContainer = VariableContainer{Float64}

const RealVectorContainer{Dimension} = VariableContainer{RealVector{Dimension}}
const Vector2DContainer = RealVectorContainer{Dimension2}
const Vector3DContainer = RealVectorContainer{Dimension3}

@inline function Base.resize!(
    real_vector_container::RealVectorContainer{Dimension},
    new_capacity::Int64,
)::Nothing where {Dimension}
    if new_capacity > capacity(real_vector_container)
        resize!(real_vector_container.variable_list_, new_capacity)
        for i in (length(real_vector_container) + 1):new_capacity
            real_vector_container.variable_list_[i] = Vector0(Dimension)
        end
    end
    return nothing
end

@inline function Base.push!(
    real_vector_container::RealVectorContainer{Dimension},
    value::RealVector{Dimension},
)::Nothing where {Dimension}
    if length(real_vector_container) == capacity(real_vector_container)
        resize!(real_vector_container, capacityExpandPolicy(capacity(real_vector_container)))
    end
    setlength!(real_vector_container, length(real_vector_container) + 1)
    # setindex!(real_vector_container, value, length(real_vector_container))
    real_vector_container.variable_list_[length(real_vector_container)] .= value
    return nothing
end
