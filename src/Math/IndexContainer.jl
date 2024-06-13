#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/08 18:15:09
  @ license: MIT
  @ description:
 =#

@inline capacityExpandPolicy(n::Int64)::Int64 = 2 * n + 1

# ! this is not a thread-safe container, work with it coupled with a `Base.Threads.ReentrantLock` in a multi-threaded environment
mutable struct IndexContainer
    index_list_::Vector{Int64}
    length_::Int64 # length is not the same as the length of index_list_
end

@inline function IndexContainer()::IndexContainer
    return IndexContainer(Vector{Int64}(), 0)
end

@inline Base.length(index_container::IndexContainer)::Int64 = index_container.length_
@inline Base.getindex(index_container::IndexContainer, i::Int64)::Int64 = index_container.index_list_[i]
@inline Base.getindex(index_container::IndexContainer, range::UnitRange{Int64})::Vector{Int64} =
    index_container.index_list_[range]
@inline Base.setindex!(index_container::IndexContainer, value::Int64, i::Int64) = index_container.index_list_[i] = value
@inline Base.eachindex(index_container::IndexContainer)::UnitRange{Int64} = Base.OneTo(length(index_container))
@inline function Base.show(io::IO, index_container::IndexContainer)::Nothing
    print(io, "IndexContainer(")
    for i in 1:length(index_container)
        print(io, index_container[i])
        if i < length(index_container)
            print(io, ", ")
        end
    end
    print(io, ")\n")
    return nothing
end
@inline capacity(index_container::IndexContainer)::Int64 = length(index_container.index_list_)
@inline setlength!(index_container::IndexContainer, length::Int64) = index_container.length_ = length

@inline function reset!(index_container::IndexContainer)::Nothing
    # * reset all item in index_list_ to 0, instead of removing them
    @simd for i in eachindex(index_container)
        @inbounds index_container[i] = 0
    end
    setlength!(index_container, 0)
    return nothing
end

@inline function Base.resize!(index_container::IndexContainer, new_capacity::Int64)::Nothing
    # * resize the index_list_ to the given capacity
    if new_capacity > capacity(index_container)
        resize!(index_container.index_list_, new_capacity) # resize doesn't remove exsisiting items, and fill the new space with 0
    end
    return nothing
end

@inline function Base.push!(index_container::IndexContainer, value::Int64)::Nothing
    if length(index_container) == capacity(index_container)
        resize!(index_container, capacityExpandPolicy(capacity(index_container)))
    end
    setlength!(index_container, length(index_container) + 1)
    setindex!(index_container, value, length(index_container))
    return nothing
end
