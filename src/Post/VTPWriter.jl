#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
  @ date: 2024/06/12 21:39:08
  @ license: MIT
  @ description:
 =#

const kFileExtension = ".vtp"
const kWallTimeFormat = "yyyy_mm_dd_HH_MM_SS.SS"
const kDensityString = "Density"
const kMassString = "Mass"
const kTypeString = "Type"
const kVelocityString = "Velocity"

@inline function assurePathExist(path::String)::Nothing
    if !isdir(path)
        @info "Create directory: $path"
        mkpath(path)
    else
        @info "Remove all files in directory: $path"
        for file in readdir(path)
            rm(joinpath(path, file))
        end
    end
    return nothing
end

@inline function getWallTime()::String
    return Dates.format(now(), kWallTimeFormat)
end

@kwdef mutable struct VTPWriter
    output_count_::Int64 = 0
    step_digit_::Int64 = 4 # 0001, 0002, ..., 9999
    file_name_::String = "result"
    output_path_::String = "example/results"
    scalar_name_list_::Vector{String} = String[kMassString, kDensityString]
    getScalarFunctions_::Vector{Function} = Function[getMass, getDensity]
    vector_name_list_::Vector{String} = String[]
    getVectorFunctions_::Vector{Function} = Function[]
end

Base.show(io::IO, vtp_writer::VTPWriter) = print(
    io,
    "VTPWriter:\n",
    "    step digit: $(vtp_writer.step_digit_)\n",
    "    file name: $(vtp_writer.file_name_)\n",
    "    output path: $(vtp_writer.output_path_)\n",
    "    scalar name list: $(vtp_writer.scalar_name_list_)\n",
    "    vector name list: $(vtp_writer.vector_name_list_)\n",
    "    current writen time: $(vtp_writer.output_count_)\n",
)

@inline assurePathExist(vtp_writer::VTPWriter)::Nothing = assurePathExist(vtp_writer.output_path_)

@inline function getOutputFileName(vtp_writer::VTPWriter)::String
    return joinpath(
        vtp_writer.output_path_,
        string(vtp_writer.file_name_, string(vtp_writer.output_count_, pad = vtp_writer.step_digit_), kFileExtension),
    )
end

@inline function addScalar!(vtp_writer::VTPWriter, scalar_name::String, scalarFunction::Function)::Nothing
    push!(vtp_writer.scalar_name_list_, scalar_name)
    push!(vtp_writer.getScalarFunctions_, scalarFunction)
    return nothing
end

@inline function addVector!(vtp_writer::VTPWriter, vector_name::String, vectorFunction::Function)::Nothing
    push!(vtp_writer.vector_name_list_, vector_name)
    push!(vtp_writer.getVectorFunctions_, vectorFunction)
    return nothing
end

@inline function saveVTP(
    vtp_writer::VTPWriter,
    particle_system::ParticleSystem{Dimension, ParticleType},
    step::Int64,
    simulation_time::Float64,
)::Nothing where {Dimension, ParticleType <: AbstractParticle{Dimension}}
    n_particles = length(particle_system)
    type = zeros(Int64, n_particles)
    positions = zeros(Float64, Dimension, n_particles)
    scalars_list = [zeros(Float64, n_particles) for _ in 1:length(vtp_writer.scalar_name_list_)]
    vectors_list = [zeros(Float64, Dimension, n_particles) for _ in 1:length(vtp_writer.vector_name_list_)]
    Threads.@threads for i in eachindex(particle_system)
        @inbounds type[i] = getType(particle_system[i])
        @inbounds positions[:, i] = getPosition(particle_system[i])
        for j in 1:length(vtp_writer.scalar_name_list_)
            @inbounds scalars_list[j][i] = vtp_writer.getScalarFunctions_[j](particle_system[i])
        end
        for j in 1:length(vtp_writer.vector_name_list_)
            @inbounds vectors_list[j][:, i] .= vtp_writer.getVectorFunctions_[j](particle_system[i])
        end
    end
    file_name = getOutputFileName(vtp_writer)
    cells = [MeshCell(PolyData.Verts(), [i]) for i in 1:n_particles]
    vtp_file = vtk_grid(file_name, positions, cells)
    vtp_file["TMSTEP"] = step
    vtp_file["TimeValue"] = simulation_time
    vtp_file["WallTime"] = getWallTime()
    vtp_file["Type"] = type
    for i in eachindex(vtp_writer.scalar_name_list_)
        @inbounds vtp_file[vtp_writer.scalar_name_list_[i]] = scalars_list[i]
    end
    for i in eachindex(vtp_writer.vector_name_list_)
        @inbounds vtp_file[vtp_writer.vector_name_list_[i]] = vectors_list[i]
    end
    vtk_save(vtp_file)
    vtp_writer.output_count_ += 1
    return nothing
end
