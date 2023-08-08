using StaticArrays: MMatrix

export Lattice, basisvectors

"""
    AbstractLattice{T} <: AbstractMatrix{T}

Represent the real lattices and the reciprocal lattices.
"""
abstract type AbstractLattice{T} <: AbstractMatrix{T} end
struct Lattice{T} <: AbstractLattice{T}
    data::MMatrix{3,3,T,9}
    function Lattice{T}(data) where {T}
        @assert !iszero(det3x3(data)) "lattice is singular!"
        return new(data)
    end
end
"""
    Lattice(data::AbstractMatrix)

Construct a `Lattice` from a matrix.

!!! note
    The basis vectors of the matrix are stored as columns.
"""
Lattice(data::AbstractMatrix{T}) where {T} = Lattice{T}(MMatrix{3,3,T}(data))
"""
    Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Construct a `Lattice` from three basis vectors.
"""
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
"""
    Lattice(data)

Construct a `Lattice` from, e.g., a vector of three basis vectors.
"""
function Lattice(data)
    if length(data) == 9  # Works for `NTuple{9}` and generators
        return Lattice(MMatrix{3,3}(data))
    elseif length(data) == 3  # Works for `NTuple{3}` & vector of vectors
        return Lattice(mapreduce(collect, hcat, data))
    else
        throw(ArgumentError("`data` shape is not recognized!"))
    end
end

"""
    basisvectors(lattice::Lattice)

Get the three primitive vectors from a `lattice`.
"""
basisvectors(lattice::Lattice) = lattice[:, 1], lattice[:, 2], lattice[:, 3]

Base.size(::AbstractLattice) = (3, 3)

Base.parent(lattice::AbstractLattice) = lattice.data

Base.getindex(lattice::AbstractLattice, i) = getindex(parent(lattice), i)

Base.setindex!(lattice::AbstractLattice, v, i) = setindex!(parent(lattice), v, i)

Base.IndexStyle(::Type{<:AbstractLattice}) = IndexLinear()

Base.BroadcastStyle(::Type{<:Lattice}) = Broadcast.ArrayStyle{Lattice}()
Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Lattice}}, ::Type{S}
) where {S} = similar(Lattice{S}, axes(bc))
Lattice{S}(::UndefInitializer, dims) where {S} = Lattice(Array{S,length(dims)}(undef, dims))

# `LinearAlgebra.det` is much slower and more inaccurate than my simple `det3x3`.
function det3x3(matrix)
    a, d, g, b, e, h, c, f, i = matrix  # Only works for 3×3 matrices
    return a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h
end
