using StaticArrays: MMatrix

export Lattice, basisvectors

"""
    AbstractLattice{T} <: AbstractMatrix{T}

Represent the real lattices and the reciprocal lattices.
"""
abstract type AbstractLattice{T} <: AbstractMatrix{T} end
struct Lattice{T} <: AbstractLattice{T}
    data::MMatrix{3,3,T,9}
end
"""
    Lattice(data::AbstractMatrix)

Construct a `Lattice` from a matrix.

!!! note
    The basis vectors of the matrix are stored as columns.
"""
Lattice(data::AbstractMatrix) = Lattice(MMatrix{3,3}(data))
"""
    Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Construct a `Lattice` from three basis vectors.
"""
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
"""
    Lattice(basisvectors::AbstractVector{<:AbstractVector})

Construct a `Lattice` from a vector of three basis vectors.
"""
Lattice(basisvectors::AbstractVector{<:AbstractVector}) =
    Lattice(reduce(hcat, basisvectors))

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
