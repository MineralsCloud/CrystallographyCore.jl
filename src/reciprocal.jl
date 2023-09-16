export ReciprocalLattice

struct ReciprocalLattice{T} <: AbstractLattice{T}
    data::MMatrix{3,3,T,9}
end
ReciprocalLattice(data::AbstractMatrix) = ReciprocalLattice(MMatrix{3,3}(data))

"""
    basisvectors(lattice::ReciprocalLattice)

Get the three basis vectors from a reciprocal lattice.
"""
basisvectors(lattice::ReciprocalLattice) = Tuple(eachbasisvector(lattice))

"""
    eachbasisvector(lattice::ReciprocalLattice)

Iterate over the three basis vectors of a reciprocal lattice.
"""
eachbasisvector(lattice::ReciprocalLattice) = eachcol(lattice)

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L130-L131
Base.one(::Type{ReciprocalLattice{T}}) where {T} =
    ReciprocalLattice(MMatrix{3,3}(SDiagonal(one(T), one(T), one(T))))
Base.one(lattice::ReciprocalLattice) = one(typeof(lattice))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L132-L133
Base.oneunit(::Type{ReciprocalLattice{T}}) where {T} =
    ReciprocalLattice(MMatrix{3,3}(SDiagonal(oneunit(T), oneunit(T), oneunit(T))))
Base.oneunit(lattice::ReciprocalLattice) = oneunit(typeof(lattice))

Base.parent(lattice::ReciprocalLattice) = lattice.data

Base.size(::ReciprocalLattice) = (3, 3)

Base.getindex(lattice::ReciprocalLattice, i::Int) = getindex(parent(lattice), i)
Base.getindex(lattice::ReciprocalLattice, I...) = getindex(parent(lattice), I...)

Base.setindex!(lattice::ReciprocalLattice, v, i::Int) = setindex!(parent(lattice), v, i)
Base.setindex!(lattice::ReciprocalLattice, X, I...) = setindex!(parent(lattice), X, I...)

Base.IndexStyle(::Type{ReciprocalLattice{T}}) where {T} = IndexLinear()

# Customizing broadcasting
# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L397-L398
# and https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/structuredbroadcast.jl#L7-L14
struct ReciprocalLatticeStyle <: Broadcast.AbstractArrayStyle{2} end
ReciprocalLatticeStyle(::Val{2}) = ReciprocalLatticeStyle()
ReciprocalLatticeStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

Base.BroadcastStyle(::Type{<:ReciprocalLattice}) = ReciprocalLatticeStyle()

Base.similar(::Broadcast.Broadcasted{ReciprocalLatticeStyle}, ::Type{T}) where {T} =
    similar(ReciprocalLattice{T}, 3, 3)
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta2/base/abstractarray.jl#L839
function Base.similar(lattice::ReciprocalLattice, ::Type{T}, dims::Dims) where {T}
    if dims == size(lattice)
        return ReciprocalLattice(similar(Matrix{T}, dims))
    else
        throw(ArgumentError("invalid dims `$dims` for `Lattice`!"))
    end
end
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/base/abstractarray.jl#L874
function Base.similar(::Type{ReciprocalLattice{T}}, dims::Dims) where {T}
    if dims == (3, 3)
        return ReciprocalLattice(similar(Matrix{T}, dims))
    else
        throw(ArgumentError("invalid dims `$dims` for `ReciprocalLattice`!"))
    end
end

Base.:*(::ReciprocalLattice, ::ReciprocalLattice) =
    error("undefined operation `*` for `ReciprocalLattice`s!")

Base.:/(::ReciprocalLattice, ::ReciprocalLattice) =
    error("undefined operation `/` for `ReciprocalLattice`s!")
