using StaticArrays: StaticMatrix, MMatrix, Size

"""
    AbstractLattice{T} <: AbstractMatrix{T}

Represent the real lattices and the reciprocal lattices.
"""
abstract type AbstractLattice{T} <: StaticMatrix{3,3,T} end
mutable struct Lattice{T} <: AbstractLattice{T}
    data::MMatrix{3,3,T,9}
end
# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L195-L198
Lattice{T}(::UndefInitializer) where {T} = Lattice(MMatrix{3,3,T,9}(undef))
"""
    Lattice(data::AbstractMatrix)

Construct a `Lattice` from a matrix.

!!! note
    The basis vectors of the matrix are stored as columns.

# Examples
```jldoctest
julia> Lattice([
    1.2 4.5 7.8
    2.3 5.6 8.9
    3.4 6.7 9.1
])
3×3 Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1
```
"""
Lattice(data::AbstractMatrix) = Lattice(MMatrix{3,3}(data))
"""
    Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Construct a `Lattice` from three basis vectors.

# Examples
```jldoctest
julia> 𝐚, 𝐛, 𝐜 = [1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.10];

julia> Lattice(𝐚, 𝐛, 𝐜)
3×3 Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1
```
"""
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
"""
    Lattice(data)

Construct a `Lattice` from, e.g., a vector of three basis vectors.

# Examples
```jldoctest
julia> Lattice([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.10]])
3×3 Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1

julia> Lattice(((1.1, 2.2, 3.1), (4.4, 5.5, 6.5), (7.3, 8.8, 9.9)))
3×3 Lattice{Float64}
 1.1  4.4  7.3
 2.2  5.5  8.8
 3.1  6.5  9.9

julia> Lattice((1.1, 2.2, 3.1, 4.4, 5.5, 6.5, 7.3, 8.8, 9.9))
3×3 Lattice{Float64}
 1.1  4.4  7.3
 2.2  5.5  8.8
 3.1  6.5  9.9

julia> Lattice(i * 1.1 for i in 1:9)
3×3 Lattice{Float64}
 1.1  4.4  7.700000000000001
 2.2  5.5  8.8
 3.3000000000000003  6.6000000000000005  9.9

julia> using Unitful, UnitfulAtomic

julia> Lattice([
    [4u"nm", 0u"m", 0.0u"cm"],
    [0u"cm", 180.0u"bohr", 0u"m"],
    [0u"bohr", 0u"nm", (3//1) * u"angstrom"],
])
3×3 Lattice{Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}
 4.0e-9 m  0.0 m  0.0 m
 0.0 m  9.525189796254e-9 m  0.0 m
 0.0 m  0.0 m  3.0e-10 m
```
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

Base.parent(lattice::Lattice) = lattice.data

Base.getindex(op::Lattice, i::Int) = getindex(parent(op), i)

Base.setindex!(op::Lattice, v, i::Int) = setindex!(parent(op), v, i)

# Customizing broadcasting
# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L397-L398
# and https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/structuredbroadcast.jl#L7-L14
struct LatticeStyle <: Broadcast.AbstractArrayStyle{2} end
LatticeStyle(::Val{2}) = LatticeStyle()
LatticeStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

Base.BroadcastStyle(::Type{<:Lattice}) = LatticeStyle()

Base.similar(::Broadcast.Broadcasted{LatticeStyle}, ::Type{T}) where {T} =
    similar(Lattice{T})
# Override https://github.com/JuliaArrays/StaticArrays.jl/blob/v1.6.2/src/abstractarray.jl#L129
function Base.similar(op::Lattice, ::Type{T}, _size::Size) where {T}
    if _size == size(op)
        Lattice{T}(undef)
    else
        return similar(Array(op), T, _size)
    end
end
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta2/base/abstractarray.jl#L839
function Base.similar(op::Lattice, ::Type{T}, dims::Dims) where {T}
    if dims == size(op)
        Lattice{T}(undef)
    else
        return similar(Array(op), T, dims)
    end
end
function Base.similar(::Type{<:Lattice}, ::Type{T}, s::Size) where {T}
    if s == (3, 3)
        Lattice{T}(undef)
    else
        return Array{T}(undef, Tuple(s))
    end
end
function Base.similar(::Type{<:Lattice}, ::Type{T}, dim, dims...) where {T}
    if (dim, dims...) == (3, 3)
        Lattice{T}(undef)
    else
        return Array{T}(undef, dim, dims...)
    end
end