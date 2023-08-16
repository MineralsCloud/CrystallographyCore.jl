using StaticArrays: MMatrix

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
