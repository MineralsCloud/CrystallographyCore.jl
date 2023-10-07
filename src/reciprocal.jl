using LinearAlgebra: I, det, cross

import Base: *, /

export ReciprocalLattice, reciprocal, isreciprocal

"""
    ReciprocalLattice(data::AbstractMatrix)

Construct a `ReciprocalLattice` from a matrix.

!!! note
    The basis vectors of the matrix are stored as columns.

!!! warning
    Avoid using this constructor directly. Use `reciprocal` instead.
"""
@struct_hash_equal_isequal_isapprox struct ReciprocalLattice{T} <: AbstractLattice{T}
    data::MMatrix{3,3,T,9}
end
ReciprocalLattice(data::AbstractMatrix) = ReciprocalLattice(MMatrix{3,3}(data))

"""
    basisvectors(lattice::ReciprocalLattice)

Get the three basis vectors from a reciprocal lattice.
"""
basisvectors(lattice::ReciprocalLattice) =
    lattice[begin:(begin + 2)], lattice[(begin + 3):(begin + 5)], lattice[(begin + 6):end]

"""
    reciprocal(lattice::Lattice)
    reciprocal(lattice::ReciprocalLattice)

Get the reciprocal of a `Lattice` or a `ReciprocalLattice`.
"""
function reciprocal(lattice::Lattice)
    Ω = det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚, 𝐛, 𝐜 = basisvectors(lattice)
    return inv(Ω) * ReciprocalLattice(hcat(cross(𝐛, 𝐜), cross(𝐜, 𝐚), cross(𝐚, 𝐛)))
end
function reciprocal(lattice::ReciprocalLattice)
    Ω⁻¹ = det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚⁻¹, 𝐛⁻¹, 𝐜⁻¹ = basisvectors(lattice)
    return inv(Ω⁻¹) * Lattice(hcat(cross(𝐛⁻¹, 𝐜⁻¹), cross(𝐜⁻¹, 𝐚⁻¹), cross(𝐚⁻¹, 𝐛⁻¹)))
end

isreciprocal(a::ReciprocalLattice, b::Lattice) = parent(a)' * parent(b) ≈ I
isreciprocal(a::Lattice, b::ReciprocalLattice) = isreciprocal(b, a)

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

Base.IndexStyle(::Type{<:ReciprocalLattice}) = IndexLinear()

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

for op in (:*, :/)
    for S in (:Lattice, :ReciprocalLattice)
        for T in (:Lattice, :ReciprocalLattice)
            @eval $op(::$S, ::$T) = error("undefined operation!")
        end
    end
end
