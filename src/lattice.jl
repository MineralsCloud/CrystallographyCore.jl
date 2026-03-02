using ConstructionBase: constructorof
using StaticArrays: SMatrix, SDiagonal
using StructEquality: @struct_hash_equal_isequal_isapprox

export Lattice, basisvectors

"""
    AbstractLattice{T}

Represent the real lattices and the reciprocal lattices.
"""
abstract type AbstractLattice{T} end
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
Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1
```
"""
@struct_hash_equal_isequal_isapprox struct Lattice{T} <: AbstractLattice{T}
    data::SMatrix{3,3,T,9}
end
Lattice(data::AbstractMatrix) = Lattice(SMatrix{3,3}(data))
"""
    Lattice(a, b, c, α, β, γ; axis = :a)

Construct a `Lattice` from the six cell parameters.

The default convention we used here is that edge vector 𝐚 in the positive x-axis direction,
edge vector 𝐛 in the x-y plane with a positive y-axis component,
and edge vector 𝐜 with a positive z-axis component in the Cartesian system.
See [Wikipedia](https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography).
You can also choose `axis = :c`.
"""
function Lattice(a, b, c, α, β, γ; axis=:a)
    Ω =
        a *
        b *
        c *
        sqrt(sind(α)^2 - cosd(β)^2 - cosd(γ)^2 + 2 * cosd(α) * cosd(β) * cosd(γ))
    if axis == :a  # See https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography
        sinγ, cosγ, cosα, cosβ, 𝟎 = sind(γ), cosd(γ), cosd(α), cosd(β), zero(a)
        return Lattice(
            [a, 𝟎, 𝟎],
            [b * cosγ, b * sinγ, zero(b)],
            [c * cosβ, c * (cosα - cosβ * cosγ) / sinγ, Ω / (a * b * sinγ)],
        )
    elseif axis == :c  # See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L356-L391
        sinα, cosα, sinβ, cosβ, 𝟎 = sind(α), cosd(α), sind(β), cosd(β), zero(c)
        x = Ω / (b * c * sinα)
        cos′ = (cosα * cosβ - cosd(γ)) / (sinα * sinβ)
        sin′ = sqrt(1 - cos′^2)
        return Lattice(
            [x, -x * cos′ / sin′, a * cosβ], [zero(b), b * sinα, b * cosα], [𝟎, 𝟎, c]
        )
    else
        error("aligning `$axis` axis is not supported!")
    end
end
"""
    Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Construct a `Lattice` from three basis vectors.

# Examples
```jldoctest
julia> 𝐚, 𝐛, 𝐜 = [1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.10];

julia> Lattice(𝐚, 𝐛, 𝐜)
Lattice{Float64}
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
Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1

julia> Lattice(((1.1, 2.2, 3.1), (4.4, 5.5, 6.5), (7.3, 8.8, 9.9)))
Lattice{Float64}
 1.1  4.4  7.3
 2.2  5.5  8.8
 3.1  6.5  9.9

julia> Lattice((1.1, 2.2, 3.1, 4.4, 5.5, 6.5, 7.3, 8.8, 9.9))
Lattice{Float64}
 1.1  4.4  7.3
 2.2  5.5  8.8
 3.1  6.5  9.9

julia> Lattice(i * 1.1 for i in 1:9)
Lattice{Float64}
 1.1  4.4  7.700000000000001
 2.2  5.5  8.8
 3.3000000000000003  6.6000000000000005  9.9

julia> using Unitful, UnitfulAtomic

julia> Lattice(
           [4u"nm", 0u"m", 0.0u"cm"],
           [0u"cm", 180.0u"bohr", 0u"m"],
           [0u"bohr", 0u"nm", (3//1) * u"angstrom"],
       )
Lattice{Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}
 4.0e-9 m  0.0 m  0.0 m
 0.0 m  9.525189796254e-9 m  0.0 m
 0.0 m  0.0 m  3.0e-10 m
```
"""
Lattice(data::NTuple{9}) = Lattice(SMatrix{3,3}(data))
Lattice(data::NTuple{3,NTuple{3}}) = Lattice(mapreduce(collect, hcat, data))
Lattice(data::AbstractVector{<:AbstractVector}) = Lattice(hcat(data[1], data[2], data[3]))  # This is faster
Lattice(iter::Base.Generator) = Lattice(SMatrix{3,3}(iter))
function Lattice(data::Union{AbstractVector,Tuple})
    if length(data) == 9
        return Lattice(SMatrix{3,3}(data))
    elseif length(data) == 3
        return Lattice(mapreduce(collect, hcat, data))
    else
        throw(DimensionMismatch("The length of the tuple must be 3 or 9."))
    end
end

"""
    basisvectors(lattice::Lattice)

Get the three basis vectors from a lattice.
"""
basisvectors(lattice::Lattice) =
    lattice[begin:(begin + 2)], lattice[(begin + 3):(begin + 5)], lattice[(begin + 6):end]

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L130-L131
Base.one(::Type{Lattice{T}}) where {T} = Lattice(SDiagonal(one(T), one(T), one(T)))
Base.one(lattice::Lattice) = one(typeof(lattice))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L132-L133
Base.oneunit(::Type{Lattice{T}}) where {T} =
    Lattice(SDiagonal(oneunit(T), oneunit(T), oneunit(T)))
Base.oneunit(lattice::Lattice) = oneunit(typeof(lattice))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L134-L135
Base.zero(::Type{Lattice{T}}) where {T} = Lattice(zeros(T, 3, 3))
Base.zero(lattice::Lattice) = zero(typeof(lattice))

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
Base.iterate(iter::Lattice) = iterate(parent(iter))
Base.iterate(iter::Lattice, state) = iterate(parent(iter), state)

Base.IteratorSize(::Type{<:Lattice}) = Base.HasShape{2}()

Base.eltype(::Type{Lattice{T}}) where {T} = T

Base.length(::Lattice) = 9

Base.size(::Lattice) = (3, 3)
# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L74
Base.size(lattice::Lattice, dim) = size(parent(lattice), dim)  # Here, `parent(A)` is necessary to avoid `StackOverflowError`.

Base.parent(lattice::Lattice) = lattice.data

Base.getindex(lattice::Lattice, i...) = getindex(parent(lattice), i...)

Base.firstindex(::Lattice) = 1

Base.lastindex(::Lattice) = 9

Base.IndexStyle(::Type{<:Lattice}) = Base.IndexCartesian()
Base.IndexStyle(::Lattice) = Base.IndexCartesian()

# You need this to let the broadcasting work.
Base.:*(lattice::Lattice, x::Number) = Lattice(parent(lattice) * x)
Base.:*(x::Number, lattice::Lattice) = lattice * x
"""
    *(R::AbstractMatrix, lattice::Lattice)

Left-multiply a lattice by a matrix and return a `Lattice`.

If the lattice matrix is ``\\mathrm{A} = [\\mathbf{a}\\ \\mathbf{b}\\ \\mathbf{c}]``
(basis vectors as columns), this computes ``\\mathrm{R}\\mathrm{A}``.
This corresponds to an active (reverse/alibi) transformation, such as
a rigid Cartesian rotation applied to the crystal basis vectors.

Only ``3×3`` matrices are supported.

See also "[Left and right matrix actions on a lattice](@ref lattice_matrix_actions)".
"""
function Base.:*(R::AbstractMatrix, lattice::Lattice)
    size(R) == (3, 3) || throw(
        DimensionMismatch(
            "matrix size $(size(R)) is invalid: only 3×3 matrices can multiply a lattice",
        ),
    )
    return Lattice(R * parent(lattice))
end
"""
    *(lattice::Lattice, P::AbstractMatrix)

Right-multiply a lattice by a matrix and return a `Lattice`.

If the lattice matrix is ``\\mathrm{A} = [\\mathbf{a}\\ \\mathbf{b}\\ \\mathbf{c}]``
(basis vectors as columns), this computes ``\\mathrm{A}\\mathrm{P}``.
This corresponds to a passive (forward/alias) change of basis.

Only ``3×3`` matrices are supported.

See also "[Left and right matrix actions on a lattice](@ref lattice_matrix_actions)".
"""
function Base.:*(lattice::Lattice, P::AbstractMatrix)
    size(P) == (3, 3) || throw(
        DimensionMismatch(
            "matrix size $(size(P)) is invalid: only 3×3 matrices can multiply a lattice",
        ),
    )
    return Lattice(parent(lattice) * P)
end

# You need this to let the broadcasting work.
Base.:/(lattice::Lattice, x::Number) = Lattice(parent(lattice) / x)
Base.:/(::Number, ::Lattice) =
    throw(ArgumentError("you cannot divide a number by a lattice!"))

Base.:+(lattice::Lattice) = lattice
# You need this to let the broadcasting work.
Base.:+(lattice::Lattice, x::Number) = Lattice(parent(lattice) .+ x)
Base.:+(x::Number, lattice::Lattice) = lattice + x

Base.:-(lattice::Lattice) = Lattice(-parent(lattice))
# You need this to let the broadcasting work.
Base.:-(lattice::Lattice, x::Number) = Lattice(parent(lattice) .- x)
Base.:-(x::Number, lattice::Lattice) = -lattice + x

Base.convert(::Type{Lattice{T}}, lattice::Lattice{T}) where {T} = lattice
Base.convert(::Type{Lattice{T}}, lattice::Lattice{S}) where {S,T} =
    Lattice(convert(SMatrix{3,3,T,9}, parent(lattice)))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta3/base/refpointer.jl#L95-L96
Base.ndims(::Type{<:Lattice}) = 2
Base.ndims(::Lattice) = 2

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L741
Base.broadcastable(lattice::Lattice) = lattice

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L49
Base.BroadcastStyle(::Type{<:Lattice}) = Broadcast.Style{Lattice}()

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L135
Base.BroadcastStyle(::Broadcast.AbstractArrayStyle{0}, b::Broadcast.Style{Lattice}) = b

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L1114-L1119
Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{Lattice}}) = Lattice(x for x in bc)  # For uniary and binary functions

Base.broadcasted(::typeof(/), ::Number, ::Lattice) =
    throw(ArgumentError("you cannot divide a number by a lattice!"))

Base.permutedims(lattice::Lattice) = Lattice(permutedims(parent(lattice)))

function Base.transpose(lattice::Lattice)
    throw(
        DomainError(
            "transpose($lattice)",
            "`transpose` is intended for linear algebra usage, use `permutedims` instead.",
        ),
    )
end

function Base.convert(::Type{T}, lattice::Lattice) where {T<:AbstractMatrix}
    eltype_T = eltype(T)
    eltype_lattice = eltype(lattice)
    T_promoted = promote_type(eltype_T, eltype_lattice)
    T_intersect = typeintersect(T_promoted, eltype_lattice)
    if T_intersect !== Base.Bottom
        T_promoted = T_intersect
    end
    if !(T_promoted <: eltype_T)
        throw(TypeError(:convert, "promoted type", eltype_T, T_promoted))
    end
    C = constructorof(T)
    return C{T_promoted}(parent(lattice))
end
