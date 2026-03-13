using StaticArrays: FieldVector, MVector, Size, SVector
using StructEquality: @struct_hash_equal_isequal

import StaticArrays: similar_type

export ReducedCoordinates, CrystalCoordinates, Cell, natoms, atomtypes, atomcounts

struct ReducedCoordinates{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end
const CrystalCoordinates = ReducedCoordinates

abstract type AbstractCell end
@struct_hash_equal_isequal mutable struct Cell{L,P,T} <: AbstractCell
    lattice::Lattice{L}
    positions::Vector{ReducedCoordinates{P}}
    atoms::Vector{T}
end
"""
    Cell(lattice, positions, atoms)

Create a new cell.

Argument `lattice` is a [`Lattice`](@ref) type.
Fractional atomic positions `positions` are given
by a vector of ``N`` vectors with floating point values, where ``N`` is the number of atoms.
Argument `atoms` is a list of ``N`` values, where the same kind of atoms
need to be the same type.
"""
function Cell(lattice, positions, atoms)
    if !(lattice isa Lattice)
        lattice = Lattice(lattice)
    end
    if length(positions) != length(atoms)
        throw(DimensionMismatch("the lengths of atomic positions and atoms are different!"))
    end
    P = reduce(promote_type, eltype.(positions))
    positions = map(ReducedCoordinates{P} ∘ collect, positions)
    L, T = eltype(lattice), eltype(atoms)
    return Cell{L,P,T}(lattice, positions, atoms)
end

"""
    natoms(cell::Cell)

Return the number of atoms contained in `cell`.
"""
natoms(cell::Cell) = length(cell.atoms)

"""
    atomtypes(cell::Cell)

Return the distinct atom values present in `cell`.
"""
atomtypes(cell::Cell) = unique(cell.atoms)

"""
    atomcounts(cell::Cell)

Return a dictionary mapping each distinct atom value in `cell.atoms`
to the number of occurrences in the cell.

Examples
```julia
julia> cell = Cell(rand(3, 3), [rand(3) for _ in 1:4], [:C, :O, :C, :H])

julia> atomcounts(cell)
Dict{Symbol, Int64} with 3 entries:
  :H => 1
  :O => 1
  :C => 2
```
"""
function atomcounts(cell::Cell)
    counts = Dict{eltype(cell.atoms),Int64}()
    for atom in cell.atoms
        counts[atom] = get(counts, atom, 0) + 1
    end
    return counts
end

"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = cell.lattice

# See https://juliaarrays.github.io/StaticArrays.jl/dev/pages/api/#StaticArraysCore.FieldVector
similar_type(::Type{<:ReducedCoordinates}, ::Type{T}, s::Size{(3,)}) where {T} =
    ReducedCoordinates{T}

struct DynamicCell{C<:Cell} <: AbstractCell
    cell::C
    lattice_velocities::SVector{3,MVector{3,Float64}}
    ion_velocities::Vector{MVector{3,Float64}}
end
