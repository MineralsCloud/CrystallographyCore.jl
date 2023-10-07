using StaticArrays: FieldVector
using StructEquality: @struct_hash_equal_isequal

export ReducedCoordinates, CrystalCoordinates, Cell, natoms, atomtypes

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

natoms(cell::Cell) = length(cell.atoms)

atomtypes(cell::Cell) = unique(cell.atoms)

"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = cell.lattice
