using StaticArrays: SVector, FieldVector
using StructEquality: @struct_hash_equal_isequal

export Cell, natoms, atomtypes

struct CrystalCoordinates{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

abstract type AbstractCell end
@struct_hash_equal_isequal struct Cell{N,L,P,T} <: AbstractCell
    lattice::Lattice{L}
    positions::SVector{N,CrystalCoordinates{P}}
    atoms::SVector{N,T}
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
    typeassert(positions, AbstractVector{<:AbstractVector})
    if length(positions) != length(atoms)
        throw(DimensionMismatch("the lengths of atomic positions and atoms are different!"))
    end
    N = length(positions)
    P = reduce(promote_type, eltype.(positions))
    positions = map(CrystalCoordinates{P}, positions)
    L, T = eltype(lattice), eltype(atoms)
    return Cell{N,L,P,T}(lattice, positions, atoms)
end

natoms(::Cell{N}) where {N} = N

atomtypes(cell::Cell) = unique(cell.atoms)

"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = cell.lattice
