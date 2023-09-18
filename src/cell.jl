using StaticArrays: MVector
using StructEquality: @struct_hash_equal_isequal

export Cell, natoms, atomtypes

abstract type AbstractCell end
"""
    Cell(lattice, positions, atoms)

Create a new cell.

Argument `lattice` is a [`Lattice`](@ref) type.
Fractional atomic positions `positions` are given
by a vector of ``N`` vectors with floating point values, where ``N`` is the number of atoms.
Argument `atoms` is a list of ``N`` values, where the same kind of atoms
need to be the same type.
"""
@struct_hash_equal_isequal struct Cell{L,P,T} <: AbstractCell
    lattice::Lattice{L}
    positions::Vector{MVector{3,P}}
    atoms::Vector{T}
end
function Cell(lattice, positions, atoms)
    if !(lattice isa Lattice)
        lattice = Lattice(lattice)
    end
    typeassert(positions, AbstractVector{<:AbstractVector})
    P = eltype(Base.promote_typeof(positions...))
    positions = collect(map(MVector{3,P}, positions))
    for position in positions
        if any(abs.(position) .> oneunit(eltype(position)))
            throw(ArgumentError("atomic positions must be fractional!"))
        end
    end
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
