module SpglibExt

using CrystallographyCore: Cell

import Spglib: SpglibCell, unwrap_convert

function unwrap_convert(cell::Cell)
    lattice, positions, atoms = cell.lattice, cell.positions, cell.atoms
    clattice = Base.cconvert(Matrix{Cdouble}, permutedims(lattice))
    cpositions = Base.cconvert(Matrix{Cdouble}, reduce(hcat, positions))
    atomtypes = unique(atoms)
    catoms = collect(Cint, findfirst(==(atom), atomtypes) for atom in atoms)  # Mapping between unique atom types and atom indices
    return clattice, cpositions, catoms
end

SpglibCell(cell::Cell) = SpglibCell(cell.lattice, cell.positions, cell.atoms)

end
