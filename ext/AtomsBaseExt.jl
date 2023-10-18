module AtomsBaseExt

using AtomsBase:
    AbstractSystem,
    Atom,
    bounding_box,
    periodicity,
    isinfinite,
    element_symbol,
    periodic_system
using CrystallographyCore: basisvectors, eachatom

import AtomsBase: FlexibleSystem
import CrystallographyCore: Lattice, Cell

Lattice(system::AbstractSystem) = Lattice(bounding_box(system))

function Cell(system::AbstractSystem)
    if !all(periodicity(system)) || isinfinite(system)
        error("the system is not periodic!")
    end
    lattice = Lattice(system)
    atoms = element_symbol(system)
    positions = position(system)
    reduced = inv(lattice).(positions)
    return Cell(lattice, reduced, atoms)
end

function FlexibleSystem(cell::Cell)
    lattice = Lattice(cell)
    box = collect(basisvectors(lattice))
    atomicpositions = map(eachatom(cell)) do (atom, position)
        Atom(atom, lattice(position))
    end
    return periodic_system(atomicpositions, box; fractional=true)
end

end
