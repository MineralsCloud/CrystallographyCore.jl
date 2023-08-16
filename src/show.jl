function Base.show(io::IO, ::MIME"text/plain", lattice::Lattice)
    summary(io, lattice)
    println(io)
    join(io, ' ' * join(row, "  ") * '\n' for row in eachrow(lattice.data))
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    summary(io, cell)
    println(io)
    println(io, " lattice:")
    for row in eachrow(cell.lattice)
        println(io, "   ", join(row, "  "))
    end
    num_atom = natoms(cell)
    println(io, " $num_atom atomic positions:")
    for position in cell.positions
        println(io, "   ", join(position, "  "))
    end
    println(io, " $num_atom atoms:")
    println(io, "   ", join(cell.atoms, "  "))
    return nothing
end
