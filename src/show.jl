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
    for row in eachrow(cell.lattice.data)
        println(io, "  ", join(row, "  "))
    end
    N = natoms(cell)
    println(io, " $N atomic positions:")
    for pos in cell.positions
        println(io, "  ", pos)
    end
    println(io, " $N atoms:")
    println(io, "  ", cell.atoms)
    return nothing
end
