export eachatom

struct EachAtom{A,B}
    atoms::Vector{A}
    positions::Vector{B}
end

"""
    eachatom(cell::Cell)

Create a generator that iterates over the atoms in a `Cell`.
"""
eachatom(cell::Cell) = EachAtom(cell.atoms, cell.positions)

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
function Base.iterate(iter::EachAtom, state=1)
    if state > length(iter)
        return nothing
    else
        return (iter.atoms[state], iter.positions[state]), state + 1
    end
end

Base.eltype(::EachAtom{A,B}) where {A,B} = Tuple{A,B}

Base.length(iter::EachAtom) = length(iter.atoms)

Base.IteratorSize(::Type{<:EachAtom}) = Base.HasLength()
