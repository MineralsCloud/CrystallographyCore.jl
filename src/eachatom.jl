export eachatom

struct EachAtom{N,A,B}
    atoms::NTuple{N,A}
    positions::NTuple{N,B}
end

"""
    eachatom(cell::Cell)

Create a generator that iterates over the atoms in a `Cell`.
"""
eachatom(cell::Cell) = EachAtom(Tuple(cell.atoms), Tuple(cell.positions))

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
function Base.iterate(iter::EachAtom, state=1)
    if state > length(iter)
        return nothing
    else
        return (iter.atoms[state], iter.positions[state]), state + 1
    end
end

Base.eltype(::Type{EachAtom{N,A,B}}) where {N,A,B} = Tuple{A,B}

Base.length(::EachAtom{N}) where {N} = N

Base.IteratorSize(::Type{<:EachAtom}) = Base.HasLength()
