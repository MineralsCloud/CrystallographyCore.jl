module UnitfulExt

using CrystallographyCore: Lattice, ReciprocalLattice
using Unitful: Quantity, DimensionError, Units, dimension

import Unitful: uconvert, ustrip

# See https://github.com/PainterQubits/Unitful.jl/blob/v1.17.0/src/conversion.jl#L102-L111
Base.convert(
    ::Type{Lattice{Quantity{T,D,U}}}, lattice::Lattice{Quantity{T,D,U}}
) where {T,D,U} = lattice
Base.convert(
    ::Type{ReciprocalLattice{Quantity{T,D,U}}}, lattice::ReciprocalLattice{Quantity{T,D,U}}
) where {T,D,U} = lattice
function Base.convert(::Type{Lattice{Quantity{T,D,U}}}, lattice::Lattice) where {T,D,U}
    if dimension(eltype(lattice)) == D
        return Lattice(map(Base.Fix1(convert, Quantity{T,D,U}), parent(lattice)))
    else
        throw(DimensionError(U(), first(lattice)))
    end
end
function Base.convert(
    ::Type{ReciprocalLattice{Quantity{T,D,U}}}, lattice::ReciprocalLattice
) where {T,D,U}
    if dimension(eltype(lattice)) == D
        return ReciprocalLattice(map(Base.Fix1(convert, Quantity{T,D,U}), parent(lattice)))
    else
        throw(DimensionError(U(), first(lattice)))
    end
end

# See https://github.com/PainterQubits/Unitful.jl/blob/v1.17.0/src/conversion.jl#L66-L74
uconvert(u::Units, lattice::Lattice) = Lattice(map(Base.Fix1(uconvert, u), parent(lattice)))
uconvert(u::Units, lattice::ReciprocalLattice) =
    ReciprocalLattice(map(Base.Fix1(uconvert, u), parent(lattice)))

ustrip(lattice::Lattice) = Lattice(map(ustrip, parent(lattice)))
ustrip(lattice::ReciprocalLattice) = ReciprocalLattice(map(ustrip, parent(lattice)))

end
