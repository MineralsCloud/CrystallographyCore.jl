struct Inverted{T<:AbstractLattice}
    lattice::T
end

Base.inv(lattice::AbstractLattice) = Inverted(lattice)
Base.inv(inverted::Inverted) = inverted.lattice

(inverted::Inverted)(cartesian::AbstractVector) = parent(inverted.lattice) \ cartesian

(lattice::AbstractLattice)(reduced::AbstractVector) = parent(lattice) * reduced
