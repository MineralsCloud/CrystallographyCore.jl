struct Inverted{T<:AbstractLattice}
    lattice::T
end

Base.inv(lattice::AbstractLattice) = Inverted(lattice)
Base.inv(inverted::Inverted) = inverted.lattice

(inverted::Inverted)(reduced::AbstractVector) = parent(inverted.lattice) \ reduced

(lattice::AbstractLattice)(reduced::AbstractVector) = parent(lattice) * reduced
