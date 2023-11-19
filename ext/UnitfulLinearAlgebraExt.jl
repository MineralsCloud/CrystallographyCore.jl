module UnitfulLinearAlgebraExt

using CrystallographyCore: Inverted, AbstractLattice
using Unitful: AbstractQuantity, NoUnits
using UnitfulLinearAlgebra: UnitfulMatrix

# See https://github.com/ggebbie/UnitfulLinearAlgebra.jl/issues/92
(inverted::Inverted{<:AbstractLattice{<:AbstractQuantity}})(
    cartesian::AbstractVector{<:AbstractQuantity}
) = collect(
    NoUnits.(
        vec(UnitfulMatrix(parent(inverted.lattice)) \ UnitfulMatrix(collect(cartesian)))
    ),
)

end
