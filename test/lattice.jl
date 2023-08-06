using Unitful, UnitfulAtomic

@testset "Test creating `Lattice`s with units" begin
    lattice = Lattice(
        [
            4u"nm" 0u"m" 0.0u"cm"
            0u"cm" 180.0u"bohr" 0u"m"
            0u"bohr" 0u"nm" (3//1)*u"angstrom"
        ]
    )
    @test lattice.data == [
        4u"nm" 0u"m" 0.0u"cm"
        0u"cm" 180.0u"bohr" 0u"m"
        0u"bohr" 0u"nm" (3//1)*u"angstrom"
    ]
end
