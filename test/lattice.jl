using CrystallographyCore: Lattice, basisvectors
using StaticArrays: SMatrix
using Unitful, UnitfulAtomic

@testset "Constructing `Lattice`s" begin
    @testset "with matrix" begin
        # General 3x3 matrix
        mat = [
            1.2 4.5 7.8
            2.3 5.6 8.9
            3.4 6.7 9.1
        ]
        @test Lattice(mat) == Lattice(SMatrix{3,3}(mat))
        @test basisvectors(Lattice(mat)) ==
            ([1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.1])
        @test all(basisvectors(Lattice(mat)) .== basisvectors(Lattice(mat)))
        # Rectangular matrix
        @test_throws DimensionMismatch Lattice([
            1 2
            3 4
            5 6
        ])
        # Complex numbers
        mat = [
            1.2+2im 4.5 7.8
            2.3 5.6 8.9
            3.4 6.7 9.1
        ]
        @test Lattice(mat) == Lattice(SMatrix{3,3}(mat))
        # Ragged array
        @test_throws DimensionMismatch Lattice([[1, 2], [3, 4, 5], [6]])
    end
    # Column vectors
    a = [1.1, 2.2, 3.1]
    b = [4.4, 5.5, 6.5]
    c = [7.3, 8.8, 9.9]
    @test Lattice(a, b, c) == Lattice([
        1.1 4.4 7.3
        2.2 5.5 8.8
        3.1 6.5 9.9
    ])
    @test basisvectors(Lattice(a, b, c)) == (a, b, c)
    @testset "with general iterables" begin
        # Tuple with 9 values
        vals = (1.1, 2.2, 3.1, 4.4, 5.5, 6.5, 7.3, 8.8, 9.9)
        @test Lattice(vals) == Lattice([
            1.1 4.4 7.3
            2.2 5.5 8.8
            3.1 6.5 9.9
        ])
        # Tuple of tuples
        vals = ((1.1, 2.2, 3.1), (4.4, 5.5, 6.5), (7.3, 8.8, 9.9))
        @test Lattice(vals) == Lattice([
            1.1 4.4 7.3
            2.2 5.5 8.8
            3.1 6.5 9.9
        ])
        # Vector of vectors
        vals = [[1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.10]]
        @test Lattice(vals) == Lattice([
            1.2 4.5 7.8
            2.3 5.6 8.9
            3.4 6.7 9.1
        ])
        @test basisvectors(Lattice(vals)) ==
            ([1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.1])
        # Vector of tuples
        @test basisvectors(Lattice([(4.01, 0, 0), (0, 4, 0), (0, 0, 3.99)])) ==
            ([4.01, 0, 0], [0, 4, 0], [0, 0, 3.99])
        # Generator of 9 values
        vals = (i * 1.1 for i in 1:9)
        @test Lattice(vals) == Lattice([
            1 4 7
            2 5 8
            3 6 9
        ] * 1.1)
    end
    @testset "with units" begin
        lattice = Lattice([
            [4u"nm", 0u"m", 0.0u"cm"],
            [0u"cm", 180.0u"bohr", 0u"m"],
            [0u"bohr", 0u"nm", (3//1) * u"angstrom"],
        ])
        @test lattice.data == [
            4u"nm" 0u"m" 0.0u"cm"
            0u"cm" 180.0u"bohr" 0u"m"
            0u"bohr" 0u"nm" (3//1)*u"angstrom"
        ]
        @test basisvectors(lattice) == (
            [4u"nm", 0u"m", 0.0u"cm"],
            [0u"cm", 180.0u"bohr", 0u"m"],
            [0u"bohr", 0u"nm", (3//1) * u"angstrom"],
        )
        @test inv(inv(lattice)) == lattice
    end
end

@testset "Test broadcasting for lattices" begin
    lattice = Lattice([1, 0, 0], [0, 1, 0], [0, 0, 1])
    @test lattice .* 4 == 4 .* lattice == Lattice([4, 0, 0], [0, 4, 0], [0, 0, 4])
    @test lattice .* 4.0 == 4.0 .* lattice == Lattice([4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0])
    @test lattice .* 4//1 == 4//1 .* lattice == Lattice([4, 0, 0], [0, 4, 0], [0, 0, 4])
    @test lattice ./ 4 == Lattice([1//4, 0, 0], [0, 1//4, 0], [0, 0, 1//4])
    @test lattice .* u"nm" ==
        u"nm" .* lattice ==
        Lattice(
            [1u"nm", 0u"nm", 0u"nm"], [0u"nm", 1u"nm", 0u"nm"], [0u"nm", 0u"nm", 1u"nm"]
        )
    @test lattice .* 1u"nm" ==
        1u"nm" .* lattice ==
        Lattice(
            [1u"nm", 0u"nm", 0u"nm"], [0u"nm", 1u"nm", 0u"nm"], [0u"nm", 0u"nm", 1u"nm"]
        )
    @test_throws ArgumentError 4 / lattice
    @test_throws ArgumentError 4.0 ./ lattice
end

@testset "Test broadcasting for reciprocal lattices" begin
    a, b, c = 4, 3, 5
    lattice = Lattice([a, -b, 0] / 2, [a, b, 0] / 2, [0, 0, c])
    @test reciprocal(lattice .* 4) == reciprocal(lattice) ./ 4
    @test_throws ArgumentError 4 / reciprocal(lattice)
    @test_throws ArgumentError 4.0 ./ reciprocal(lattice)
end
