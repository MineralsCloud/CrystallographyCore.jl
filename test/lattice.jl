using CrystallographyCore: Lattice, basisvectors
using StaticArrays: SMatrix
using Unitful, UnitfulAtomic

@testset "Test `isbits`" begin
    @test isbitstype(Lattice{Int})
    @test isbits(Lattice([1, 0, 0], [0, 1, 0], [0, 0, 1]))
end

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
    @testset "Test creating `Lattice`s from 6 lattice constants" begin
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L96
        @test basisvectors(Lattice(2, 1, 5, 90, 90, 90; axis=:c)) ==
            ([2, 0, 0], [0, 1, 0], [0, 0, 5])  # Orthorombic
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L104
        @test basisvectors(Lattice(1, 2, 3, 90, 120, 90; axis=:c)) ==
            ([0.8660254037844387, 0, -0.5], [0, 2, 0], [0, 0, 3])  # Monoclinic
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L117
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(1, 2, 3, 75, 40, 81; axis=:c)) .≈
            ([0.641327, -0.04330811, 0.76604444], [0, 1.93185165, 0.51763809], [0, 0, 3]),
        )  # Triclinic
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L131
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(3, 4, 20, 45, 90, 126; axis=:c)) .≈ (
                [1.66767891, -2.49376163, 1.83697020e-16],
                [0, 2.82842712, 2.82842712],
                [0, 0, 20],
            ),
        )  # Triclinic
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L158
        @test all(
            basisvectors(Lattice(2, 2, 3, 90, 90, 90; axis=:c)) .≈
            ([2, 0, 0], [0, 2, 0], [0, 0, 3]),
        )  # Tetragonal
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L165
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(1, 1, 1, 87, 87, 87; axis=:c)) .≈
            ([0.99739377, 0.04966497, 0.05233596], [0, 0.99862953, 0.05233596], [0, 0, 1]),
        )  # Rhombohedral
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L173
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(1, 2, 3, 90, 115, 90; axis=:c)) .≈
            ([0.906307787, 8.71102450e-17, -0.422618262], [0, 2, 1.2246468e-16], [0, 0, 3]),
        )  # Monoclinic
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L177
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(2, 3, 1, 115, 90, 90; axis=:c)) .≈
            ([2, 1.92231042e-16, 1.22464680e-16], [0, 2.71892336, -1.26785479], [0, 0, 1]),
        )  # Monoclinic
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L181
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(3, 1, 2, 90, 90, 115; axis=:c)) .≈
            ([2.71892336, -1.26785479, 1.83697020e-16], [0, 1, 6.123234e-17], [0, 0, 2]),
        )  # Monoclinic
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L189
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(2, 2, 3, 90, 90, 120; axis=:c)) .≈
            ([1.73205081, -1, 1.22464680e-16], [0, 2, 1.2246468e-16], [0, 0, 3]),
        )  # Hexagonal
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L193
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(3, 2, 2, 120, 90, 90; axis=:c)) .≈
            ([3, 3.18172572e-16, 1.83697020e-16], [0, 1.73205081, -1], [0, 0, 2]),
        )  # Hexagonal
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L197
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(2, 3, 2, 90, 120, 90; axis=:c)) .≈
            ([1.73205081, 1.83697020e-16, -1], [0, 3, 1.8369702e-16], [0, 0, 2]),
        )  # Hexagonal
        # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L201
        # Compared with Python's results
        @test all(
            basisvectors(Lattice(2, 2, 2, 90, 120, 90; axis=:c)) .≈
            ([1.73205081, 1.83697020e-16, -1], [0, 2, 1.2246468e-16], [0, 0, 2]),
        )  # Hexagonal
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
