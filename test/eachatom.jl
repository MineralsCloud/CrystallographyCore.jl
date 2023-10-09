@testset "Test `eachatom`" begin
    lattice = Lattice(
        [
            4.59983732 0.0 0.0
            0.0 4.59983732 0.0
            0.0 0.0 2.95921356
        ] .* u"angstrom"
    )
    positions = [
        [0, 0, 0.5],
        [0.5, 0.5, 0.0],
        [0.19567869, 0.80432131, 0],
        [0.80432131, 0.19567869, 0],
        [0.30432131, 0.30432131, 0.5],
        [0.69567869, 0.69567869, 0.5],
    ]
    atoms = [22, 22, 8, 8, 8, 8]
    cell = Cell(lattice, positions, atoms)
    for ((atom, position), atom′, position′) in zip(eachatom(cell), atoms, positions)
        @test atom == atom′
        @test position == position′
    end
end
