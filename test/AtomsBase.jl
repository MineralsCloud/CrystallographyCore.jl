using AtomsBase
using Unitful: ustrip, @u_str
using UnitfulLinearAlgebra

# See https://juliamolsim.github.io/AtomsBase.jl/stable/tutorial/#System-interface-and-conventions v0.3.5
@testset "Test example from AtomsBase.jl documentation" begin
    box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"  # Note the unit!
    bc = [Periodic(), Periodic(), Periodic()]
    hydrogen = FlexibleSystem(
        [Atom(:H, [0, 0, 1.0]u"bohr"), Atom(:H, [0, 0, 3.0]u"bohr")], box, bc
    )
    cell = Cell(
        box, [[0, 0, ustrip(u"Å", 0.1u"bohr")], [0, 0, ustrip(u"Å", 0.3u"bohr")]], [:H, :H]
    )
    @test Cell(hydrogen) == cell
    @test all(
        position(FlexibleSystem(cell)) .≈
        [uconvert.(u"Å", [0, 0, 1.0]u"bohr"), uconvert.(u"Å", [0, 0, 3.0]u"bohr")],
    )
end
