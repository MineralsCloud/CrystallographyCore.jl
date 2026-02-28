using CrystallographyCore
using Test

@testset "CrystallographyCore.jl" begin
    # Write your tests here.
    include("lattice.jl")
    include("eachatom.jl")
    include("transform.jl")
end
