using JULINT
using Test

@testset "JULINT.jl" begin
    @test JULINT.JULINT_greet() == "Welcome to JULINT!"
    @test JULINT.JULINT_greet() != "Hello world!"
end
