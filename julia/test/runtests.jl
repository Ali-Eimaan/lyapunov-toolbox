# Entry point for `Pkg.test()`.
#
# Julia requires the test entry point to live at <package>/test/runtests.jl, but
# the repository specification puts the test files in a top-level tests/
# directory shared with the Python and MATLAB sides. This file bridges the two:
# it is the only thing in julia/test/, and it contains no assertions of its own.
# Add new test files to tests/ and include them here.

using Test
using LyapunovToolbox

const TESTS = joinpath(@__DIR__, "..", "..", "tests")

@testset "LyapunovToolbox.jl" begin
    include(joinpath(TESTS, "test_polynomial_utils.jl"))
    include(joinpath(TESTS, "test_sos_lyapunov.jl"))
    include(joinpath(TESTS, "test_roa_estimation.jl"))
end
