using Test
using FusionRings

include("test_utils.jl")
include("test_fixtures.jl")

@testset "FusionRings.jl" begin
    include("creation_tests.jl")
    include("properties_tests.jl")
    include("operations_tests.jl")
    include("subring_tests.jl")
    include("automorphism_tests.jl")
    include("numeric_tests.jl")
    include("import_export_tests.jl")
    include("formatting_tests.jl")
end