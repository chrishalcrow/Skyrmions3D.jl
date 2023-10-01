using Skyrmions3D
using Test

@testset "Differentiation tests" begin include("do_differentiation_tests.jl") end
@testset "Initialise tests" begin include("do_initialise_tests.jl") end
@testset "Initialise tests" begin include("do_initialise_skyrmions_tests.jl") end
@testset "Properties tests" begin include("do_properties_tests.jl") end
@testset "Transform tests" begin include("do_transform_tests.jl") end
@testset "Plotting tests" begin include("do_plotting_tests.jl") end
@testset "Plotting tests" begin include("do_diff_tests.jl") end


