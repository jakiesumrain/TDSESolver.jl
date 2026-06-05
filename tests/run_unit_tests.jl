"""
Master test runner for all unit tests.

Runs unit tests for all foundational modules and reports results.
"""

using Test

println("=" ^ 80)
println("Running Unit Tests for TDSE Solver Foundational Modules")
println("=" ^ 80)
println()

# Track test results
test_files = [
    "unit/test_physical_units.jl",
    "unit/test_config_parser.jl",
    "unit/test_validation.jl",
    "unit/test_potential.jl",
    "unit/test_gps_grid.jl",
    "unit/test_angular_grid.jl",
    "unit/test_results_io.jl"
]

module_names = [
    "PhysicalUnits",
    "ConfigParser",
    "Validation",
    "Potential",
    "GPSGrid",
    "AngularGrid",
    "ResultsIO"
]

passed = 0
failed = 0
errors = []

println("Testing $(length(test_files)) foundational modules:\n")

for (i, test_file) in enumerate(test_files)
    module_name = module_names[i]
    print("[$i/$(length(test_files))] Testing $module_name module... ")

    try
        include(test_file)
        println("✓ PASSED")
        global passed += 1
    catch e
        println("✗ FAILED")
        global failed += 1
        push!(errors, (module_name, e))
    end
end

println()
println("=" ^ 80)
println("Test Summary")
println("=" ^ 80)
println("Passed: $passed / $(length(test_files))")
println("Failed: $failed / $(length(test_files))")

if failed > 0
    println()
    println("Failed Tests:")
    for (module_name, error) in errors
        println("  - $module_name:")
        println("    $(sprint(showerror, error))")
    end
    println()
    error("❌ Some tests failed!")
else
    println()
    println("✓ All foundational module tests passed!")
    println()
    println("Phase 2 (Foundational Infrastructure) is validated and ready.")
    println("You can now proceed to Phase 3 (User Story 1 - MVP Implementation).")
end

println("=" ^ 80)
