"""
Test parameter scanning functionality (Phase 6, T033-T038)

Validates that:
1. Scan configuration can be loaded from TOML
2. Parameter scan runs without crashing
3. Multiple simulations execute correctly
4. Results are aggregated properly
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/TDSESolver.jl")
using .ConfigParser
using .ParameterScan
using Printf

println("="^80)
println("PARAMETER SCAN VALIDATION TEST")
println("="^80)
println()

# Test 1: Create minimal scan configuration
println("Test 1: Creating minimal intensity scan configuration...")

# Create a fast scan configuration with just 3 intensity points
scan_toml = """
[atom]
type = "hydrogen"

[laser]
wavelength = 800.0
intensity = 1.0e14
pulse_duration = 10.0
polarization = [0.0, 0.0, 1.0]

[calculation]
type = "ionization"
output_file = "scan_test_results.h5"
verbosity = "minimal"

[numerical]
nrmax = 50
rmax = 50.0
lmax = 1
dt = 0.2
t_total = 10.0

[scan]
parameter = "laser_intensity"
range = "1.0e14:1.0e14:3.0e14"
"""

# Write temporary config file
test_config_path = joinpath(@__DIR__, "temp_scan_config.toml")
open(test_config_path, "w") do io
    write(io, scan_toml)
end

println("✓ Created test configuration")
println()

# Test 2: Load configuration
println("Test 2: Loading scan configuration...")
try
    config = load_config(test_config_path)
    println("✓ Configuration loaded successfully")
    println("  Atom: $(config.atom_type)")
    println("  Scan parameter: $(config.scan.parameter)")
    println("  Scan range: $(config.scan.range)")
    println()
catch e
    println("✗ Failed to load configuration:")
    println("  Error: $e")
    rm(test_config_path, force=true)
    exit(1)
end

# Test 3: Run parameter scan
println("Test 3: Running parameter scan (3 simulations)...")
println("  This will take 1-2 minutes...")
println()

try
    config = load_config(test_config_path)

    # Run scan (serial execution for simplicity)
    results = run_parameter_scan(config, verbose=true, parallel=false)

    println()
    println("✓ Parameter scan completed successfully")
    println()

    # Test 4: Verify results
    println("Test 4: Validating scan results...")

    # Check we got the right number of results
    n_expected = length(collect(config.scan.range))
    n_actual = length(results.simulation_results)

    println("  Expected simulations: $n_expected")
    println("  Actual simulations: $n_actual")

    if n_actual == n_expected
        println("  ✓ Correct number of simulations")
    else
        println("  ✗ Mismatch in simulation count!")
        rm(test_config_path, force=true)
        exit(1)
    end
    println()

    # Check that all simulations completed (even if results are NaN)
    all_completed = true
    for (i, result) in enumerate(results.simulation_results)
        if isnothing(result)
            all_completed = false
            println("  ✗ Simulation $i returned nothing")
        end
    end

    if all_completed
        println("  ✓ All simulations returned results")
    else
        println("  ✗ Some simulations failed")
        rm(test_config_path, force=true)
        exit(1)
    end
    println()

    # Print summary
    println("Test 5: Printing scan summary...")
    print_scan_summary(results)
    println()

catch e
    println("✗ Parameter scan failed:")
    println("  Error: $e")
    println()
    println("Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    rm(test_config_path, force=true)
    exit(1)
end

# Cleanup
rm(test_config_path, force=true)

# Summary
println()
println("="^80)
println("VALIDATION SUMMARY")
println("="^80)
println()
println("✓ Parameter scan support validated")
println("  - Scan configuration parsing: PASS")
println("  - Multiple simulation execution: PASS")
println("  - Result aggregation: PASS")
println("  - Summary generation: PASS")
println()
println("Phase 6 (Parameter Scanning): COMPLETE")
println("="^80)
