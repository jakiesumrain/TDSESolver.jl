"""
Validate that all example configurations can be loaded successfully
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/TDSESolver.jl")
using .ConfigParser

println("="^80)
println("EXAMPLE CONFIGURATION VALIDATION")
println("="^80)
println()

examples_dir = joinpath(@__DIR__, "..", "examples")
toml_files = filter(f -> endswith(f, ".toml"), readdir(examples_dir))

println("Found $(length(toml_files)) example configurations:")
for file in toml_files
    println("  - $file")
end
println()

all_valid = true

for file in toml_files
    print("Testing $file... ")
    try
        config = load_config(joinpath(examples_dir, file))
        println("✓ PASS")

        # Print brief summary
        println("  Atom: $(config.atom_type)")
        println("  Wavelength: $(config.wavelength_nm) nm")
        println("  Intensity: $(config.intensity_W_cm2) W/cm²")
        if !isnothing(config.scan)
            println("  Scan: $(config.scan.parameter)")
        end
        println()
    catch e
        println("✗ FAIL")
        println("  Error: $e")
        println()
        all_valid = false
    end
end

println("="^80)
if all_valid
    println("✓ All example configurations are valid")
else
    println("✗ Some configurations failed validation")
    exit(1)
end
println("="^80)
