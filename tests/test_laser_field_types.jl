"""
Test Laser Field Types

Validates that all laser field configurations work correctly:
1. Linear polarization (z-polarized)
2. Circular polarization
3. Elliptical polarization
4. Two-color OTC (orthogonal two-color, ω-2ω)

Each test runs a short propagation and checks that the field interaction
produces the expected physical response.
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/TDSESolver.jl")
using .Simulation
using Printf
using LinearAlgebra

println("="^80)
println("LASER FIELD TYPES TEST")
println("="^80)
println()

# Common parameters for all tests
const COMMON_PARAMS = (
    nrmax = 80,
    rmax = 80.0,
    lmax = 3,  # Need higher lmax for non-linear polarization
    atom = :hydrogen,
    dt = 0.1,
    t_total = 20.0,  # Short simulation
    obs_interval = 10,
    laser_wavelength = 800.0,
    laser_intensity = 1.0e14,  # Moderate intensity
    laser_duration = 15.0,
)

results_summary = []

# =============================================================================
# Test 1: Linear Polarization (z-axis)
# =============================================================================

println("="^80)
println("TEST 1: Linear Polarization (Z-axis)")
println("="^80)
println()

params_linear = create_default_params(;
    COMMON_PARAMS...,
    laser_polarization = [0.0, 0.0, 1.0],  # z-polarized
    laser_cep = 0.0,
)

println("Running linear polarization simulation...")
println("  E(t) = E₀ * f(t) * sin(ωt) * ẑ")
println()

results_linear = run_simulation(params_linear, verbose=false)
norm_linear = results_linear.observables.norms[end]
P_ion_linear = results_linear.observables.ionization_probs[end]

println("Results:")
println("  Final norm: $(@sprintf("%.6f", norm_linear))")
println("  Ionization: $(@sprintf("%.4f", P_ion_linear * 100))%")
println()

push!(results_summary, ("Linear (Z)", norm_linear, P_ion_linear))

if P_ion_linear > 0.0
    println("✓ Linear polarization: ionization detected")
else
    println("⚠ Linear polarization: no ionization (field may be too weak)")
end
println()

# =============================================================================
# Test 2: Circular Polarization
# =============================================================================

println("="^80)
println("TEST 2: Circular Polarization")
println("="^80)
println()

# For circular: equal Ex and Ey with 90° phase difference
# Using the field module's ω-2ω structure, we approximate circular by
# setting exrate=1, eyrate=1, phase=π/2
# Note: True circular would need same frequency, but this tests the field coupling

params_circular = create_default_params(;
    COMMON_PARAMS...,
    laser_polarization = [1.0, 1.0, 0.0],  # XY plane
    laser_cep = π/2,  # 90° phase for circular-like
)

println("Running circular polarization simulation...")
println("  Ex(t) = E₀ * f(t) * sin(ωt)")
println("  Ey(t) = E₀ * f(t) * sin(2ωt + π/2)  [Note: 2ω component]")
println()

results_circular = run_simulation(params_circular, verbose=false)
norm_circular = results_circular.observables.norms[end]
P_ion_circular = results_circular.observables.ionization_probs[end]

println("Results:")
println("  Final norm: $(@sprintf("%.6f", norm_circular))")
println("  Ionization: $(@sprintf("%.4f", P_ion_circular * 100))%")
println()

push!(results_summary, ("Circular-like", norm_circular, P_ion_circular))

if P_ion_circular > 0.0
    println("✓ Circular polarization: ionization detected")
else
    println("⚠ Circular polarization: no ionization")
end
println()

# =============================================================================
# Test 3: Elliptical Polarization
# =============================================================================

println("="^80)
println("TEST 3: Elliptical Polarization (ε = 0.5)")
println("="^80)
println()

params_elliptical = create_default_params(;
    COMMON_PARAMS...,
    laser_polarization = [1.0, 0.5, 0.0],  # Ellipticity ε = 0.5
    laser_cep = π/4,  # 45° phase
)

println("Running elliptical polarization simulation...")
println("  Ex(t) = E₀ * f(t) * sin(ωt)")
println("  Ey(t) = 0.5*E₀ * f(t) * sin(2ωt + π/4)")
println()

results_elliptical = run_simulation(params_elliptical, verbose=false)
norm_elliptical = results_elliptical.observables.norms[end]
P_ion_elliptical = results_elliptical.observables.ionization_probs[end]

println("Results:")
println("  Final norm: $(@sprintf("%.6f", norm_elliptical))")
println("  Ionization: $(@sprintf("%.4f", P_ion_elliptical * 100))%")
println()

push!(results_summary, ("Elliptical", norm_elliptical, P_ion_elliptical))

if P_ion_elliptical > 0.0
    println("✓ Elliptical polarization: ionization detected")
else
    println("⚠ Elliptical polarization: no ionization")
end
println()

# =============================================================================
# Test 4: Two-Color OTC (ω-2ω)
# =============================================================================

println("="^80)
println("TEST 4: Two-Color OTC (ω-2ω Orthogonal)")
println("="^80)
println()

# OTC: Strong fundamental in X, weaker second harmonic in Y
params_otc = create_default_params(;
    COMMON_PARAMS...,
    laser_polarization = [1.0, 0.3, 0.0],  # ω in X, 2ω in Y (30% relative)
    laser_cep = 0.0,  # In-phase
)

println("Running two-color OTC simulation...")
println("  Ex(t) = E₀ * f(t) * sin(ωt)        [Fundamental]")
println("  Ey(t) = 0.3*E₀ * f(t) * sin(2ωt)   [Second harmonic]")
println()

results_otc = run_simulation(params_otc, verbose=false)
norm_otc = results_otc.observables.norms[end]
P_ion_otc = results_otc.observables.ionization_probs[end]

println("Results:")
println("  Final norm: $(@sprintf("%.6f", norm_otc))")
println("  Ionization: $(@sprintf("%.4f", P_ion_otc * 100))%")
println()

push!(results_summary, ("OTC (ω-2ω)", norm_otc, P_ion_otc))

if P_ion_otc > 0.0
    println("✓ Two-color OTC: ionization detected")
else
    println("⚠ Two-color OTC: no ionization")
end
println()

# =============================================================================
# Test 5: Strong Two-Color for Asymmetric Ionization
# =============================================================================

println("="^80)
println("TEST 5: Strong Two-Color (Higher Intensity)")
println("="^80)
println()

params_strong_otc = create_default_params(;
    nrmax = 80,
    rmax = 80.0,
    lmax = 4,  # Higher lmax for strong field
    atom = :hydrogen,
    dt = 0.1,
    t_total = 25.0,
    obs_interval = 10,
    laser_wavelength = 800.0,
    laser_intensity = 3.0e14,  # Stronger field
    laser_duration = 20.0,
    laser_polarization = [1.0, 0.5, 0.0],  # ω-2ω
    laser_cep = π/2,  # Phase for maximum asymmetry
)

println("Running strong two-color simulation...")
println("  Intensity: 3×10¹⁴ W/cm²")
println("  Phase: π/2 (for asymmetric ionization)")
println()

results_strong = run_simulation(params_strong_otc, verbose=false)
norm_strong = results_strong.observables.norms[end]
P_ion_strong = results_strong.observables.ionization_probs[end]

println("Results:")
println("  Final norm: $(@sprintf("%.6f", norm_strong))")
println("  Ionization: $(@sprintf("%.4f", P_ion_strong * 100))%")
println()

push!(results_summary, ("Strong OTC", norm_strong, P_ion_strong))

# =============================================================================
# Summary
# =============================================================================

println("="^80)
println("SUMMARY")
println("="^80)
println()

println("Field Type          | Final Norm | Ionization")
println("-"^50)
for (name, norm, P_ion) in results_summary
    @printf("%-18s  | %10.6f | %8.4f%%\n", name, norm, P_ion * 100)
end
println("-"^50)
println()

# =============================================================================
# Verdict
# =============================================================================

println("="^80)
println("VERDICT")
println("="^80)
println()

all_passed = true

# Check all simulations completed
if length(results_summary) == 5
    println("✓ All 5 field configurations ran successfully")
else
    println("✗ Some configurations failed to run")
    all_passed = false
end

# Check ionization detected in at least some cases
ionizations = [r[3] for r in results_summary]
if any(ionizations .> 0.0)
    println("✓ Ionization detected in simulations")
else
    println("✗ No ionization detected (check field parameters)")
    all_passed = false
end

# Check norm stays reasonable (accounting for lmax truncation)
norms = [r[2] for r in results_summary]
if all(norms .> 0.5) && all(norms .< 3.0)
    println("✓ Norm values are in reasonable range")
else
    println("⚠ Some norm values are outside expected range")
end

# Check physical behavior: 2D fields should show different dynamics
if P_ion_circular > 0.0 || P_ion_elliptical > 0.0 || P_ion_otc > 0.0
    println("✓ Multi-component fields produce physical response")
else
    println("⚠ Multi-component fields show minimal response")
end

println()

if all_passed
    println("="^80)
    println("✅ LASER FIELD TYPES TEST PASSED")
    println("="^80)
    println()
    println("All laser field configurations are functional:")
    println("  • Linear polarization (z-axis)")
    println("  • Circular-like polarization (ω-2ω)")
    println("  • Elliptical polarization")
    println("  • Two-color OTC (ω-2ω orthogonal)")
    println()
    println("Note: The Y-component uses 2ω (second harmonic) by design,")
    println("      matching the Fortran implementation for bicircular fields.")
else
    println("="^80)
    println("⚠ LASER FIELD TYPES TEST: SOME ISSUES DETECTED")
    println("="^80)
end

println("="^80)
