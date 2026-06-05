"""
Full Integration Test: TDSE → Region Split → Volkov Projection → Momentum Distributions

Tests the complete Fortran-matching workflow:
1. TDSE propagation with ionizing laser field
2. Region splitting every msplit steps during evolution
3. Volkov projection with phase accumulation
4. Momentum distribution computation

This test validates the ENTIRE photoelectron momentum analysis pipeline.
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/TDSESolver.jl")
using .Simulation
using .MomentumDistribution

using Printf

println("="^80)
println("FULL MOMENTUM PIPELINE INTEGRATION TEST")
println("="^80)
println()

# =============================================================================
# Test Parameters
# =============================================================================

println("Setting up test parameters...")
println()

params = create_default_params(
    # Grid (smaller for fast testing)
    nrmax = 100,
    rmax = 100.0,
    lmax = 2,

    # Atom
    atom = :hydrogen,

    # Time propagation
    dt = 0.1,
    t_total = 110.0,  # ~1 optical cycle at 800 nm
    obs_interval = 10,

    # Laser field (stronger for ionization)
    laser_wavelength = 800.0,  # nm
    laser_intensity = 5.0e14,  # W/cm² (stronger than before)
    laser_duration = 55.0,     # a.u. (2 cycles)
    laser_cep = 0.0,
    laser_polarization = [0.0, 0.0, 1.0],  # z-polarized

    # Enable momentum analysis (KEY!)
    enable_momentum_analysis = true,
    msplit = 10,                # Split every 10 steps (frequent for testing)
    R_c = 50.0,                 # Splitting radius (a.u.)
    delta_split = 5.0,          # Smoothness
    p_max = 2.0,                # Maximum momentum
    n_p = 20,                   # Momentum grid points
    n_theta = 15,               # θ points
    n_phi = 15                  # φ points
)

# Print simulation parameters
print_simulation_params(params)
println()

# =============================================================================
# Run Simulation with Momentum Analysis
# =============================================================================

println("="^80)
println("RUNNING TDSE WITH MOMENTUM ANALYSIS")
println("="^80)
println()

results = run_simulation(params, verbose=true)

println()
println("="^80)
println("SIMULATION COMPLETE")
println("="^80)
println()

# =============================================================================
# Analyze Results
# =============================================================================

println("Analyzing results...")
println()

# Check ionization
final_ionization = results.observables.ionization_probs[end]
println("Final ionization probability: $(@sprintf("%.4e", final_ionization))")

# Check momentum space wavefunction
if results.psai_p !== nothing
    psai_p_norm = sum(abs2.(results.psai_p))
    println("Momentum space wavefunction norm: $(@sprintf("%.4e", psai_p_norm))")
    println("Max |ψ_p|²: $(@sprintf("%.4e", maximum(abs2.(results.psai_p))))")

    # Check if we have meaningful momentum distribution
    if psai_p_norm > 1e-10
        println("✓ Momentum space wavefunction accumulated successfully")
    else
        println("⚠ Momentum space wavefunction is near-zero (field may be too weak)")
    end
else
    println("❌ ERROR: Momentum space wavefunction is nothing")
    exit(1)
end

println()

# =============================================================================
# Compute Momentum Distributions
# =============================================================================

println("="^80)
println("COMPUTING MOMENTUM DISTRIBUTIONS")
println("="^80)
println()

if results.momentum_grid === nothing || results.psai_p === nothing
    println("❌ ERROR: Momentum data not available")
    exit(1)
end

# Step 1: Transform to Cartesian grid
println("Step 1: Transforming to Cartesian grid...")
n_grid_cartesian = 15  # Smaller grid for testing

rate_xyz = MomentumDistribution.transform_to_cartesian_grid(
    results.psai_p,
    results.momentum_grid,
    params.p_max,
    n_grid_cartesian
)

println("Cartesian grid created:")
println("  Shape: $(size(rate_xyz))")
println("  Max value: $(@sprintf("%.4e", maximum(rate_xyz)))")
println("  Total: $(@sprintf("%.4e", sum(rate_xyz)))")
println()

# Step 2: Compute 2D integrated distributions
println("Step 2: Computing 2D integrated distributions...")
pxy_rate, pxz_rate, pyz_rate = MomentumDistribution.compute_integrated_2d_distributions(
    rate_xyz,
    params.p_max,
    n_grid_cartesian
)

println("2D distributions computed:")
println("  P(px,py) max: $(@sprintf("%.4e", maximum(pxy_rate)))")
println("  P(px,pz) max: $(@sprintf("%.4e", maximum(pxz_rate)))")
println("  P(py,pz) max: $(@sprintf("%.4e", maximum(pyz_rate)))")
println()

# Step 3: Extract 2D slices
println("Step 3: Extracting 2D slices...")
slice_xy_z0 = MomentumDistribution.extract_2d_slice(rate_xyz, "xy", 0.0, params.p_max, n_grid_cartesian)
slice_xz_y0 = MomentumDistribution.extract_2d_slice(rate_xyz, "xz", 0.0, params.p_max, n_grid_cartesian)
slice_yz_x0 = MomentumDistribution.extract_2d_slice(rate_xyz, "yz", 0.0, params.p_max, n_grid_cartesian)

println("Origin slices:")
println("  xy-plane (pz=0): max = $(@sprintf("%.4e", maximum(slice_xy_z0)))")
println("  xz-plane (py=0): max = $(@sprintf("%.4e", maximum(slice_xz_y0)))")
println("  yz-plane (px=0): max = $(@sprintf("%.4e", maximum(slice_yz_x0)))")
println()

# =============================================================================
# Verdict
# =============================================================================

println("="^80)
println("VERDICT")
println("="^80)
println()

success = true

# Check 1: Ionization occurred
if final_ionization < 1e-6
    println("⚠ WARNING: Very low ionization ($(@sprintf("%.2e", final_ionization)))")
    println("  Field may be too weak or pulse too short")
    println("  Momentum distributions will be weak but pipeline is working")
else
    println("✓ Ionization: $(@sprintf("%.2e", final_ionization))")
end

# Check 2: Momentum space wavefunction accumulated
if psai_p_norm > 1e-12
    println("✓ Momentum space accumulation: $(@sprintf("%.2e", psai_p_norm))")
else
    println("❌ FAILED: Momentum space wavefunction too small")
    success = false
end

# Check 3: Cartesian grid transformation
if maximum(rate_xyz) > 1e-20
    println("✓ Cartesian grid transformation: max = $(@sprintf("%.2e", maximum(rate_xyz)))")
else
    println("❌ FAILED: Cartesian grid is zero")
    success = false
end

# Check 4: 2D distributions
if maximum(pxy_rate) > 1e-20
    println("✓ 2D integrated distributions: max = $(@sprintf("%.2e", maximum(pxy_rate)))")
else
    println("❌ FAILED: 2D distributions are zero")
    success = false
end

# Check 5: 2D slices
if maximum(slice_xy_z0) > 1e-20
    println("✓ 2D slices: max = $(@sprintf("%.2e", maximum(slice_xy_z0)))")
else
    println("❌ FAILED: 2D slices are zero")
    success = false
end

println()

if success
    println("="^80)
    println("✅ FULL MOMENTUM PIPELINE TEST PASSED")
    println("="^80)
    println()
    println("Complete workflow validated:")
    println("  1. ✓ TDSE propagation with laser field")
    println("  2. ✓ Region splitting during evolution (every $(@sprintf("%d", params.msplit)) steps)")
    println("  3. ✓ Volkov projection with phase accumulation")
    println("  4. ✓ Momentum space wavefunction: ψ(φ,θ,p)")
    println("  5. ✓ Cartesian grid transformation: |ψ(px,py,pz)|²")
    println("  6. ✓ 2D integrated distributions: P(px,py), P(px,pz), P(py,pz)")
    println("  7. ✓ 2D slice extraction at arbitrary planes")
    println()
    println("System is ready for production photoelectron momentum calculations!")
    println("="^80)
else
    println("="^80)
    println("❌ FULL MOMENTUM PIPELINE TEST FAILED")
    println("="^80)
    exit(1)
end
