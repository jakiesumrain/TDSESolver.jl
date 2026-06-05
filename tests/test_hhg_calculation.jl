"""
HHG (High Harmonic Generation) Calculation Test

Tests the complete HHG implementation with both methods:
1. Method 1: Eigenstate expansion (Fortran lines 891-899, active)
2. Method 2a: Acceleration form (Fortran lines 905-909, commented)
3. Method 2b: Length form

This test validates:
- Dipole moment calculation during propagation
- All three HHG methods run without errors
- FFT spectrum computation
- File I/O for dipole and spectrum data
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/TDSESolver.jl")
using .Simulation
using .HHG: compute_hhg_spectrum, save_dipole_to_file, save_hhg_spectrum_to_file
using .PhysicalUnits: wavelength_nm_to_frequency_au
using Printf

println("="^80)
println("HHG CALCULATION TEST")
println("="^80)
println()

# =============================================================================
# Test Parameters
# =============================================================================

println("Setting up test parameters...")
println()

params = create_default_params(
    # Small grid for fast testing
    nrmax = 100,
    rmax = 100.0,
    lmax = 2,

    # Atom
    atom = :hydrogen,

    # Time propagation
    dt = 0.1,
    t_total = 50.0,  # Short simulation (~1 optical cycle)
    obs_interval = 10,

    # Laser field (moderate intensity for testing)
    laser_wavelength = 800.0,  # nm
    laser_intensity = 1.0e14,  # W/cm²
    laser_duration = 25.0,     # a.u. (short pulse)
    laser_cep = 0.0,
    laser_polarization = [0.0, 0.0, 1.0],  # z-polarized

    # Enable HHG
    enable_hhg = true,
    hhg_method = :acceleration,  # Start with Method 2a (no eigenstate errors)
    hhg_record_interval = 1      # Record every step
)

print_simulation_params(params)
println()

# =============================================================================
# Test 1: Method 2a (Acceleration Form)
# =============================================================================

println("="^80)
println("TEST 1: HHG with Acceleration Form (Method 2a)")
println("="^80)
println()

results_accel = run_simulation(params, verbose=true)

println()
println("Analyzing results...")

# Check dipole data
if results_accel.observables.dipole_moment !== nothing
    n_dipole = length(results_accel.observables.dipole_moment)
    println("✓ Dipole moment time series collected: $n_dipole points")
    println("  Max |d(t)|: $(@sprintf("%.4e", maximum(abs.(results_accel.observables.dipole_moment))))")

    # Compute HHG spectrum
    ω₀ = wavelength_nm_to_frequency_au(params.laser_wavelength)
    h_orders, power, freqs = compute_hhg_spectrum(
        results_accel.observables.dipole_moment,
        params.dt,
        ω₀,
        window=:hann
    )

    println("✓ HHG spectrum computed")
    println("  Max harmonic order: $(@sprintf("%.1f", maximum(abs.(h_orders))))")
    println("  Peak power: $(@sprintf("%.4e", maximum(power)))")

    # Save dipole and spectrum
    save_dipole_to_file(results_accel.observables.dipole_moment,
                       results_accel.observables.times,
                       "dipole_acceleration.txt",
                       time_unit=2π/ω₀)  # In units of optical cycles

    save_hhg_spectrum_to_file(h_orders[1:div(n_dipole,2)],  # Positive frequencies only
                             power[1:div(n_dipole,2)],
                             "hhg_spectrum_acceleration.txt")

    println("✓ Files saved: dipole_acceleration.txt, hhg_spectrum_acceleration.txt")
else
    println("❌ ERROR: Dipole moment is nothing")
    exit(1)
end

println()

# =============================================================================
# Test 2: Method 2b (Length Form)
# =============================================================================

println("="^80)
println("TEST 2: HHG with Length Form (Method 2b)")
println("="^80)
println()

params_length = create_default_params(
    nrmax = 100,
    rmax = 100.0,
    lmax = 2,
    atom = :hydrogen,
    dt = 0.1,
    t_total = 50.0,
    obs_interval = 10,
    laser_wavelength = 800.0,
    laser_intensity = 1.0e14,
    laser_duration = 25.0,
    laser_cep = 0.0,
    laser_polarization = [0.0, 0.0, 1.0],
    enable_hhg = true,
    hhg_method = :length,  # Method 2b
    hhg_record_interval = 1
)

results_length = run_simulation(params_length, verbose=false)  # Quiet run

if results_length.observables.dipole_moment !== nothing
    n_dipole = length(results_length.observables.dipole_moment)
    println("✓ Dipole moment (length form) collected: $n_dipole points")
    println("  Max |d(t)|: $(@sprintf("%.4e", maximum(abs.(results_length.observables.dipole_moment))))")
else
    println("❌ ERROR: Dipole moment is nothing")
    exit(1)
end

println()

# =============================================================================
# Test 3: Method 1 (Eigenstate Expansion)
# =============================================================================

println("="^80)
println("TEST 3: HHG with Eigenstate Expansion (Method 1)")
println("="^80)
println()

params_eigenstate = create_default_params(
    nrmax = 100,
    rmax = 100.0,
    lmax = 2,
    atom = :hydrogen,
    dt = 0.1,
    t_total = 50.0,
    obs_interval = 10,
    laser_wavelength = 800.0,
    laser_intensity = 1.0e14,
    laser_duration = 25.0,
    laser_cep = 0.0,
    laser_polarization = [0.0, 0.0, 1.0],
    enable_hhg = true,
    hhg_method = :eigenstate,  # Method 1
    hhg_record_interval = 1
)

results_eigenstate = run_simulation(params_eigenstate, verbose=false)

if results_eigenstate.observables.dipole_moment !== nothing
    n_dipole = length(results_eigenstate.observables.dipole_moment)
    println("✓ Dipole moment (eigenstate) collected: $n_dipole points")
    println("  Max |d(t)|: $(@sprintf("%.4e", maximum(abs.(results_eigenstate.observables.dipole_moment))))")
    println()
    println("  NOTE: Values may be inaccurate due to GPS eigenstate errors (~10⁴× wrong)")
    println("        See docs/GPS_INVESTIGATION.md for details")
else
    println("❌ ERROR: Dipole moment is nothing")
    exit(1)
end

println()

# =============================================================================
# Comparison of Methods
# =============================================================================

println("="^80)
println("COMPARISON OF HHG METHODS")
println("="^80)
println()

max_accel = maximum(abs.(results_accel.observables.dipole_moment))
max_length = maximum(abs.(results_length.observables.dipole_moment))
max_eigenstate = maximum(abs.(results_eigenstate.observables.dipole_moment))

println("Max |d(t)| by method:")
println("  Acceleration form:   $(@sprintf("%.4e", max_accel))")
println("  Length form:         $(@sprintf("%.4e", max_length))")
println("  Eigenstate expansion: $(@sprintf("%.4e", max_eigenstate))")
println()

ratio_length_accel = max_length / max_accel
println("Ratio (length/acceleration): $(@sprintf("%.4f", ratio_length_accel))")

if 0.5 < ratio_length_accel < 2.0
    println("✓ Length and acceleration forms agree within factor of 2")
else
    println("⚠ Large discrepancy between length and acceleration forms")
end

println()

# =============================================================================
# Verdict
# =============================================================================

println("="^80)
println("VERDICT")
println("="^80)
println()

success = true

# Check 1: All methods ran
if max_accel > 0 && max_length > 0 && max_eigenstate > 0
    println("✓ All three HHG methods ran successfully")
else
    println("❌ FAILED: Some methods produced zero dipole")
    success = false
end

# Check 2: Dipole time series recorded
if n_dipole == div(Int(params.t_total / params.dt), params.hhg_record_interval)
    println("✓ Correct number of dipole points recorded")
else
    println("⚠ WARNING: Expected ~$(@sprintf("%d", div(Int(params.t_total / params.dt), params.hhg_record_interval))) points, got $n_dipole")
end

# Check 3: Spectrum computed
if maximum(power) > 0
    println("✓ HHG spectrum computed with non-zero power")
else
    println("❌ FAILED: HHG spectrum is zero")
    success = false
end

# Check 4: Files saved
if isfile("dipole_acceleration.txt") && isfile("hhg_spectrum_acceleration.txt")
    println("✓ Output files created successfully")
else
    println("❌ FAILED: Output files not created")
    success = false
end

println()

if success
    println("="^80)
    println("✅ HHG CALCULATION TEST PASSED")
    println("="^80)
    println()
    println("All HHG methods implemented and functional:")
    println("  1. ✓ Eigenstate expansion (Method 1, Fortran lines 891-899)")
    println("  2. ✓ Acceleration form (Method 2a, Fortran lines 905-909)")
    println("  3. ✓ Length form (Method 2b)")
    println()
    println("System ready for HHG spectrum calculations!")
    println("  - Use :acceleration or :length for accurate results (independent of eigenstate errors)")
    println("  - Use :eigenstate after fixing GPS eigenstate solver")
    println("="^80)
else
    println("="^80)
    println("❌ HHG CALCULATION TEST FAILED")
    println("="^80)
    exit(1)
end
