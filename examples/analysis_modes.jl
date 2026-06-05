"""
Quick Reference: Analysis Mode Examples

This file shows the exact parameter settings for each of the 4 operational modes.
Copy-paste the mode you need!
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/TDSESolver.jl")
using .Simulation

# =============================================================================
# Mode 1: Ionization Only (FASTEST - ~10 sec)
# =============================================================================
# Use for: Quick tests, parameter sweeps, basic ionization probability
# =============================================================================

params_mode1 = create_default_params(
    # Basic setup
    atom = :hydrogen,
    nrmax = 200,
    rmax = 100.0,
    lmax = 2,
    dt = 0.1,
    t_total = 100.0,

    # Laser
    laser_wavelength = 800.0,
    laser_intensity = 1.0e14,
    laser_duration = 50.0,

    # BOTH SWITCHES OFF
    enable_momentum_analysis = false,  # ← OFF
    enable_hhg = false                 # ← OFF
)

# results1 = run_simulation(params_mode1)
# Output: ionization_probs, norms, energies

# =============================================================================
# Mode 2: Ionization + Momentum Distributions (~20 sec)
# =============================================================================
# Use for: Photoelectron spectra, ATI peaks, angular distributions
# =============================================================================

params_mode2 = create_default_params(
    # Basic setup
    atom = :hydrogen,
    nrmax = 200,
    rmax = 100.0,
    lmax = 2,
    dt = 0.1,
    t_total = 100.0,

    # Laser (higher intensity for ionization)
    laser_wavelength = 800.0,
    laser_intensity = 5.0e14,  # ← Higher for ionization
    laser_duration = 50.0,

    # MOMENTUM ON, HHG OFF
    enable_momentum_analysis = true,   # ← ON
    msplit = 50,
    R_c = 100.0,
    delta_split = 5.0,
    p_max = 2.0,
    n_p = 30,
    n_theta = 20,
    n_phi = 20,

    enable_hhg = false                 # ← OFF
)

# results2 = run_simulation(params_mode2)
# Output: ionization_probs + momentum_grid + psai_p
#
# Post-process momentum distributions:
# using .MomentumDistribution
# rate_xyz = transform_to_cartesian_grid(results2.psai_p, results2.momentum_grid, p_max, 30)
# pxy_rate, pxz_rate, pyz_rate = compute_integrated_2d_distributions(rate_xyz, p_max, 30)

# =============================================================================
# Mode 3: Ionization + HHG Spectrum (~12 sec)
# =============================================================================
# Use for: Harmonic generation, attosecond physics, spectral analysis
# =============================================================================

params_mode3 = create_default_params(
    # Basic setup
    atom = :hydrogen,
    nrmax = 200,
    rmax = 100.0,
    lmax = 2,
    dt = 0.1,
    t_total = 100.0,

    # Laser
    laser_wavelength = 800.0,
    laser_intensity = 1.0e14,
    laser_duration = 50.0,

    # MOMENTUM OFF, HHG ON
    enable_momentum_analysis = false,  # ← OFF

    enable_hhg = true,                 # ← ON
    hhg_method = :acceleration,        # :eigenstate, :acceleration, or :length
    hhg_record_interval = 1            # Every step for accurate spectrum
)

# results3 = run_simulation(params_mode3)
# Output: ionization_probs + dipole_moment
#
# Post-process HHG spectrum:
# using .HHG
# using .PhysicalUnits: wavelength_nm_to_frequency_au
# ω₀ = wavelength_nm_to_frequency_au(params_mode3.laser_wavelength)
# h_orders, power, freqs = compute_hhg_spectrum(results3.observables.dipole_moment,
#                                                params_mode3.dt, ω₀, window=:hann)
# save_hhg_spectrum_to_file(h_orders[1:div(end,2)], power[1:div(end,2)], "hhg_spectrum.txt")

# =============================================================================
# Mode 4: Full Analysis - Ionization + Momentum + HHG (~25 sec)
# =============================================================================
# Use for: Complete study, publication-quality results, all physics
# =============================================================================

params_mode4 = create_default_params(
    # Basic setup
    atom = :hydrogen,
    nrmax = 200,
    rmax = 100.0,
    lmax = 2,
    dt = 0.1,
    t_total = 100.0,

    # Laser (high intensity for both ionization and HHG)
    laser_wavelength = 800.0,
    laser_intensity = 5.0e14,
    laser_duration = 50.0,

    # BOTH SWITCHES ON
    enable_momentum_analysis = true,   # ← ON
    msplit = 50,
    R_c = 100.0,
    delta_split = 5.0,
    p_max = 2.0,
    n_p = 30,
    n_theta = 20,
    n_phi = 20,

    enable_hhg = true,                 # ← ON
    hhg_method = :acceleration,
    hhg_record_interval = 1
)

# results4 = run_simulation(params_mode4)
# Output: ALL - ionization_probs + momentum_grid + psai_p + dipole_moment
#
# Can do BOTH post-processing analyses!

# =============================================================================
# Quick Mode Selector
# =============================================================================

println("Analysis Mode Examples")
println("="^80)
println()
println("Mode 1: Ionization only")
println("  - Switches: enable_momentum_analysis=false, enable_hhg=false")
println("  - Use: Quick tests, parameter sweeps")
println()
println("Mode 2: Ionization + Momentum")
println("  - Switches: enable_momentum_analysis=true, enable_hhg=false")
println("  - Use: Photoelectron spectra, ATI analysis")
println()
println("Mode 3: Ionization + HHG")
println("  - Switches: enable_momentum_analysis=false, enable_hhg=true")
println("  - Use: Harmonic generation, attosecond physics")
println()
println("Mode 4: Full Analysis")
println("  - Switches: enable_momentum_analysis=true, enable_hhg=true")
println("  - Use: Complete study, all physics")
println()
println("="^80)
println()
println("Uncomment the mode you want to run and execute:")
println("  julia examples/analysis_modes.jl")
println()

# Uncomment ONE of these to run:
# results = run_simulation(params_mode1)  # Ionization only
# results = run_simulation(params_mode2)  # + Momentum
# results = run_simulation(params_mode3)  # + HHG
# results = run_simulation(params_mode4)  # Full analysis
