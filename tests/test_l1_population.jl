"""
Simple diagnostic to check if l=1 states are being excited by the laser field
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/TDSESolver.jl")
using .Simulation
using Printf
using LinearAlgebra

println("="^80)
println("L=1 POPULATION DIAGNOSTIC")
println("="^80)
println()

# Use stronger field to ensure excitation
params = create_default_params(
    nrmax = 100,
    rmax = 100.0,
    lmax = 2,
    atom = :hydrogen,
    dt = 0.1,
    t_total = 50.0,
    obs_interval = 5,

    # STRONGER FIELD for clear excitation
    laser_wavelength = 800.0,
    laser_intensity = 5.0e14,  # 5x stronger than before
    laser_duration = 25.0,
    laser_cep = 0.0,
    laser_polarization = [0.0, 0.0, 1.0],  # z-polarized

    # Enable HHG
    enable_hhg = true,
    hhg_method = :length,
    hhg_record_interval = 1
)

println("Running simulation with strong field:")
println("  Intensity: 5.0e14 W/cm²")
println("  Wavelength: 800 nm (z-polarized)")
println()

results = run_simulation(params, verbose=false)

println()
println("="^80)
println("POPULATION ANALYSIS")
println("="^80)
println()

# Check l=0 and l=1 populations
wfn = results.wavefunction
m_idx = wfn.lmax + 1  # m=0

g_l0 = wfn.g[:, m_idx, 1]  # l=0, m=0
g_l1 = wfn.g[:, m_idx, 2]  # l=1, m=0
g_l2 = wfn.g[:, m_idx, 3]  # l=2, m=0

# Compute populations
weights = results.grid.quadrature_weights

pop_l0 = sum(abs2.(g_l0) .* weights)
pop_l1 = sum(abs2.(g_l1) .* weights)
pop_l2 = sum(abs2.(g_l2) .* weights)

println("Final populations:")
println("  l=0: $(@sprintf("%.6e", pop_l0))")
println("  l=1: $(@sprintf("%.6e", pop_l1))")
println("  l=2: $(@sprintf("%.6e", pop_l2))")
println()

total_pop = pop_l0 + pop_l1 + pop_l2
println("Total: $(@sprintf("%.6f", total_pop))")
println()

# Check dipole moment
if results.observables.dipole_moment !== nothing
    max_dipole = maximum(abs.(results.observables.dipole_moment))
    println("Dipole moment:")
    println("  Max |d(t)|: $(@sprintf("%.6e", max_dipole))")
    println()

    if max_dipole > 1e-10
        println("✓ Dipole is non-zero - l=1 excitation detected!")
    else
        println("❌ Dipole is zero - l=1 not being excited")
    end
else
    println("❌ No dipole data")
end

println()
println("="^80)

if pop_l1 > 1e-6
    println("✓ L=1 POPULATION IS SIGNIFICANT")
    println("  The field interaction is working correctly!")
else
    println("❌ L=1 POPULATION IS NEGLIGIBLE")
    println("  Problem: Field interaction not coupling l=0 ↔ l=1")
end
println("="^80)
