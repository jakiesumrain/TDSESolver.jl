"""
Test Region Splitting Module

Validates smooth wavefunction decomposition into inner/outer regions.

Tests:
1. Splitting function properties (sigmoid shape, smoothness)
2. Conservation: ψ_inner + ψ_outer = ψ_total
3. Population conservation: P_inner + P_outer = P_total
4. Smooth transition region
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include("../src/grid/GPSGrid.jl")
using .GPSGrid

include("../src/wavefunction/Wavefunction.jl")
using .Wavefunction

include("../src/hamiltonian/Potential.jl")
using .Potential

include("../src/hamiltonian/Hamiltonian.jl")
using .Hamiltonian

include("../src/region_split/RegionSplit.jl")
using .RegionSplit

using LinearAlgebra
using Printf

println("="^80)
println("Region Splitting Test")
println("="^80)
println()

# =============================================================================
# 1. Setup grid and wavefunction
# =============================================================================

println("STEP 1: Setup")
println("-"^80)

# Radial grid
nrmax = 100
rmax = 150.0
L = 25.0
# Note: α is computed automatically as α = 2L/rmax in create_gps_grid
lmax = 2

println("Creating GPS grid...")
radial_grid = create_gps_grid(nrmax, rmax, L=L)
actual_rmax = maximum(radial_grid.radial_grid)
println("✓ Radial grid: nrmax=$nrmax, rmax=$(round(actual_rmax, digits=2)) a.u.")

# Create wavefunction with hydrogen ground state
pot = get_potential(:hydrogen)
ham = solve_eigenstates(radial_grid, pot, lmax, n_max=nrmax, E_cutoff=Inf)

wfn = Wavefunction.create_wavefunction(radial_grid.nrmax, lmax, radial_grid.quadrature_weights)
φ_ground = ham.eigenvectors[:, 1, 1]
Wavefunction.initialize_ground_state!(wfn, φ_ground, 1, 0)

initial_norm = Wavefunction.compute_norm(wfn)
println("✓ Initial wavefunction: ‖ψ₀‖ = $(@sprintf("%.12f", initial_norm))")
println()

# =============================================================================
# 2. Create region splitter
# =============================================================================

println("="^80)
println("STEP 2: Create region splitter")
println("-"^80)

# Use 70% of actual rmax for R_c (safe within grid)
R_c = 0.7 * actual_rmax
delta = 5.0

println("Parameters:")
println("  R_c = $(@sprintf("%.1f", R_c)) a.u. (critical radius)")
println("  Δ = $delta a.u. (smoothness parameter)")
println("  Transition width ≈ $(4*delta) a.u.")
println()

splitter = create_region_splitter(radial_grid, R_c, delta)
println("✓ Splitter created")
println()

# =============================================================================
# 3. Test splitting function properties
# =============================================================================

println("="^80)
println("STEP 3: Validate splitting function")
println("-"^80)

f_split = get_split_function(splitter)

# Find indices for key regions
r = radial_grid.radial_grid
idx_inner = findfirst(x -> x > R_c - 2*delta, r)
idx_center = findfirst(x -> x > R_c, r)
idx_outer = findfirst(x -> x > R_c + 2*delta, r)

println("Splitting function values:")
println("  f_split(r=0) = $(@sprintf("%.6f", f_split[1])) (should be ≈0)")
println("  f_split(R_c - 2Δ) = $(@sprintf("%.6f", f_split[idx_inner])) (should be ≈0.12)")
println("  f_split(R_c) = $(@sprintf("%.6f", f_split[idx_center])) (should be ≈0.5)")
println("  f_split(R_c + 2Δ) = $(@sprintf("%.6f", f_split[idx_outer])) (should be ≈0.88)")
println("  f_split(rmax) = $(@sprintf("%.6f", f_split[end])) (should be ≈1)")
println()

# Check monotonicity
is_monotonic = all(diff(f_split) .>= 0)
println("Monotonicity check: $(is_monotonic ? "✓ PASS" : "✗ FAIL")")
println()

# =============================================================================
# 4. Split wavefunction
# =============================================================================

println("="^80)
println("STEP 4: Split wavefunction")
println("-"^80)

# Create outer wavefunction
wfn_outer = Wavefunction.create_wavefunction(radial_grid.nrmax, lmax, radial_grid.quadrature_weights)

# Perform splitting
split_wavefunction!(wfn_outer, wfn, splitter)
println("✓ Wavefunction split")
println()

# =============================================================================
# 5. Verify conservation
# =============================================================================

println("="^80)
println("STEP 5: Verify conservation")
println("-"^80)

# Compute populations
pop_inner, pop_outer = compute_split_populations(wfn, wfn_outer, radial_grid)
pop_total = pop_inner + pop_outer

println("Populations:")
println("  Inner region: $(@sprintf("%.8f", pop_inner)) ($(@sprintf("%.2f", 100*pop_inner))%)")
println("  Outer region: $(@sprintf("%.8f", pop_outer)) ($(@sprintf("%.2f", 100*pop_outer))%)")
println("  Total: $(@sprintf("%.12f", pop_total))")
println()

# Conservation check
conservation_error = abs(pop_total - initial_norm)
println("Conservation:")
println("  Initial norm: $(@sprintf("%.12f", initial_norm))")
println("  Split total: $(@sprintf("%.12f", pop_total))")
println("  Error: $(@sprintf("%.2e", conservation_error))")
println()

# =============================================================================
# 6. Verify ψ_inner + ψ_outer = ψ_total
# =============================================================================

println("="^80)
println("STEP 6: Verify ψ_inner + ψ_outer = ψ_total")
println("-"^80)

# Compute ψ_inner = ψ_total - ψ_outer
global max_diff = 0.0
for l in 0:lmax
    for m in -l:l
        m_idx = m + lmax + 1
        for ir in 1:nrmax
            # Inner component
            psi_inner = wfn.g[ir, m_idx, l+1] - wfn_outer.g[ir, m_idx, l+1]

            # Expected inner: (1 - f_split) * psi_total
            psi_inner_expected = (1.0 - f_split[ir]) * wfn.g[ir, m_idx, l+1]

            diff = abs(psi_inner - psi_inner_expected)
            global max_diff = max(max_diff, diff)
        end
    end
end

println("Inner region verification:")
println("  Max |ψ_inner - (1-f)ψ_total|: $(@sprintf("%.2e", max_diff))")
println()

# =============================================================================
# 7. Verdict
# =============================================================================

println("="^80)
println("VERDICT")
println("="^80)
println()

success = true

# Check splitting function
if f_split[1] < 0.01 && f_split[end] > 0.99
    println("✅ Splitting function: correct asymptotic behavior")
else
    println("❌ Splitting function: incorrect asymptotic behavior")
    success = false
end

if abs(f_split[idx_center] - 0.5) < 0.10
    println("✅ Splitting function: f(R_c) ≈ 0.5 (discrete grid)")
else
    println("❌ Splitting function: f(R_c) ≠ 0.5")
    success = false
end

if is_monotonic
    println("✅ Splitting function: monotonically increasing")
else
    println("❌ Splitting function: not monotonic")
    success = false
end

# Check conservation
if conservation_error < 1e-10
    println("✅ EXCELLENT: Population conservation < 1e-10")
elseif conservation_error < 1e-8
    println("✓ GOOD: Population conservation < 1e-8")
else
    println("❌ POOR: Population conservation error > 1e-8")
    success = false
end

# Check decomposition
if max_diff < 1e-14
    println("✅ EXCELLENT: ψ_inner + ψ_outer = ψ_total (machine precision)")
elseif max_diff < 1e-10
    println("✓ GOOD: ψ_inner + ψ_outer ≈ ψ_total")
else
    println("❌ POOR: ψ_inner + ψ_outer ≠ ψ_total")
    success = false
end

println()

if success
    println("="^80)
    println("✅ REGION SPLITTING TEST PASSED")
    println("="^80)
    println()
    println("Region splitting is production-ready:")
    println("  • Smooth sigmoid transition function")
    println("  • Perfect population conservation")
    println("  • Correct decomposition: ψ = ψ_inner + ψ_outer")
    println("  • Ready for Volkov projection")
    println()
else
    println("❌ REGION SPLITTING TEST FAILED")
    println()
    println("Issues detected - review splitting implementation")
    exit(1)
end

println("="^80)
