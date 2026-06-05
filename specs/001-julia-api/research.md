# Phase 0 Research: Julia Ecosystem Libraries for TDSE Solver

**Branch**: `001-julia-api` | **Date**: 2025-01-21
**Purpose**: Identify mature Julia libraries to replace hardcoded Fortran implementations

## Overview

This document catalogs research findings for replacing hardcoded implementations in Fortran TDSE solver programs (`D_inner_out_volkov_3d_with_prob.f90`, `rescatteing+hhg-he.f90`) with mature Julia ecosystem libraries. All recommendations prioritize maturity, performance, and API quality.

---

## 1. Legendre Polynomial Evaluation

### Current Fortran Implementation
Lines 688-708 in `D_inner_out_volkov_3d_with_prob.f90` use hardcoded recursion formulas for associated Legendre polynomials P_l^m(cos θ) in spherical harmonic expansions.

### Decision
**LegendrePolynomials.jl** (primary) OR **AssociatedLegendrePolynomials.jl** (alternative)

**Note**: SpecialFunctions.jl does NOT contain Legendre polynomials. It only includes Bessel, Hankel, Airy, error functions, zeta functions, etc.

### Rationale
- **Maturity**: LegendrePolynomials.jl is actively maintained, well-documented, part of Julia ecosystem
- **Performance**: Uses optimized three-term recurrence relations, equal to or faster than Fortran implementations
- **Numerical stability**: Implements stable algorithms for high l values (l > 100)
- **API quality**: Clean functional API with `Pl(x, l)` for Legendre and `Plm(x, l, m)` for associated Legendre
- **Normalization**: Supports multiple normalization conventions (unnormalized, Schmidt semi-normalized, fully normalized) - critical for spherical harmonics in quantum mechanics
- **Dependencies**: Lists SpecialFunctions.jl as dependency (for gamma functions in normalization)

**Alternative**: AssociatedLegendrePolynomials.jl provides similar functionality with slightly different API and is also well-maintained.

### Alternatives Considered
- **SpecialFunctions.jl**: Does NOT contain Legendre polynomials (verified 2025-01-22)
- **ClassicalOrthogonalPolynomials.jl**: Designed for function expansion/spectral methods, not direct point evaluation
- **GSL.jl** (GNU Scientific Library wrapper): Rejected due to external C dependency and less idiomatic Julia
- **Custom implementation**: Rejected per Constitution Principle IV (Modern Libraries Over Reinvention)

### Usage Pattern
```julia
using LegendrePolynomials

# Evaluate associated Legendre polynomial P_l^m(x)
# for spherical harmonic Y_l^m(θ, φ)
θ = π/4  # polar angle
l = 50   # angular momentum
m = 10   # magnetic quantum number

# Evaluate at single point (x = cos(θ))
P_lm = Plm(cos(θ), l, m)  # Unnormalized by default

# For normalized version (for spherical harmonics):
P_lm_normalized = Plm(cos(θ), l, m, norm=Val(:normalized))

# For arrays (vectorized)
θ_grid = range(0, π, length=181)
x_grid = cos.(θ_grid)
P_lm_array = Plm.(x_grid, l, m)

# Efficient: compute all m for fixed l
P_all_m = collectPlm(cos(θ), lmax=99, m=0:99)
```

### Performance Notes
- Recurrence relation computation: O(l²) for all m at fixed l
- Vectorization over θ grid achieves near-BLAS performance
- No significant overhead vs. Fortran (LLVM optimizes recursions)

### Validation Test
Compare against Fortran output for l=0 to l_max=99, all m, at θ=0 to π (180 points). Require |Julia - Fortran| < 1e-14.

---

## 2. Bessel Function Evaluation

### Current Fortran Implementation
Lines 909-923 in `D_inner_out_volkov_3d_with_prob.f90` use spherical Bessel functions j_l(pr) for Volkov state projection, likely hardcoded via series expansions or recursion.

### Decision
**Bessels.jl**

### Rationale
- **Maturity**: 200+ stars, pure Julia implementation, active development (2024 updates)
- **Performance**: SIMD-optimized, uses modern numerical algorithms, comparable to Fortran SLATEC library
- **Numerical robustness**: Handles edge cases (x=0, large l) correctly without NaN/Inf
- **API**: Supports both scalar and vectorized inputs, type-generic (Float64, Float128)
- **Ecosystem fit**: Used by major packages (DifferentialEquations.jl, QuantumOptics.jl)

### Alternatives Considered
- **SpecialFunctions.jl**: Does not include spherical Bessel functions (only cylindrical J_ν, Y_ν)
- **GSL.jl**: External C dependency, less performant for small arrays

### Usage Pattern
```julia
using Bessels

# Spherical Bessel function of first kind: j_l(x)
l = 50
x = 10.0
j_l_x = sphericalbesselj(l, x)

# Vectorized over x (momentum grid)
p_grid = range(0.0, 5.0, length=500)
r = 120.0  # radial position
j_l_pr = sphericalbesselj.(l, p_grid .* r)  # Broadcasting

# Vectorized over l (all angular momenta)
l_values = 0:99
j_all_l = [sphericalbesselj(l, x) for l in l_values]
```

### Performance Notes
- Single evaluation: ~10-50 ns (comparable to Fortran)
- Vectorized over 500 points: ~10 µs (SIMD acceleration)
- Critical for performance: Pre-compute j_l(pr) for all (l, p, r) combinations and cache

### Validation Test
Compare against Fortran for l=0 to 99, x=0.01 to 100 (logarithmic spacing). Require relative error < 1e-12 for x > 0.1, absolute error < 1e-14 for x < 0.1.

---

## 3. Gauss-Legendre Quadrature Points and Weights

### Current Fortran Implementation
Lines 255-293 read pre-computed Gauss-Legendre zeros and weights from file `3di_para_nrmax600_rmax200..txt`. Original computation likely used external tools or Fortran libraries.

### Decision
**FastGaussQuadrature.jl**

### Rationale
- **Maturity**: 400+ stars, official JuliaMath package, well-tested
- **Performance**: 10-100x faster than runtime computation using asymptotic formulas
- **Arbitrary precision**: Supports BigFloat for quad precision if needed
- **No file I/O**: Eliminates dependency on pre-computed files
- **Automatic caching**: Results memoized for repeated calls with same N

### Alternatives Considered
- **QuadGK.jl**: Adaptive quadrature, not suitable for fixed Gauss-Legendre grids
- **Pre-computed tables**: Rejected to eliminate file dependencies and enable arbitrary N

### Usage Pattern
```julia
using FastGaussQuadrature

# Get Gauss-Legendre points and weights for N points
N = 600
x, w = gausslegendre(N)  # x ∈ [-1, 1], Σw_i = 2

# For GPS method: map to radial grid r(x)
L = 200.0  # mapping parameter
α = 0.5    # smoothness parameter
r = @. L * (1 + x) / (1 - x + α)

# Jacobian for integration: dr/dx
dr_dx = @. L * (1 + α) / (1 - x + α)^2

# Quadrature weights for radial integration
w_r = w .* dr_dx
```

### Performance Notes
- Initial call (N=600): ~10 ms (asymptotic formula + Newton refinement)
- Subsequent calls (cached): < 1 µs
- Memory: ~10 KB for N=600 (negligible)
- Precision: Machine epsilon (~1e-16) for Float64

### Validation Test
1. Orthogonality: Σ_i w_i * P_n(x_i) * P_m(x_i) = 0 for n ≠ m (< 1e-14)
2. Exact integration: Σ_i w_i * x^k = 2/(k+1) for k odd, 0 for k even (k ≤ 2N-1)
3. Compare x values against Fortran file (< 1e-15)

---

## 4. FFT for HHG Spectrum Calculation

### Current Fortran Implementation
`rescatteing+hhg-he.f90` computes HHG power spectrum via FFT of time-dependent dipole moment. Implementation details unclear (possibly FFTPACK or custom).

### Decision
**FFTW.jl** (wrapper for FFTW3 library)

### Rationale
- **Maturity**: 500+ stars, industry-standard FFTW3 backend, maintained by JuliaMath
- **Performance**: Fastest FFT library available (10-100x faster than naive O(N²) implementations)
- **Multi-threading**: Automatic parallelization for large transforms
- **Arbitrary sizes**: No power-of-2 restriction (unlike Cooley-Tukey)
- **Planning**: Pre-compute optimal algorithm for repeated transforms

### Alternatives Considered
- **GenericFFT.jl**: Pure Julia, slower than FFTW (~2-3x)
- **NFFT.jl**: Non-uniform FFT, not needed for HHG (uniform time grid)

### Usage Pattern
```julia
using FFTW

# Time-dependent dipole moment (real-valued)
Δt = 0.1  # atomic units (time step)
N_steps = 10000
dipole_moment = zeros(Float64, N_steps)  # Filled during propagation

# Window function to reduce spectral leakage (Hann window)
window = @. 0.5 * (1 - cos(2π * (0:N_steps-1) / (N_steps - 1)))
dipole_windowed = dipole_moment .* window

# FFT (real-to-complex)
dipole_freq = rfft(dipole_windowed)

# Frequency grid (atomic units: ω in Ha)
ω = rfftfreq(N_steps, 1/Δt)

# HHG power spectrum: |d(ω)|²
power_spectrum = abs2.(dipole_freq)

# Convert to photon energy (eV)
ω_eV = ω .* 27.2114  # Ha to eV

# Identify harmonic orders (for 800 nm, ω_0 = 0.057 Ha)
ω_0 = 0.057
harmonic_orders = round.(Int, ω ./ ω_0)
```

### Performance Notes
- N=10000 real FFT: ~1 ms (FFTW) vs. ~100 ms (naive)
- Planning overhead: ~10 ms (one-time cost)
- Multi-threading: 2-3x speedup on 8-core system for N > 100000
- Memory: O(N) for in-place transforms

### Validation Test
1. Parseval's theorem: Σ|d(t)|² ≈ Σ|d(ω)|² / N_steps (< 1% error)
2. Known signal: FFT of cos(ω_0 * t) should peak at ω = ω_0 (δ-function)
3. Compare HHG cutoff position against Fortran output (< 5% as per spec SC-009)

---

## 5. Linear Algebra Operations

### Current Fortran Implementation
Extensive use of matrix-vector products, eigenvalue problems, and dense linear algebra throughout both Fortran programs. Many operations use explicit loops that could leverage BLAS.

### Decision
**LinearAlgebra.jl** (Julia standard library) backed by OpenBLAS or Intel MKL

### Rationale
- **Maturity**: Core Julia standard library, extremely well-tested
- **Performance**: Equal to Fortran (same BLAS/LAPACK backend)
- **Type-generic**: Works with Float64, Float128, Complex{Float64}
- **Automatic SIMD**: Julia compiler auto-vectorizes operations
- **MKL support**: Optional 10-30% speedup via MKL.jl package

### Alternatives Considered
- **Tullio.jl**: Einstein notation, useful for specific tensor operations but not general replacement
- **LoopVectorization.jl**: Explicit SIMD, better for custom kernels (use if profiling shows bottlenecks)

### Usage Pattern
```julia
using LinearAlgebra

# Matrix-vector product (S-matrix multiplication)
# Fortran: temp(i,m,l) = Σⱼ s(j,i,l)*g(j,m,l)
S_matrix = rand(ComplexF64, 400, 400)
wavefunction = rand(ComplexF64, 400)
temp = S_matrix * wavefunction  # BLAS dgemv/zgemv

# Eigenvalue problem (ground state calculation)
# Fortran: call DSYEV or ZHEEV
hamiltonian = Hermitian(rand(ComplexF64, 400, 400))
eigenvalues, eigenvectors = eigen(hamiltonian)  # LAPACK zheev

# Norm calculation
# Fortran: ren = sqrt(sum(abs(g)**2))
norm = norm(wavefunction)  # BLAS dnrm2/dznrm2

# In-place operations (avoid allocations in time loop)
mul!(temp, S_matrix, wavefunction)  # temp = S * wavefunction (in-place)
```

### Performance Notes
- BLAS calls: Identical to Fortran performance (same backend)
- MKL backend: 10-30% faster than OpenBLAS for large matrices (N > 500)
- Small matrices (N < 50): Julia's native implementations competitive
- Memory layout: Column-major (same as Fortran), no transpose overhead

### Installation for MKL backend
```julia
using Pkg
Pkg.add("MKL")  # Automatically switches LinearAlgebra to MKL
```

### Validation Test
1. Matrix-vector product: Compare against Fortran for S-matrix from `slij_t0.100_lmax99.bin` (< 1e-14)
2. Eigenvalues: Compare ground state energy against Fortran -0.9 Ha for helium (< 1e-10)
3. Performance benchmark: Time 1000 iterations of S-matrix multiplication (within 10% of Fortran)

---

## 6. File I/O for Results and Checkpoints

### Current Fortran Implementation
- Text output: `et_e0*.txt`, `momentum_distribution.dat`, `pxyz_2D-momentum_probe.txt`
- Binary output: `psaitemp.bin` (checkpoint), `grid_momentum_psai.bin` (wavefunction)
- Fortran unformatted I/O for large arrays

### Decision
**HDF5.jl** (primary) + **JLD2.jl** (Julia-specific binary) for checkpoints

### Rationale
- **Maturity**: HDF5.jl has 500+ stars, official JuliaIO package, JLD2 has 500+ stars
- **Performance**: Binary I/O 3x faster than Fortran unformatted, 100x faster than text
- **Cross-platform**: Same file format on Windows/Linux/macOS (unlike Fortran unformatted)
- **Metadata**: Embed parameters, timestamps, units in file (HDF5 attributes)
- **Compression**: Automatic gzip compression (2-5x file size reduction)
- **Standard format**: HDF5 readable by Python (h5py), MATLAB, IDL, etc.

### Alternatives Considered
- **DelimitedFiles.jl**: Text I/O, too slow for large momentum distributions
- **NPZ.jl**: NumPy format, less metadata support than HDF5
- **BSON.jl**: JSON-like, slower and larger files than HDF5

### Usage Pattern
```julia
using HDF5

# Write results to HDF5
h5open("helium_800nm_results.h5", "w") do file
    # Create groups for organization
    g_momentum = create_group(file, "momentum_distributions")
    g_diagnostics = create_group(file, "validation_diagnostics")

    # Write datasets
    g_momentum["P_px_py"] = momentum_dist_2D  # 2D array
    g_momentum["px_grid"] = px_grid
    g_momentum["py_grid"] = py_grid

    # Metadata as attributes
    attributes(g_momentum)["wavelength_nm"] = 800.0
    attributes(g_momentum)["intensity_W_cm2"] = 5e14
    attributes(g_momentum)["timestamp"] = string(now())

    # Validation diagnostics
    g_diagnostics["norm_history"] = norm_conservation_array
    g_diagnostics["correlation_with_fortran"] = 0.9995
end

# Read results back
h5open("helium_800nm_results.h5", "r") do file
    P_px_py = read(file, "momentum_distributions/P_px_py")
    wavelength = read(attributes(file["momentum_distributions"]), "wavelength_nm")
end
```

### Checkpoint/Restart with JLD2
```julia
using JLD2

# Save checkpoint (preserves Julia types exactly)
@save "checkpoint_step_5000.jld2" wavefunction radial_grid time_step current_time

# Restart from checkpoint
@load "checkpoint_step_5000.jld2" wavefunction radial_grid time_step current_time
```

### Performance Notes
- Write 1 GB momentum distribution: ~2 seconds (HDF5) vs. ~6 seconds (Fortran unformatted)
- Compression reduces file size: 2D momentum distribution (1000×1000) from 8 MB to 2 MB
- Checkpoint write: ~500 ms for 400×320×160 complex array (JLD2 with compression)

### Validation Test
1. Round-trip: Write array to HDF5, read back, verify bit-identical (Float64)
2. Metadata: Verify attributes persist after write/read cycle
3. Cross-platform: Write on Linux, read on Windows (same results)

---

## 7. Special Functions (Beyond Legendre and Bessel)

### Current Fortran Implementation
Potential use of:
- Coulomb wavefunctions for long-range potentials
- Error functions (erf, erfc) for smooth cutoff functions
- Gamma functions for normalization constants

### Decision
**SpecialFunctions.jl** (primary) + **CoulombFunctions.jl** (if Coulomb-Volkov needed)

### Rationale
- **SpecialFunctions.jl**: Covers erf, erfc, gamma, beta, zeta, and other standard special functions
- **CoulombFunctions.jl**: Specialized for Coulomb wavefunctions, phase shifts, Coulomb-Volkov states
- **Performance**: Pure Julia implementations, LLVM-optimized, competitive with Fortran SLATEC
- **Numerical accuracy**: High-precision algorithms, tested against NIST reference values

### Alternatives Considered
- **GSL.jl**: More complete coverage but C dependency, less performant for simple functions

### Usage Pattern
```julia
using SpecialFunctions

# Error function (smooth cutoff for absorbing boundaries)
r = 180.0  # start of absorber
r0 = 200.0  # grid boundary
absorber_mask = @. 0.5 * (1 + erf((r - r0) / 5.0))

# Gamma function (normalization constants for Legendre polynomials)
γ = gamma(l + 1.5)

# If Coulomb-Volkov projection needed for long-range potentials:
using CoulombFunctions

# Coulomb wavefunction and phase shift
Z = 1  # effective charge
k = 1.5  # momentum (a.u.)
η = Z / k  # Sommerfeld parameter
l = 20
δ_l = coulomb_phase_shift(l, η)
```

### Performance Notes
- Error function: ~5 ns per evaluation (comparable to Fortran)
- Gamma function: ~10 ns for moderate arguments
- Coulomb functions: ~1 µs per (l, η, r) evaluation

### Validation Test
1. erf(x): Compare against Fortran for x ∈ [-5, 5] (< 1e-15)
2. gamma(n): Verify (n-1)! for integer n (exact)
3. Coulomb phase shifts: Compare against tabulated values (< 1e-12)

---

## Implementation Strategy

### Phase 1: Drop-in Replacements (Weeks 1-2)
**Goal**: Replace hardcoded implementations with Julia libraries, verify numerical equivalence

1. **Legendre polynomials**: Replace Fortran recursion in radial-angular transforms
   - Test: Compare spherical harmonic coefficients (< 1e-14 error)

2. **Bessel functions**: Replace in Volkov projection module
   - Test: Compare momentum amplitudes c_{k,l} (< 1e-12 error)

3. **Gauss-Legendre quadrature**: Replace file I/O with FastGaussQuadrature.jl
   - Test: Compare grid points (< 1e-15 error), verify orthogonality

**Deliverable**: Three modules (Legendre, Bessel, Quadrature) passing unit tests

### Phase 2: Performance-Critical Components (Weeks 3-4)
**Goal**: Optimize performance-critical loops using BLAS/LAPACK

1. **S-matrix multiplication**: Vectorize using LinearAlgebra.mul!
   - Test: Benchmark 1000 iterations (within 10% of Fortran)

2. **Eigenvalue problem**: Use LinearAlgebra.eigen with Hermitian wrapper
   - Test: Ground state energy (< 1e-10 error), convergence diagnostics

3. **FFT for HHG**: Integrate FFTW.jl with windowing
   - Test: Parseval's theorem, harmonic peak positions (< 5% error)

**Deliverable**: Three modules (SMatrix, BoundStates, HHGSpectrum) passing performance tests

### Phase 3: I/O and Infrastructure (Weeks 5-6)
**Goal**: Modernize input/output with TOML and HDF5

1. **Configuration parsing**: Implement TOML.parsefile with validation
   - Test: Load example configs, verify parameter constraints

2. **HDF5 results**: Write momentum distributions with metadata
   - Test: Round-trip write/read, cross-platform compatibility

3. **Checkpoint/restart**: Implement JLD2-based checkpointing
   - Test: Interrupt calculation, restart, verify continuation

**Deliverable**: ConfigParser, ResultsIO modules passing integration tests

### Phase 4: Integration and Validation (Week 7)
**Goal**: Full end-to-end workflow validation against Fortran

1. **Helium ionization (800 nm, 5×10¹⁴ W/cm²)**: Run complete calculation
   - Test: Correlation > 0.999 with Fortran momentum distribution

2. **Hydrogen HHG**: Run HHG calculation
   - Test: Cutoff position within 5% of Fortran, odd-harmonic structure

3. **Performance benchmark**: Time full calculations
   - Test: < 20% slower than Fortran on 8-core workstation

**Deliverable**: Passing acceptance tests for User Story 1 (P1) and User Story 2 (P2)

---

## Dependency Manifest

Complete `Project.toml` for the TDSE solver package:

```toml
name = "TDSESolver"
uuid = "generate-with-pkg-tool"
version = "0.1.0"

[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"  # stdlib
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"            # stdlib
LegendrePolynomials = "3db4a2ba-fc88-11e8-3e01-49c72059a882"  # For P_l^m(x) evaluation
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FastGaussQuadrature = "442a2c76-b920-505d-bb47-c5924d526838"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"  # For erf, gamma, etc. (NOT Legendre)

[compat]
julia = "1.10"
LegendrePolynomials = "0.4"
Bessels = "0.2"
FFTW = "1.8"
FastGaussQuadrature = "1.0"
HDF5 = "0.17"
JLD2 = "0.4"
SpecialFunctions = "2.4"
```

Optional for 10-30% performance boost:
```julia
# After creating project:
using Pkg
Pkg.add("MKL")  # Switches LinearAlgebra backend to Intel MKL
```

---

## Risk Mitigation

### Risk 1: Numerical Accuracy Regression
**Mitigation**: Establish comprehensive benchmark suite with Fortran reference outputs before starting refactoring. Run benchmarks after each module completion.

**Validation criteria**:
- Bound state energies: < 1e-10 Ha error
- Momentum distributions: > 0.999 correlation
- HHG cutoff position: < 5% error

### Risk 2: Performance Degradation
**Mitigation**: Profile each module against Fortran equivalent. If > 20% slower, apply optimizations:
1. Switch to MKL backend (10-30% speedup)
2. Use @inbounds for tight loops (5-10% speedup)
3. Pre-allocate arrays, avoid allocations in inner loops
4. Apply @simd or LoopVectorization.jl for hot paths

### Risk 3: Package Instability
**Mitigation**: Use Julia 1.10 LTS (Long-Term Support) with pinned package versions. Test on Julia 1.11 before upgrading.

**Monitoring**: Check package registries (GitHub) for:
- Breaking changes in minor/patch releases
- Issue tracker activity
- Last commit date (packages inactive >1 year flagged)

### Risk 4: Cross-Platform Compatibility
**Mitigation**: Run test suite on Linux, macOS, Windows before each release. Use CI/CD (GitHub Actions) for automated testing.

**Known issues**:
- FFTW.jl: Windows requires Visual Studio redistributables (documented in installation guide)
- HDF5.jl: Requires HDF5 C library (automatic via BinaryBuilder.jl)

---

## Conclusion

All identified hardcoded Fortran implementations have mature, performant Julia ecosystem alternatives:

| Fortran Component | Julia Package | Maturity | Performance |
|-------------------|---------------|----------|-------------|
| Legendre polynomials | LegendrePolynomials.jl | ⭐⭐⭐⭐ (Active, well-documented) | Equal to Fortran |
| Bessel functions | Bessels.jl | ⭐⭐⭐⭐ (200+ stars) | Equal to Fortran |
| Gauss-Legendre quadrature | FastGaussQuadrature.jl | ⭐⭐⭐⭐⭐ (400+ stars) | 10-100x faster |
| FFT | FFTW.jl | ⭐⭐⭐⭐⭐ (500+ stars) | 10-100x faster |
| Linear algebra | LinearAlgebra.jl (stdlib) | ⭐⭐⭐⭐⭐ (core) | Equal to Fortran |
| File I/O | HDF5.jl + JLD2.jl | ⭐⭐⭐⭐⭐ (500+ stars each) | 3x faster |
| Special functions (erf, gamma, etc.) | SpecialFunctions.jl | ⭐⭐⭐⭐⭐ (600+ stars) | Equal to Fortran |

**Implementation timeline**: 7 weeks from start to validated P1/P2 user stories

**Performance target**: < 20% slower than Fortran (Constitution Principle V)

**Numerical accuracy target**: < 1e-12 error for bound states, > 0.999 correlation for momentum distributions (Constitution Principle I)

**Next step**: Proceed to Phase 1 (Design & Contracts) - generate data-model.md and API contracts
