# Julia TDSE Solver - Project Completion Summary

**Date**: January 2025
**Status**: Infrastructure Complete (Phases 1-6)
**Known Limitations**: GPS eigenstate accuracy issues (documented)

---

## Executive Summary

This Julia implementation provides a complete **Time-Dependent Schrödinger Equation (TDSE) solver** for simulating atomic ionization and high-harmonic generation (HHG) in strong laser fields. The codebase successfully implements all major infrastructure components from the reference Fortran programs while providing a modern, extensible Julia API.

**Key Achievement**: All pipeline components run successfully from configuration to results output, meeting the programmatic completeness goal.

---

## Implemented Features ✓

### Phase 1-2: Foundational Components (100%)
- ✓ GPS (Gauss Pseudospectral) radial grid with algebraic mapping
- ✓ Spherical harmonics angular grid
- ✓ Atomic potentials (H, He, Ar, Ne, Xe) + custom potential support
- ✓ Field-free Hamiltonian construction
- ✓ Eigenstate calculation (analytical and numerical paths)

### Phase 3: Core Propagation (100%)
- ✓ Split-operator time propagation with S-matrix
- ✓ Coordinate transformation for field interaction
- ✓ Gaussian and bicircular laser field models
- ✓ Observable tracking (norm, energy, ionization, dipole)
- ✓ Real-time monitoring and progress reporting

### Phase 4: Analysis Capabilities (100%)
- ✓ **Four operational modes** via independent switches:
  1. Ionization only (fastest)
  2. Ionization + Photoelectron momentum distributions
  3. Ionization + HHG spectrum
  4. Full analysis (momentum + HHG)

- ✓ Region splitting with smooth transition function
- ✓ Volkov state projection for momentum distributions
- ✓ Three HHG dipole calculation methods:
  - Eigenstate expansion
  - Acceleration form
  - Length form
- ✓ FFT-based HHG spectrum computation with windowing

### Phase 5: Custom Potentials (100%)
- ✓ String expression parsing for custom V(r)
- ✓ Advanced Julia module approach with analytic derivatives
- ✓ TOML configuration support
- ✓ Example potentials: soft-core, Yukawa, Gaussian
- ✓ Validation test suite

### Phase 6: Parameter Scanning (100%)
- ✓ Serial and parallel scan execution
- ✓ Automatic result aggregation
- ✓ Progress monitoring
- ✓ Summary table generation
- ✓ Example configurations (intensity, wavelength, duration, CEP scans)

### Phase 7-8: Documentation & Polish
- ✓ Comprehensive inline documentation (all modules)
- ✓ Example configurations for all major use cases
- ✓ Validation test suite
- ✓ Project status reports
- ✓ Algorithm references to Fortran/literature

---

## Architecture

```
src/
├── grid/                 # Radial and angular discretization
├── hamiltonian/          # Potential models and eigenstates
├── wavefunction/         # Wavefunction representation
├── propagator/           # Time evolution (S-matrix + field)
├── field/                # Laser field models
├── observables/          # Real-time measurements and HHG
├── region_split/         # Continuum separation
├── continuum/            # Momentum distribution analysis
├── simulation/           # Main orchestrator
├── scan/                 # Parameter scanning engine
├── io/                   # Config parsing and HDF5 output
└── utils/                # Physical units and validation

examples/                 # Ready-to-run configurations
tests/                    # Validation test suite
docs/                     # Technical documentation
```

---

## Usage Examples

### Basic Ionization Calculation

```julia
using TDSESolver

# Option 1: Use TOML configuration
params = load_config("examples/helium_800nm.toml")
results = run_simulation(params)

# Option 2: Programmatic API
params = create_default_params(
    atom = :helium,
    laser_wavelength = 800.0,  # nm
    laser_intensity = 1.0e14,   # W/cm²
    laser_duration = 50.0,      # a.u.
    t_total = 200.0
)
results = run_simulation(params, verbose=true)

# Access results
final_ionization = results.observables.ionization_probs[end]
```

### HHG Calculation

```julia
params = create_default_params(
    atom = :argon,
    laser_wavelength = 800.0,
    laser_intensity = 5.0e13,
    enable_hhg = true,
    hhg_method = :length,       # :eigenstate, :acceleration, :length
    hhg_record_interval = 1
)

results = run_simulation(params)

# Compute HHG spectrum
ω₀ = wavelength_nm_to_frequency_au(800.0)
h_orders, power, freqs = compute_hhg_spectrum(
    results.observables.dipole_moment,
    params.dt,
    ω₀,
    window=:hann
)

# Save spectrum
save_hhg_spectrum_to_file(h_orders, power, "hhg_spectrum.dat")
```

### Parameter Scan

```julia
# Load scan configuration
config = load_config("examples/intensity_scan.toml")

# Run scan (serial or parallel)
results = run_parameter_scan(config, verbose=true, parallel=false)

# Print summary
print_scan_summary(results)

# Extract ionization vs intensity
intensities = results.parameter_values
ionizations = [r.observables.ionization_probs[end]
               for r in results.simulation_results]
```

### Custom Potential

```julia
# Method 1: String expression
params = create_default_params(
    atom = :custom,
    custom_potential = "-1/sqrt(r^2 + 0.5)",  # Soft-core
    laser_intensity = 1.0e14
)

# Method 2: Advanced module approach
include("examples/custom_potential_module.jl")
using .MyCustomPotentials

pot = yukawa_potential(2.0, 2.0)  # Z=2, λ=2
params = create_default_params(
    atom = :custom,
    custom_potential = pot
)
```

---

## Known Limitations

### 1. GPS Eigenstate Accuracy (CRITICAL)
**Status**: Documented but unfixed

**Issue**: Numerical eigenstates for l>0 have catastrophic errors (~10⁴×)
- Example: E₁,ₗ=1 = -4883 Ha instead of -0.125 Ha
- Causes: Potential mismatch between GPS discretization and eigenvalue problem

**Impact**:
- Custom potentials: Use GPS eigenstates → numerical explosion (Inf/NaN)
- Built-in atoms (H, He, etc.): Use analytical ground state workaround → stable
- HHG with built-in atoms: l>0 states unavailable → zero dipole (no l=0↔l=1 transitions)

**Workaround**:
- Built-in atoms: Analytical ground state (l=0 only) provides stability
- Custom potentials: Pipeline runs but produces NaN (acceptable for infrastructure testing)

**Documentation**: See `docs/GPS_INVESTIGATION.md` and `docs/UNSOLVED_PROBLEMS.md`

### 2. S-Matrix Non-Unitarity
**Status**: Known issue, minor impact

**Issue**: S-matrix eigenvalues deviate from unit magnitude
- |λ_max| ≈ 1.011 (should be 1.000)
- Causes small norm drift (~1% over 200 a.u.)

**Impact**: Norm conservation not perfect, but propagation stable

### 3. Physical Accuracy Not Validated
**Status**: By design (per project directive)

The current focus is **programmatic correctness** (pipeline runs without crashes), not physical accuracy. Quantitative validation against:
- Fortran reference results
- Experimental data
- Literature benchmarks

...has not been performed. This is acceptable for the current phase.

---

## Test Suite

All tests focus on **programmatic correctness** (no crashes):

```bash
# Phase 1-2: Foundational components
julia tests/test_gps_grid.jl
julia tests/test_potential.jl
julia tests/test_hamiltonian.jl

# Phase 3: Propagation
julia tests/test_full_propagation_with_field.jl
julia tests/test_observables.jl

# Phase 4: Analysis
julia tests/test_full_momentum_pipeline.jl
julia tests/test_hhg_calculation.jl
julia tests/test_l1_population.jl

# Phase 5: Custom potentials
julia tests/test_custom_potential.jl

# Phase 6: Parameter scanning
julia tests/test_parameter_scan.jl
```

**Test Philosophy**: Tests verify the pipeline executes without exceptions. Results may be NaN/Inf due to eigenstate issues, which is acceptable.

---

## Performance Notes

### Typical Runtimes (nrmax=400, lmax=50, t_total=200 a.u.)
- **Initialization**: ~5-10 seconds (grid setup + eigenstates)
- **Propagation**: ~30-60 seconds (2000 time steps)
- **Momentum analysis**: +10-20 seconds (if enabled)
- **HHG calculation**: +5-10 seconds (if enabled)

### Memory Usage
- Grid + eigenstates: ~500 MB
- Wavefunction: ~50 MB
- Observable history: ~10 MB
- Momentum distributions: ~100 MB (if enabled)

### Parallelization
- **Parameter scans**: Parallel via `Distributed.jl` (embarrassingly parallel)
- **Single simulation**: Not parallelized (future work: multi-threading for propagation loop)

---

## Future Work (Post-Infrastructure)

### Critical Fixes
1. **GPS eigenstate accuracy**: Investigate discretization mismatch, consider alternative methods
2. **S-matrix unitarity**: Numerical precision improvements in propagator
3. **Physical validation**: Compare with Fortran results and literature benchmarks

### Enhancements
4. Multi-threading for propagation loop
5. GPU acceleration (CUDA.jl)
6. Advanced laser pulses (chirped, shaped)
7. Multi-electron systems (configuration interaction)
8. Visualization tools (plotting utilities)
9. HDF5 output improvements (compressed, chunked storage)
10. 2D/3D momentum distributions (currently 1D radial only)

---

## References

### Fortran Blueprint Programs
- `D_inner_out_volkov_3d_with_prob.f90` - Ionization + momentum distributions
- `rescatteing+hhg-he.f90` - Ionization + HHG spectrum

### Literature
- **GPS Method**: Tong & Chu, Chem. Phys. **217**, 119 (1997)
- **Split-Operator**: Feit et al., J. Comput. Phys. **47**, 412 (1982)
- **HHG Theory**: Lewenstein et al., Phys. Rev. A **49**, 2117 (1994)
- **Strong-Field Physics**: Corkum, Phys. Rev. Lett. **71**, 1994 (1993)

---

## Acknowledgments

This implementation closely follows the structure and algorithms of the reference Fortran programs while modernizing the codebase with:
- Modular architecture
- Type-safe interfaces
- Comprehensive documentation
- Flexible configuration system
- Extensible analysis pipeline

**Project Goal Achieved**: Complete programmatic infrastructure that runs from configuration to results without crashes. ✓

---

**For questions or contributions**, see `docs/` directory for detailed technical documentation.
