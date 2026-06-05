# Analysis Modes: Ionization, Momentum Distributions, and HHG

This document shows how to configure simulations for different physics analysis using independent control switches.

## Control Switches

Two independent Boolean flags control analysis capabilities:

```julia
enable_momentum_analysis::Bool = false  # Photoelectron momentum distributions
enable_hhg::Bool = false                # High harmonic generation spectrum
```

Both can be enabled/disabled independently, giving **4 operational modes**.

---

## Mode 1: Ionization Only (Fastest)

**Use case:** Basic ionization probability, no spectral analysis

```julia
params = create_default_params(
    # Standard parameters
    atom = :hydrogen,
    laser_intensity = 1.0e14,
    laser_wavelength = 800.0,
    t_total = 100.0,

    # Analysis switches (both OFF)
    enable_momentum_analysis = false,  # No momentum distributions
    enable_hhg = false                 # No HHG spectrum
)

results = run_simulation(params)

# Available results:
# - results.observables.ionization_probs
# - results.observables.norms
# - results.observables.energies
```

**Output:**
- Ionization probability P_ion(t)
- Norm conservation checks
- Energy expectation values

**Typical runtime:** ~10 seconds (small grid)

---

## Mode 2: Ionization + Momentum Distributions

**Use case:** Photoelectron momentum analysis, ATI spectra

```julia
params = create_default_params(
    # Standard parameters
    atom = :hydrogen,
    laser_intensity = 5.0e14,  # Higher for ionization
    laser_wavelength = 800.0,
    t_total = 100.0,

    # Momentum analysis ON
    enable_momentum_analysis = true,
    msplit = 50,                # Split every 50 steps
    R_c = 100.0,                # Splitting radius
    delta_split = 5.0,          # Smoothness
    p_max = 2.0,                # Max momentum
    n_p = 30,                   # Momentum grid resolution
    n_theta = 20,
    n_phi = 20,

    # HHG OFF
    enable_hhg = false
)

results = run_simulation(params)

# Available results:
# - results.observables.ionization_probs
# - results.momentum_grid              # Momentum space grid
# - results.psai_p                     # ψ(φ,θ,p) momentum wavefunction

# Post-process momentum distributions:
using .MomentumDistribution

rate_xyz = transform_to_cartesian_grid(results.psai_p,
                                       results.momentum_grid,
                                       params.p_max,
                                       30)  # Cartesian grid size

pxy_rate, pxz_rate, pyz_rate = compute_integrated_2d_distributions(
    rate_xyz, params.p_max, 30
)
```

**Output:**
- All Mode 1 outputs
- **Momentum distributions:** P(px,py,pz), P(px,py), P(px,pz), P(py,pz)
- **ATI spectra:** Photoelectron energy and angular distributions

**Typical runtime:** ~20 seconds (more projections)

---

## Mode 3: Ionization + HHG Spectrum

**Use case:** High harmonic generation, attosecond physics

```julia
params = create_default_params(
    # Standard parameters
    atom = :hydrogen,
    laser_intensity = 1.0e14,
    laser_wavelength = 800.0,
    t_total = 100.0,

    # Momentum analysis OFF
    enable_momentum_analysis = false,

    # HHG ON
    enable_hhg = true,
    hhg_method = :acceleration,     # :eigenstate, :acceleration, :length
    hhg_record_interval = 1         # Every step for accurate spectrum
)

results = run_simulation(params)

# Available results:
# - results.observables.ionization_probs
# - results.observables.dipole_moment   # d(t) time series

# Post-process HHG spectrum:
using .HHG
using .PhysicalUnits: wavelength_nm_to_frequency_au

ω₀ = wavelength_nm_to_frequency_au(params.laser_wavelength)
h_orders, power, freqs = compute_hhg_spectrum(
    results.observables.dipole_moment,
    params.dt,
    ω₀,
    window=:hann
)

# Save spectrum
save_hhg_spectrum_to_file(h_orders[1:div(end,2)],
                          power[1:div(end,2)],
                          "hhg_spectrum.txt")
```

**Output:**
- All Mode 1 outputs
- **Dipole moment:** d(t) time series (complex)
- **HHG spectrum:** S(ω) = |d(ω)|², harmonic orders

**Typical runtime:** ~12 seconds (dipole calculation adds minimal overhead)

---

## Mode 4: Full Analysis (Ionization + Momentum + HHG)

**Use case:** Complete strong-field ionization study

```julia
params = create_default_params(
    # Standard parameters
    atom = :hydrogen,
    laser_intensity = 5.0e14,
    laser_wavelength = 800.0,
    t_total = 100.0,

    # BOTH ON
    enable_momentum_analysis = true,
    msplit = 50,
    R_c = 100.0,
    delta_split = 5.0,
    p_max = 2.0,
    n_p = 30,
    n_theta = 20,
    n_phi = 20,

    enable_hhg = true,
    hhg_method = :acceleration,
    hhg_record_interval = 1
)

results = run_simulation(params)

# Available results:
# - Ionization probability
# - Momentum distributions (photoelectron spectra)
# - HHG spectrum (attosecond emission)
```

**Output:**
- All Mode 1, 2, and 3 outputs combined
- Complete strong-field ionization analysis

**Typical runtime:** ~25 seconds (both analyses enabled)

---

## Quick Reference Table

| Mode | `enable_momentum_analysis` | `enable_hhg` | Use Case | Runtime |
|------|---------------------------|--------------|----------|---------|
| **1** | `false` | `false` | Ionization only | Fastest |
| **2** | `true`  | `false` | Ionization + photoelectron spectra | Medium |
| **3** | `false` | `true`  | Ionization + HHG spectrum | Fast |
| **4** | `true`  | `true`  | Full analysis (all physics) | Slowest |

---

## HHG Method Selection

When `enable_hhg = true`, choose dipole calculation method:

```julia
hhg_method = :eigenstate     # Method 1 (Fortran active, eigenstate expansion)
hhg_method = :acceleration   # Method 2a (Fortran commented, acceleration form)
hhg_method = :length         # Method 2b (length gauge)
```

**Recommendation:**
- **For production:** Use `:acceleration` or `:length` (immune to eigenstate errors)
- **For testing:** Use `:eigenstate` (once GPS eigenstate solver is fixed)

---

## Performance Considerations

### Memory Usage

- **Mode 1 (Ionization):** ~100 MB
- **Mode 2 (+ Momentum):** ~500 MB (stores ψ(φ,θ,p))
- **Mode 3 (+ HHG):** ~150 MB (stores d(t) time series)
- **Mode 4 (Both):** ~550 MB

### Computational Cost

**Dominant costs:**
1. **Time propagation:** O(nrmax² × lmax × n_steps) - Always required
2. **Region splitting:** O(nrmax × lmax × n_steps/msplit) - Only if `enable_momentum_analysis=true`
3. **Volkov projection:** O(nrmax × n_p × lmax) - Only if `enable_momentum_analysis=true`
4. **Dipole calculation:** O(nrmax × lmax) - Only if `enable_hhg=true`, negligible cost

**Rule of thumb:**
- Momentum analysis adds ~50% overhead
- HHG adds <5% overhead
- Both together adds ~55% overhead

---

## Example: Switching Between Modes

```julia
# Define base parameters once
base_params = (
    atom = :hydrogen,
    nrmax = 200,
    rmax = 100.0,
    lmax = 2,
    dt = 0.1,
    t_total = 100.0,
    laser_wavelength = 800.0,
    laser_intensity = 5.0e14,
    laser_duration = 50.0
)

# Run Mode 1: Ionization only (quick check)
params1 = create_default_params(;
    base_params...,
    enable_momentum_analysis = false,
    enable_hhg = false
)
results1 = run_simulation(params1)

# Run Mode 2: Add momentum analysis
params2 = create_default_params(;
    base_params...,
    enable_momentum_analysis = true,
    enable_hhg = false
)
results2 = run_simulation(params2)

# Run Mode 3: Add HHG instead
params3 = create_default_params(;
    base_params...,
    enable_momentum_analysis = false,
    enable_hhg = true,
    hhg_method = :acceleration
)
results3 = run_simulation(params3)

# Run Mode 4: Full analysis
params4 = create_default_params(;
    base_params...,
    enable_momentum_analysis = true,
    enable_hhg = true,
    hhg_method = :acceleration
)
results4 = run_simulation(params4)
```

---

## Typical Workflow

### Research Workflow

1. **Exploratory run (Mode 1):**
   - Quick parameter sweep
   - Check ionization yields
   - Verify field strength

2. **Detailed photoelectron analysis (Mode 2):**
   - Compute momentum distributions
   - Extract ATI peaks
   - Analyze angular distributions

3. **HHG spectrum (Mode 3):**
   - Compute harmonic spectrum
   - Identify cutoff energies
   - Analyze phase matching

4. **Publication-quality (Mode 4):**
   - Run with both analyses
   - Generate all figures
   - Cross-validate results

### Development Workflow

1. **Fast testing (Mode 1):** Verify code changes quickly
2. **Feature testing (Mode 2 or 3):** Test specific analysis
3. **Integration testing (Mode 4):** Verify complete pipeline

---

## Summary

✅ **Two independent switches:** `enable_momentum_analysis` and `enable_hhg`

✅ **Four operational modes:** Combine switches for different physics

✅ **Runtime control:** No code changes needed, just parameter flags

✅ **Flexible:** Enable only the analysis you need

✅ **Efficient:** Only pay for what you use

The system is designed for maximum flexibility - you control exactly what physics you want to calculate!
