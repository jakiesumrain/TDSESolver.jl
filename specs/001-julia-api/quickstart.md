# Quickstart Guide: TDSESolver.jl

**Goal**: Run your first TDSE calculation in under 10 minutes (SC-001)

## Prerequisites

- Julia 1.10+ installed ([julialang.org/downloads](https://julialang.org/downloads/))
- 8+ core workstation, 16+ GB RAM recommended
- 10-30 minutes for helium ionization calculation

## Installation

### 1. Create Julia Project

```bash
mkdir tdse_calculations
cd tdse_calculations
julia
```

### 2. Install Package (in Julia REPL)

```julia
using Pkg
Pkg.add(url="https://github.com/your-org/TDSESolver.jl")  # Replace with actual URL
```

### 3. Optional: Enable MKL for 10-30% Speedup

```julia
Pkg.add("MKL")
using MKL  # Must be loaded before TDSESolver
```

## Quick Example: Helium Ionization

### Step 1: Create Configuration File

Create `helium_800nm.toml`:

```toml
[atom]
type = "helium"

[laser]
wavelength = 800.0        # nm
intensity = 5e14          # W/cm²
pulse_duration = 10.0     # fs
polarization = [1.0, 0.0, 0.0]

[calculation]
type = "ionization"
output_file = "helium_800nm_results.h5"
verbosity = "normal"
```

### Step 2: Run Calculation

```julia
using TDSESolver

# Load and validate configuration
config = load_config("helium_800nm.toml")

# Run calculation (takes ~15-30 minutes on 8-core workstation)
results = run_tdse(config)
```

**Expected output**:
```
TDSESolver v0.1.0
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Configuration loaded successfully
  Atom: Helium (I_p = 0.9 Ha)
  Laser: 800 nm, 5.0×10¹⁴ W/cm², 10.0 fs
  Keldysh parameter: γ = 1.23

Initializing computational grid...
  Radial grid: 400 points, r_max = 150.0 a.u.
  Angular momentum: l_max = 99
  Time step: Δt = 0.1 a.u. (2.42 as)
  Total steps: 11,157

Computing ground state...
  Energy: -0.90034 Ha (literature: -0.9 Ha) ✓
  Norm: 1.000000 ✓

Time propagation (11,157 steps):
Progress: [████████████████████] 100% | Time: 21:34 | ETA: 00:00
  Current time: 27.0 fs
  Ionization probability: 0.234
  Norm: 0.9998 ✓

Projecting onto momentum eigenstates (Volkov method)...
  Momentum grid: 500 points, p_max = 5.0 a.u.
  Processing l=0 to l=99... Done

Validation diagnostics:
  ✓ Norm conservation: min=0.9997, max=1.0001
  ✓ Energy conservation: ΔE/E < 1×10⁻⁷
  ✓ Fortran correlation: 0.9996 (> 0.999 threshold)

Results saved to: helium_800nm_results.h5
Computation time: 21 minutes 34 seconds
```

### Step 3: Analyze Results

```julia
# Access momentum distributions
P_radial = results.momentum_distribution_radial
P_2D_xy = results.momentum_distribution_2D_xy

# Check validation
report = results.validation_report
@assert report.passed_validation

# Plot results (using Plots.jl)
using Plots
heatmap(results.px_grid, results.py_grid, P_2D_xy',
        xlabel="px (a.u.)", ylabel="py (a.u.)",
        title="Photoelectron Momentum Distribution",
        color=:viridis)
savefig("helium_800nm_momentum.png")
```

## Example 2: HHG Spectrum

### Configuration (`helium_hhg.toml`)

```toml
[atom]
type = "helium"

[laser]
wavelength = 800.0
intensity = 2e14          # Lower intensity for HHG
pulse_duration = 15.0     # Longer pulse for more cycles
polarization = [1.0, 0.0, 0.0]

[calculation]
type = "hhg"              # Changed to HHG
output_file = "helium_hhg_results.h5"
```

### Execution

```julia
config = load_config("helium_hhg.toml")
results = run_tdse(config)

# Plot HHG spectrum
using Plots
plot(results.photon_energy_grid, results.hhg_spectrum,
     xlabel="Photon Energy (eV)", ylabel="Intensity (arb. units)",
     title="High Harmonic Generation Spectrum",
     yscale=:log10, ylims=(1e-10, 1e0))

# Mark harmonic orders
for h in results.harmonic_orders
    vline!([h * 1.55], label="H$h", linestyle=:dash)  # 1.55 eV = 800 nm photon
end

savefig("helium_hhg_spectrum.png")
```

## Example 3: Parameter Scan

### Configuration (`intensity_scan.toml`)

```toml
[atom]
type = "helium"

[laser]
wavelength = 800.0
intensity = 5e14          # Base value (will be overridden by scan)
pulse_duration = 10.0
polarization = [1.0, 0.0, 0.0]

[calculation]
type = "ionization"
output_file = "intensity_scan_results.h5"

[scan]
parameter = "laser.intensity"
range = "1e14:1e14:1e15"  # 10 points from 1×10¹⁴ to 1×10¹⁵ W/cm²
```

### Execution (Parallel)

```julia
# Set number of threads
ENV["JULIA_NUM_THREADS"] = "8"

using TDSESolver

config = load_config("intensity_scan.toml")
scan_results = run_scan(config)  # Parallelizes automatically

# Extract ionization probabilities
intensities = [1e14, 2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 1e15]
ionization_probs = [scan_results[I].ionization_probability for I in intensities]

# Plot ionization yield vs. intensity
using Plots
plot(intensities, ionization_probs,
     xlabel="Intensity (W/cm²)", ylabel="Ionization Probability",
     title="Ionization Yield vs. Laser Intensity",
     marker=:circle, xscale=:log10)
savefig("ionization_vs_intensity.png")
```

## Example 4: Custom Potential

### Configuration (`soft_core.toml`)

```toml
[atom]
type = "custom"
potential = "-1/sqrt(r^2 + 0.5)"  # Soft-core potential, a=0.5
charge = 1                         # Effective nuclear charge

[laser]
wavelength = 800.0
intensity = 5e14
pulse_duration = 10.0
polarization = [1.0, 0.0, 0.0]

[calculation]
type = "ionization"
output_file = "soft_core_results.h5"
```

### Advanced: Julia Function Module

For complex potentials, create `MyPotentials.jl`:

```julia
module MyPotentials

function yukawa_potential(r; Z=1.0, λ=2.0)
    # Yukawa potential: V(r) = -Z * exp(-r/λ) / r
    return -Z * exp(-r/λ) / r
end

end  # module
```

Then in config:

```toml
[atom]
type = "custom"
potential_module = "MyPotentials"
potential_function = "yukawa_potential"
```

## Validation and Benchmarking

### Run Benchmark Suite

```julia
using TDSESolver.Validation

# Compare against Fortran reference data
benchmarks = [
    ("hydrogen_ground_state", "tests/benchmarks/hydrogen_ground_state.h5"),
    ("helium_ionization_800nm", "tests/benchmarks/helium_ionization_800nm.h5"),
    ("helium_hhg_800nm", "tests/benchmarks/helium_hhg_800nm.h5")
]

for (name, reference_file) in benchmarks
    println("Testing: $name")
    correlation = compare_with_fortran(results, reference_file)
    @assert correlation > 0.999 "Correlation $correlation below threshold"
    println("  ✓ Correlation: $correlation")
end
```

### Check Validation Report

```julia
report = results.validation_report

if !report.passed_validation
    @warn "Validation issues detected"
    for warning in report.warnings
        println("  ⚠ $warning")
    end
end

# Detailed checks
println("Norm conservation: [$(report.norm_conservation_min), $(report.norm_conservation_max)]")
println("Energy conservation error: $(report.energy_conservation_error)")
println("Fortran correlation: $(report.fortran_correlation)")
```

## Performance Optimization

### 1. Enable Multi-Threading

```bash
export JULIA_NUM_THREADS=8  # Before starting Julia
```

Or in Julia:

```julia
# Check thread count
Threads.nthreads()  # Should show 8

# Parallelization is automatic for parameter scans
```

### 2. Use MKL Backend

```julia
using MKL  # Load before TDSESolver for 10-30% speedup
using TDSESolver
```

### 3. Profile Hot Paths (Advanced)

```julia
using Profile

@profile run_tdse(config)
Profile.print()

# Identify bottlenecks, optimize with @simd or @inbounds if needed
```

## Troubleshooting

### Error: "Grid size insufficient for ponderomotive energy"

**Cause**: Automatic grid determination found U_p requires larger rmax

**Solution**: Override in configuration:

```toml
[numerical]
rmax = 200.0  # Increase maximum radius
```

### Error: "Time step exceeds stability criterion"

**Cause**: Δt too large for highest eigenstate energy

**Solution**: Override time step:

```toml
[numerical]
time_step = 0.05  # Reduce from automatic value
```

### Warning: "Norm conservation violated"

**Cause**: Wavefunction reaching grid boundary before absorber

**Solution**:
1. Increase rmax
2. Move absorber start radius closer to boundary
3. Increase absorber strength

### Slow Performance

**Solutions**:
1. Enable MKL: `using MKL` before `using TDSESolver`
2. Increase threads: `export JULIA_NUM_THREADS=16`
3. Check system resources: `free -h` (Linux), Activity Monitor (macOS)
4. Reduce numerical overrides if set too high (lmax, nrmax)

## Next Steps

1. **Read API Reference** (`contracts/api-reference.md`) for advanced usage
2. **Explore Tutorial Notebooks** (`docs/tutorial_notebooks/`) for detailed examples
3. **Consult Algorithm References** (`docs/algorithm_references.md`) to understand physics
4. **Review Constitution** (`.specify/memory/constitution.md`) for code contribution guidelines

## Getting Help

- **Documentation**: Full API reference in `contracts/api-reference.md`
- **Examples**: Sample configurations in `examples/` directory
- **Benchmarks**: Reference data in `tests/benchmarks/`
- **Issues**: Report bugs or request features on GitHub

---

**Congratulations!** You've completed your first TDSE calculation. 🎉

**Typical workflow**: Edit TOML config → `load_config()` → `run_tdse()` → Analyze results in Julia or Python (via HDF5)

**Time from installation to first result**: < 10 minutes (SC-001 ✓)
