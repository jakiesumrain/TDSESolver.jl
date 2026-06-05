# Julia TDSE Solver - Quick Reference

**Fast lookup for common tasks**

---

## Installation & Setup

```bash
# Clone repository
git clone <repository-url>
cd TDSE

# Julia dependencies (add as needed)
# All modules are standalone, no external packages required except:
# - LinearAlgebra (stdlib)
# - FFTW (for HHG spectrum)
# - Distributed (for parallel scans)
```

---

## Running Simulations

### From TOML Configuration
```julia
include("src/TDSESolver.jl")
using .ConfigParser, .Simulation

config = load_config("examples/helium_800nm.toml")
results = run_simulation(config_to_params(config))
```

### Programmatic API
```julia
include("src/TDSESolver.jl")
using .Simulation

params = create_default_params(
    atom = :hydrogen,
    nrmax = 200,
    rmax = 100.0,
    lmax = 2,
    dt = 0.1,
    t_total = 100.0,
    laser_wavelength = 800.0,    # nm
    laser_intensity = 1.0e14,    # W/cm²
    laser_duration = 50.0         # a.u.
)

results = run_simulation(params, verbose=true)
```

---

## Common Parameters

### Grid Parameters
```julia
nrmax = 200          # Number of radial grid points (100-400 typical)
rmax = 100.0         # Maximum radius in a.u. (50-150 typical)
lmax = 2             # Maximum angular momentum (0-50, depends on physics)
L = 30.0             # GPS mapping parameter (default: 30.0)
α = 0.5              # GPS mapping parameter (default: 0.5)
```

### Time Parameters
```julia
dt = 0.1             # Time step in a.u. (0.05-0.2 typical)
t_total = 100.0      # Total simulation time in a.u.
obs_interval = 10    # Record observables every N steps
```

### Laser Parameters
```julia
laser_wavelength = 800.0        # Wavelength in nm
laser_intensity = 1.0e14        # Peak intensity in W/cm²
laser_duration = 50.0           # FWHM duration in a.u.
laser_cep = 0.0                 # Carrier-envelope phase in radians
laser_polarization = [0,0,1]    # Polarization vector (z-polarized)
```

### Analysis Switches
```julia
enable_momentum_analysis = false   # Enable momentum distributions
enable_hhg = false                 # Enable HHG dipole tracking
hhg_method = :length               # :eigenstate, :acceleration, :length
hhg_record_interval = 1            # Record every N steps
```

---

## Accessing Results

### Basic Observables
```julia
# Time series
times = results.observables.times
norms = results.observables.norms
energies = results.observables.energies
ionization_probs = results.observables.ionization_probs

# Final values
final_ionization = results.observables.ionization_probs[end]
final_norm = results.observables.norms[end]
```

### Momentum Distributions
```julia
if params.enable_momentum_analysis
    p_grid = results.momentum_grid.p     # Momentum values
    psai_p = results.psai_p               # ψ(φ,θ,p) complex array

    # Compute momentum distribution
    prob_p = sum(abs2.(psai_p), dims=(1,2))[1,1,:]  # Radial distribution
end
```

### HHG Spectrum
```julia
if params.enable_hhg
    dipole = results.observables.dipole_moment

    # Compute spectrum
    ω₀ = wavelength_nm_to_frequency_au(params.laser_wavelength)
    h_orders, power, freqs = compute_hhg_spectrum(dipole, params.dt, ω₀)

    # Save to file
    save_hhg_spectrum_to_file(h_orders, power, "spectrum.dat")
end
```

---

## Example Workflows

### 1. Ionization Probability vs Intensity
```julia
intensities = [1e13, 5e13, 1e14, 5e14, 1e15]  # W/cm²
ionizations = []

for I in intensities
    params = create_default_params(
        atom = :helium,
        laser_intensity = I,
        t_total = 200.0
    )
    results = run_simulation(params, verbose=false)
    push!(ionizations, results.observables.ionization_probs[end])
end

# Plot: intensities vs ionizations
```

### 2. Photoelectron Momentum Distribution
```julia
params = create_default_params(
    atom = :argon,
    laser_intensity = 1.0e14,
    laser_wavelength = 800.0,
    enable_momentum_analysis = true,
    msplit = 50,              # Split every 50 steps
    R_c = 100.0,              # Splitting radius
    p_max = 2.0,              # Maximum momentum
    n_p = 50                  # Momentum grid points
)

results = run_simulation(params)

# Access momentum distribution
p = results.momentum_grid.p
prob_p = sum(abs2.(results.psai_p), dims=(1,2))[1,1,:]

# Plot: p vs prob_p
```

### 3. HHG Spectrum Calculation
```julia
params = create_default_params(
    atom = :neon,
    laser_intensity = 5.0e13,
    laser_wavelength = 800.0,
    t_total = 300.0,
    enable_hhg = true,
    hhg_method = :acceleration,
    hhg_record_interval = 1
)

results = run_simulation(params)

# Compute spectrum
ω₀ = wavelength_nm_to_frequency_au(800.0)
h_orders, power, freqs = compute_hhg_spectrum(
    results.observables.dipole_moment,
    params.dt,
    ω₀,
    window=:hann
)

# Plot: h_orders vs log10(power)
```

### 4. Custom Soft-Core Potential
```julia
params = create_default_params(
    atom = :custom,
    custom_potential = "-1/sqrt(r^2 + 0.5)",  # a = 0.5
    laser_intensity = 5.0e14,
    t_total = 100.0
)

results = run_simulation(params)
```

### 5. Parameter Scan
```julia
# Create scan TOML (see examples/intensity_scan.toml)
config = load_config("intensity_scan.toml")

# Run scan
results = run_parameter_scan(config, verbose=true, parallel=false)

# Extract data
intensities = results.parameter_values
ionizations = [r.observables.ionization_probs[end]
               for r in results.simulation_results]
```

---

## Common Pitfalls

### 1. Grid Too Small
**Symptom**: Wavefunction reflects from boundary, norm > 1.0
**Fix**: Increase `rmax` or reduce `laser_intensity`

### 2. Time Step Too Large
**Symptom**: Norm increases dramatically, instability
**Fix**: Reduce `dt` (try 0.05 instead of 0.1)

### 3. Custom Potential → NaN
**Symptom**: Results are NaN/Inf
**Cause**: GPS eigenstates for custom potentials have numerical errors
**Status**: Known limitation, pipeline runs but results invalid

### 4. HHG Dipole is Zero
**Symptom**: `dipole_moment` all zeros
**Cause**: No l=1 eigenstates available (l=0 only mode for stability)
**Status**: Known limitation with analytical ground state workaround

---

## File Formats

### TOML Configuration
```toml
[atom]
type = "hydrogen"  # or "helium", "argon", "neon", "xenon", "custom"

[laser]
wavelength = 800.0
intensity = 1.0e14
pulse_duration = 50.0
polarization = [0.0, 0.0, 1.0]

[calculation]
type = "ionization"
output_file = "results.h5"
verbosity = "normal"

[numerical]
nrmax = 200
rmax = 100.0
lmax = 2
dt = 0.1
t_total = 200.0

[scan]  # Optional
parameter = "laser_intensity"
range = "1e13:1e13:1e15"
```

### HDF5 Output Structure
```
results.h5
├── /observables
│   ├── times [vector]
│   ├── norms [vector]
│   ├── energies [vector]
│   └── ionization_probs [vector]
├── /wavefunction
│   └── g_final [3D array: radial × m × l]
├── /momentum (if enabled)
│   ├── p_grid [vector]
│   └── psai_p [3D array: φ × θ × p]
└── /parameters
    └── [all simulation parameters as attributes]
```

---

## Performance Tips

1. **Start small**: Test with `nrmax=100, lmax=2, t_total=50` before full runs
2. **Disable analysis**: Set `enable_momentum_analysis=false, enable_hhg=false` for faster testing
3. **Parameter scans**: Use `parallel=true` with `addprocs(n)` for multi-core
4. **Memory**: Reduce `lmax` if running out of memory (scales as O(lmax²))
5. **Accuracy vs Speed**:
   - Fast: `nrmax=100, lmax=2, dt=0.2`
   - Balanced: `nrmax=200, lmax=10, dt=0.1`
   - Accurate: `nrmax=400, lmax=50, dt=0.05`

---

## Testing

```bash
# Quick tests (< 1 minute each)
julia tests/test_gps_grid.jl
julia tests/test_potential.jl
julia tests/test_hamiltonian.jl

# Full pipeline tests (2-5 minutes each)
julia tests/test_full_propagation_with_field.jl
julia tests/test_full_momentum_pipeline.jl
julia tests/test_custom_potential.jl
julia tests/test_parameter_scan.jl
```

---

## Troubleshooting

### Module Not Found
```julia
# Ensure correct path
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
include("src/TDSESolver.jl")
```

### Simulation Hangs
- Check `t_total` is reasonable (50-500 a.u. typical)
- Reduce `nrmax` for testing
- Check `dt` not too small

### Results Look Wrong
- **By design**: Physical accuracy not validated yet
- Check docs/UNSOLVED_PROBLEMS.md for known issues
- For built-in atoms (H, He): Should be reasonable for ionization
- For custom potentials: Expect NaN (GPS eigenstate issues)

---

## Getting Help

1. Check `docs/PROJECT_COMPLETION_SUMMARY.md` for overview
2. See `docs/UNSOLVED_PROBLEMS.md` for known limitations
3. Review example files in `examples/`
4. Examine test files in `tests/` for usage patterns
5. Read inline documentation in source files

---

**Last Updated**: January 2025
