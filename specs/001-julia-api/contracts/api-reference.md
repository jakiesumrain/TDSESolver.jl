# API Reference: TDSESolver.jl

**Version**: 0.1.0 | **Julia**: 1.10+

## High-Level API (User-Facing)

### Configuration Loading

#### `load_config(filepath::String) -> SimulationConfiguration`

Parses and validates a TOML configuration file.

**Parameters**:
- `filepath`: Path to TOML configuration file

**Returns**: `SimulationConfiguration` struct

**Validation** (fail-fast at load time):
- All physical parameters > 0
- Polarization vector normalized
- Calculation type valid (:ionization or :hhg)
- Custom potential specified if atom_type == :custom
- Parameter scan ranges valid

**Throws**: `ArgumentError` with descriptive message if validation fails

**Example**:
```julia
using TDSESolver

config = load_config("helium_800nm.toml")
```

---

### Calculation Execution

#### `run_tdse(config::SimulationConfiguration) -> CalculationResults`

Executes TDSE calculation with given configuration.

**Algorithm** (from Fortran blueprints):
1. Initialize physical system (load potential, compute ground state)
2. Setup computational grid (GPS method, angular grids)
3. Compute S-matrix for time evolution
4. Time propagation loop:
   - Split-operator propagation (3-step)
   - Apply absorbing boundaries
   - Monitor diagnostics (norm, ionization probability)
   - Checkpoint if configured
5. Project onto momentum eigenstates (Volkov) or compute HHG spectrum
6. Validate results and generate report

**Returns**: `CalculationResults` struct

**Progress Monitoring**: Outputs to terminal (default) and optional log file

**Example**:
```julia
results = run_tdse(config)
```

---

### Parameter Scans

#### `run_scan(config::SimulationConfiguration) -> Dict{Any, CalculationResults}`

Executes parameter scan specified in configuration.

**Parallelization**: Automatically distributes independent calculations across available threads

**Returns**: Dictionary mapping parameter values to results

**Example**:
```julia
scan_results = run_scan(config)

# Access specific result
results_at_5e14 = scan_results[5e14]
```

---

### Result Loading

#### `load_result(filepath::String) -> CalculationResults`

Loads calculation results from HDF5 file.

**Parameters**:
- `filepath`: Path to HDF5 results file

**Returns**: `CalculationResults` struct (lazy-loads large arrays)

**Example**:
```julia
results = load_result("helium_800nm_results.h5")
P_px_py = results.momentum_distribution_2D_xy
```

---

## Module-Level API (Advanced Users)

### Grid Module (`TDSESolver.Grid`)

#### `create_gps_grid(nrmax, rmax, L, α) -> GPSGrid`

Creates Generalized Pseudospectral grid with algebraic mapping.

**Parameters**:
- `nrmax::Int`: Number of radial points
- `rmax::Float64`: Maximum radius (a.u.)
- `L::Float64`: Mapping parameter
- `α::Float64`: Smoothness parameter

**Returns**: `GPSGrid` struct containing:
- `radial_grid`: r_i = L(1+x_i)/(1-x_i+α)
- `radial_derivative`: dr/dx
- `quadrature_weights`: Gauss-Legendre weights

**Algorithm**: From Tong & Chu 1997, uses FastGaussQuadrature.jl

---

#### `create_angular_grid(nthmax, nphimax) -> AngularGrid`

Creates θ, φ grids for 3D calculations.

**Parameters**:
- `nthmax::Int`: Number of θ points
- `nphimax::Int`: Number of φ points

**Returns**: `AngularGrid` struct with grids and weights

---

### Hamiltonian Module (`TDSESolver.Hamiltonian`)

#### `get_potential(atom_type::Symbol) -> Function`

Returns potential function V(r) for built-in atoms.

**Supported atoms**:
- `:hydrogen`: V(r) = -1/r
- `:helium`: V(r) = -Z_eff/r (with screening)
- `:argon`, `:neon`, `:xenon`: Effective potentials

**Returns**: Function `V(r::Float64) -> Float64`

---

#### `parse_custom_potential(expression::String) -> Function`

Parses Julia expression string into potential function.

**Example**:
```julia
# Soft-core potential
V = parse_custom_potential("-1/sqrt(r^2 + 0.5)")
```

---

### Eigenstate Module (`TDSESolver.Eigenstate`)

#### `compute_ground_state(grid, potential, n, l) -> (energy, wavefunction)`

Computes bound state using imaginary time propagation or diagonalization.

**Parameters**:
- `grid::GPSGrid`: Radial grid
- `potential::Function`: V(r)
- `n::Int`: Principal quantum number
- `l::Int`: Angular momentum

**Returns**: Tuple of (energy in Ha, wavefunction on grid)

**Validation**: Checks energy against literature values

---

#### `compute_smatrix(grid, eigenvalues, eigenvectors, Δt, lmax) -> Array{ComplexF64, 3}`

Computes time evolution operator in energy representation.

**Formula** (from Fortran, lines 255-293):
```
S_{l,ij} = Σₙ φ_{n,l}(xi) φ_{n,l}(xj) exp(-iE_{n,l}Δt/2) * wⱼr'ⱼ
```

**Returns**: S-matrix of shape (nrmax, nrmax, 0:lmax)

**Performance**: Pre-computed once, reused throughout propagation

---

### Field Module (`TDSESolver.Field`)

#### `create_laser_pulse(wavelength, intensity, duration, polarization) -> LaserPulse`

Creates time-dependent electric field function.

**Envelope**: Sin²((π/2)(t/T_pulse)) for smooth turn-on/off

**Returns**: `LaserPulse` struct with methods:
- `electric_field(t)`: E(t) vector (a.u.)
- `vector_potential(t)`: A(t) vector (a.u.)

---

### Propagator Module (`TDSESolver.Propagator`)

#### `propagate_split_operator!(ψ, S_matrix, E_field, grid, Δt)`

Applies one time step of split-operator propagation (in-place).

**Algorithm** (3-step, from Fortran lines 606-736):
1. Apply S-matrix: ψ → exp(-iĤ₀Δt/2) ψ
2. Transform to angular rep, apply field interaction: exp(-ir·E Δt)
3. Transform back, apply S-matrix again

**Modifies**: `ψ` in-place

---

### Continuum Module (`TDSESolver.Continuum`)

#### `project_volkov(ψ_outer, grid, A_field, p_grid, lmax) -> momentum_amplitudes`

Projects outer wavefunction onto Volkov states.

**Algorithm** (from Fortran lines 909-923, PRA 74.031405):
1. Bessel expansion: C_kl(m,l,p) = Σᵣ ψ(r,m,l) * j_l(pr) * weights
2. Volkov phase accumulation: exp(-i∫(p²/2 + A²/2 - p·A)dt)

**Returns**: momentum amplitudes C(p, l, m)

---

### Observable Module (`TDSESolver.Observable`)

#### `compute_momentum_distribution(momentum_amplitudes, p_grid) -> distributions`

Computes P(p), P(px,py), P(px,pz), P(px,py,pz) from momentum amplitudes.

**Formula**: P(p) = |C(p, l, m)|²

**Returns**: Named tuple with 1D, 2D, 3D distributions

---

#### `compute_hhg_spectrum(dipole_history, Δt) -> (photon_energies, spectrum)`

Computes HHG power spectrum via FFT of dipole moment.

**Algorithm**:
1. Apply Hann window to dipole(t)
2. FFT using FFTW.jl
3. Compute |d(ω)|²
4. Identify harmonic orders and cutoff

**Returns**: Tuple of (ω in eV, S(ω))

---

### IO Module (`TDSESolver.IO`)

#### `save_results(results::CalculationResults, filepath::String)`

Saves results to HDF5 file with metadata.

**HDF5 Structure**:
```
/momentum_distributions/
    P_radial, P_px_py, P_px_pz, grids...
/hhg_spectrum/ (if applicable)
    spectrum, photon_energy, harmonic_orders...
/validation_diagnostics/
    norm_history, energy_conservation, correlation...
/metadata/ (attributes)
    wavelength, intensity, timestamp, versions...
```

---

#### `query_result(filepath::String, path::Symbol) -> Array`

Queries specific dataset from HDF5 file without loading entire file.

**Example**:
```julia
P_xy = query_result("results.h5", :momentum_px_py)
```

---

## Validation Module (`TDSESolver.Validation`)

#### `validate_calculation(results::CalculationResults) -> ValidationReport`

Generates validation report checking:
- Norm conservation (> 0.9999)
- Energy conservation (field-free periods)
- Symmetry (if applicable)
- Benchmark correlation (if Fortran data available)
- HHG cutoff law (I_p + 3.17U_p, within 5%)

**Returns**: `ValidationReport` with pass/fail status and warnings

---

## Utility Functions

#### `PhysicalUnits.au_to_eV(energy_au::Float64) -> Float64`

Converts atomic units to electron volts (× 27.2114).

#### `PhysicalUnits.intensity_SI_to_au(I_W_cm2::Float64) -> Float64`

Converts W/cm² to atomic units of intensity.

#### `PhysicalUnits.wavelength_nm_to_frequency_au(λ_nm::Float64) -> Float64`

Converts wavelength (nm) to frequency (a.u.).

---

## Error Handling

All functions throw `ArgumentError` for invalid inputs with descriptive messages. Configuration loading implements fail-fast validation per clarification Q3.

**Common errors**:
- `"Wavelength must be positive"`: λ ≤ 0
- `"Grid size insufficient for ponderomotive energy U_p = X Ha"`: nrmax too small
- `"Time step Δt = X exceeds stability criterion π/Emax = Y"`: Δt too large
- `"Custom potential required for atom_type = :custom"`: Missing potential specification

---

## Performance Optimization

For optimal performance:
1. Use MKL backend: `using MKL` before `using TDSESolver`
2. Set threads: `export JULIA_NUM_THREADS=8`
3. Pre-allocate for scans: Reuse grid/S-matrix across parameter values

**Expected performance**: < 20% slower than Fortran on same hardware (Constitution Principle V)

---

## References

- **Fortran blueprints**: `explore/D_inner_out_volkov_3d_with_prob.f90`, `explore/rescatteing+hhg-he.f90`
- **Scientific papers**: Tong & Chu 1997 (split-operator), PRA 74.031405(R)(2006) (region splitting)
- **Constitution**: `.specify/memory/constitution.md` (Core Principles I-V)
