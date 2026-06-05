# TDSE Refactoring Constitution
<!-- Refactoring Fortran TDSE Solver to Modern Julia -->

## Core Principles

### I. Numerical Accuracy First (NON-NEGOTIABLE)
All refactored code must maintain or exceed the numerical accuracy of the original Fortran implementations. Every module must:
- Validate against Fortran benchmarks with error tolerance < 1e-12 for bound states, < 1e-10 for continuum states
- Use double precision (Float64) by default; quad precision (Float128) for critical phase calculations
- Implement compensated summation (Kahan) for accumulation operations
- Document all numerical approximations and convergence criteria
- Preserve mathematical formulations exactly as described in reference papers (PRA articles)

### II. Modular Architecture
Code organization follows the 10-module architecture defined in refactoring_plan.tex:
- Each module is self-contained with clear interfaces and single responsibility
- Modules must be independently testable without requiring full simulation setup
- Public APIs documented with mathematical foundations, input/output specifications
- Private implementation details hidden; expose only necessary functions
- No circular dependencies between modules; dependency graph must be acyclic

### III. Fortran Blueprint Fidelity
Original Fortran programs (`D_inner_out_volkov_3d_with_prob.f90`, `rescatteing+hhg-he.f90`) serve as authoritative references:
- Preserve all core algorithms: split-operator method, GPS discretization, region splitting
- Maintain physical correctness of time propagation, field interactions, observable calculations
- Consult scientific articles when algorithm intent is unclear
- Do NOT "improve" physics without explicit validation against published benchmarks
- Use refactoring_plan.tex cautiously; verify against actual Fortran code when discrepancies exist

### IV. Modern Libraries Over Reinvention
Eliminate hardcoded implementations where mature Julia libraries exist:
- Linear algebra: Use LinearAlgebra.jl (BLAS/LAPACK) for matrix operations
- Special functions: Use SpecialFunctions.jl for Bessel, Legendre functions
- FFT: Use FFTW.jl for HHG spectrum calculations
- Quadrature: Use FastGaussQuadrature.jl for Gauss-Legendre grids
- I/O: Use HDF5.jl for large wavefunction storage
- Only implement custom numerics when no suitable library exists AND justify in documentation

### V. Performance Preservation
Refactored Julia code must maintain computational efficiency comparable to Fortran:
- Parallelize over independent angular momentum blocks (l values) using multi-threading
- Vectorize operations where possible; avoid scalar loops
- Pre-allocate arrays; minimize allocations in time evolution loop
- Profile performance bottlenecks; document optimization decisions
- Target: < 20% performance degradation compared to optimized Fortran (acceptable for development phase)
- Provide benchmark suite comparing timings for standard test cases

## Naming Conventions and Code Clarity (NON-NEGOTIABLE)

### Physically Meaningful Variable Names
The original Fortran codes use cryptic variable names (e.g., `g`, `s`, `plx`, `rp`, `w1`) that obscure physical meaning. The Julia refactoring **must** use descriptive names:

**Forbidden (Fortran style):**
```julia
g = zeros(ComplexF64, nr, -lmax:lmax, 0:lmax)
s = compute_matrix(...)
plx = legendre_zeros(nr)
```

**Required (Julia style):**
```julia
radial_wavefunction = zeros(ComplexF64, nr, -lmax:lmax, 0:lmax)
propagator_matrix = compute_time_evolution_matrix(...)
legendre_collocation_points = find_legendre_zeros(nr)
```

### Variable Naming Guidelines

**Wavefunction components:**
- `ψ` or `wavefunction` (not `g`)
- `radial_wavefunction` for radial part (not `g` or `temp`)
- `angular_wavefunction` for grid representation (not `psai_sp`)
- `momentum_amplitude` for momentum space (not `c_kl`)

**Grid quantities:**
- `radial_grid` (not `r`)
- `radial_derivative_mapping` (not `rp` for dr/dx)
- `quadrature_weights` (not `w1`)
- `legendre_points` (not `plx`)
- `theta_grid`, `phi_grid` (not `thi_sp`, `fi_sp`)

**Physical quantities:**
- `electric_field` (not `et`)
- `vector_potential` (not `A_et`)
- `field_amplitude` (not `e0`)
- `laser_frequency` (not `wmga`)
- `pulse_duration` (not `tp`)

**Operators and matrices:**
- `hamiltonian_matrix` (not `H`)
- `propagator_matrix` or `time_evolution_matrix` (not `s` or `S`)
- `kinetic_energy_matrix` (not generic variable)
- `potential_energy` (not `v` or `vr`)

**Eigenstate data:**
- `eigenvalues` (not `val`)
- `eigenvectors` (not `vec`)
- `bound_state_energies`, `continuum_energies` (not just `val`)

**Angular momentum:**
- `angular_momentum_l` (not `l` alone in isolated contexts)
- `magnetic_quantum_number_m` (not `m`)
- `legendre_polynomial` (not `plend` or `p`)
- `spherical_harmonic` (not `plgd`)

### Consistency Across Modules
Variable names must be **consistent** throughout the codebase:
- If radial wavefunction is `radial_wavefunction` in Grid module, use same name in Propagator module
- Do not alternate between `ψ`, `wavefunction`, `psi` for the same physical quantity
- Establish naming dictionary in project documentation
- Use module prefixes when needed: `Grid.radial_points`, `Propagator.time_step`

### Mathematical Symbols in Code
For well-established mathematical symbols, Julia allows Unicode:
- `ψ` for wavefunction (acceptable alternative to `wavefunction`)
- `θ`, `φ` for angles (acceptable alternatives to `theta`, `phi`)
- `Δt` for time step (acceptable alternative to `time_step`)
- `ω` for frequency (acceptable alternative to `frequency`)
- **But:** Always provide ASCII alias for accessibility: `const psi = ψ`

### Constants and Parameters
Physical constants must have clear names:
- `ionization_potential` (not `ip`)
- `ponderomotive_energy` (not `up`)
- `keldysh_parameter` (not `gama`)
- `MAX_ANGULAR_MOMENTUM` (not `lmax`)
- `NUM_RADIAL_POINTS` (not `nrmax`)
- `GRID_MAX_RADIUS` (not `rmax`)

### Loop Indices
Use descriptive loop indices:
```julia
# Forbidden
for i = 1:nr, j = 1:nth, k = 1:nphi
    # What do i, j, k represent?
end

# Required
for radial_idx = 1:num_radial_points
    for theta_idx = 1:num_theta_points
        for phi_idx = 1:num_phi_points
            # Clear iteration structure
        end
    end
end
```

Exception: Simple mathematical iterations where convention is clear:
```julia
for n = 1:num_eigenstates  # n is standard quantum number
    for l = 0:max_angular_momentum  # l is standard angular momentum
        for m = -l:l  # m is standard magnetic quantum number
```

### Comments and Documentation
Variable names should be self-documenting, but complex quantities require comments:
```julia
# Region splitting function F_s(r) from PRA 74, 031405(R) (2006), Eq. 5
# Smoothly separates inner (r < R_c) and outer (r > R_c) regions
split_function = 1.0 ./ (1.0 .+ exp.(-(radial_grid .- critical_radius) ./ smoothness_parameter))
```

### Enforcement
Code reviews must verify:
- No cryptic single-letter variables (except standard quantum numbers n, l, m)
- No abbreviations without clear physical meaning
- Consistent naming across all modules
- Comments explain physical interpretation where variable name alone is insufficient

## Validation Requirements

### Benchmark Test Suite (NON-NEGOTIABLE)
Every merge must pass validation against Fortran reference results:
- **Hydrogen ground state**: Energy converged to -0.5 Ha ± 1e-10
- **Helium ground state**: Energy converged to literature value ± 1e-8
- **Field-free evolution**: Norm conservation > 0.9999 after 1000 time steps
- **Strong-field ionization**: Match Fortran momentum distributions (correlation > 0.999)
- **HHG spectrum**: Harmonic peak positions and cutoff within 1% of Fortran
- Test data stored in `tests/benchmarks/` with expected outputs from Fortran runs

### Physical Consistency Checks
All simulations must satisfy:
- Norm conservation (accounting for absorbing boundaries)
- Energy conservation in field-free regions
- Symmetry preservation (e.g., cylindrical symmetry for linearly polarized fields)
- HHG cutoff law: ωmax ≈ Ip + 3.17Up
- Selection rules: only odd harmonics for symmetric systems

## Development Workflow

### Test-Driven Development
For new features or modules:
1. Write unit tests defining expected behavior (mathematical properties, known solutions)
2. Obtain user/reviewer approval of test specifications
3. Run tests → verify they fail initially (Red)
4. Implement feature until tests pass (Green)
5. Refactor for clarity/performance while maintaining passing tests

### Integration Testing Focus
Critical areas requiring integration tests:
- Module boundaries: Grid → Hamiltonian → Eigenstate → Propagator chains
- Time evolution loop: Full propagation with all modules active
- Region splitting: Inner/outer decomposition and Volkov projection
- Observable extraction: Populations, dipole moments, momentum distributions
- I/O roundtrip: Write wavefunction → read back → verify bit-identical

### Code Review Requirements
All pull requests must include:
- Evidence of benchmark validation (test output logs)
- Documentation updates (if public APIs changed)
- Performance profile comparison (if performance-critical code modified)
- Reviewer checklist: numerical accuracy, architecture compliance, test coverage

## Configuration and Extensibility

### Input Parameter System
Use structured configuration files (YAML/JSON format) instead of hardcoded constants:
- Grid parameters: Nr, Nθ, Nφ, rmax, L, α
- Field parameters: E0, ω, envelope type, polarization
- Propagation parameters: Δt, absorber settings, region split parameters
- Observable selection: which outputs to compute and save
- Example configurations provided in `examples/` directory

### Extension Points
Design for future enhancements:
- Pluggable potential models (hydrogen, helium, soft-core, custom)
- Multiple field configurations (two-color, chirped pulses)
- Alternative propagation schemes (Runge-Kutta, Magnus)
- Different continuum bases (plane waves, Coulomb-Volkov, Siegert states)
- Use Julia's multiple dispatch for algorithm variants

## Documentation Standards

### Required Documentation
Each module must provide:
- **Purpose**: One-paragraph description of module responsibility
- **Mathematical foundation**: Key equations implemented (LaTeX in docstrings)
- **Interface specification**: Function signatures, parameter descriptions, return types
- **Usage examples**: Minimal working code snippets
- **References**: Citations to scientific papers for algorithms

### Inline Documentation
- All public functions: docstrings with `@doc` macro
- Complex algorithms: explanatory comments referencing equation numbers from papers
- Non-obvious optimizations: justify why code structure differs from naive implementation
- Physical units: specify atomic units vs. SI units explicitly

## Governance

### Constitution Authority
This constitution supersedes all other development practices and style guides.

### Amendment Process
Constitution changes require:
1. Written proposal with rationale in GitHub issue
2. Discussion period (minimum 1 week)
3. Approval from project lead
4. Migration plan for existing code affected by changes
5. Version bump following semantic versioning

### Compliance Verification
All code reviews must verify:
- Numerical accuracy requirements met (Principle I)
- Modular architecture preserved (Principle II)
- Fortran blueprint fidelity maintained (Principle III)
- Appropriate library usage (Principle IV)
- Performance benchmarks pass (Principle V)
- Naming conventions followed (Naming Conventions section)
- Variable names physically meaningful and consistent across modules
- Test coverage adequate (Development Workflow)

### Complexity Justification
Any deviation from principles requires explicit justification in:
- Pull request description
- Code comments at deviation point
- Documentation update explaining rationale

**Version**: 1.1.0 | **Ratified**: 2025-01-21 | **Last Amended**: 2025-01-21
