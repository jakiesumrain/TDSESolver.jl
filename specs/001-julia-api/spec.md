# Feature Specification: Research-Grade TDSE Solver with User-Friendly Julia API

**Feature Branch**: `001-julia-api`
**Created**: 2025-01-21
**Status**: Draft
**Input**: User description: "This project will be implemented in modern julia facilities for maximum performance. The final goal is a research-grade TDSE solver. And users can use the solver with minimal setup, which means that there should exist a well-documented APIs that the user can manipulate without knowing the details of the actual algorithms behind."

## Clarifications

### Session 2025-01-21

- Q: How should users specify simulation parameters (wavelength, intensity, pulse duration, time step, etc.)? → A: Structured configuration files (TOML/YAML) that users create/edit, then loaded via API (e.g., `config = load_config("helium_800nm.toml"); run_tdse(config)`)

- Q: How should users specify custom potentials in configuration files? → A: Mathematical expressions in Julia syntax in config files, plus direct Julia function or module representing the potential that can be modified by experienced users

- Q: When should parameter validation occur - at configuration load time or just before starting calculation? → A: At configuration load time: `load_config()` parses file AND validates all parameters, failing immediately with detailed error messages if any parameter is invalid

- Q: How should users monitor calculation progress (FR-008, User Story 4)? → A: Both: terminal output (default) with optional structured log files for batch processing, user can configure verbosity level

- Q: How should parameter ranges for scans be defined in configuration files? → A: Separate scan section in config file with range notation (e.g., `[scan]` section with `parameter = "intensity"`, `range = "1e14:1e14:5e14"`) - clear separation of single vs. scan calculations

- Q: What is the role of Fortran programs and scientific papers in implementation? → A: Fortran programs (`D_inner_out_volkov_3d_with_prob.f90`, `rescatteing+hhg-he.f90`) are authoritative blueprints and starters - implementation must follow their algorithms, not create new ones. Scientific papers (Tong & Chu 1997, PRA articles) must be consulted frequently to verify correctness

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Basic Ionization Calculation (Priority: P1)

A researcher needs to compute photoelectron momentum distributions for a helium atom in a strong laser field with minimal setup time.

**Why this priority**: This represents the core value proposition - users can run production-quality TDSE calculations without understanding GPS methods, split-operator propagation, or Volkov projections. This is the minimum viable product that delivers immediate research value.

**Independent Test**: Can be fully tested by providing laser parameters (wavelength, intensity, pulse duration) and atom type, then verifying that physically correct momentum distributions are produced that match published benchmark results.

**Acceptance Scenarios**:

1. **Given** a researcher has basic knowledge of laser-atom physics, **When** they specify laser parameters (800 nm, 5×10¹⁴ W/cm², 10 fs pulse), atom type (helium), and request ionization calculation, **Then** the system computes photoelectron momentum distributions without requiring the user to understand numerical grid construction, time propagation methods, or region splitting algorithms.

2. **Given** a user provides simulation parameters, **When** computation completes, **Then** the system returns momentum distributions in standard formats (1D radial, 2D slices, 3D if requested) with clear physical units (atomic units) and coordinate systems clearly documented.

3. **Given** a user specifies calculation parameters, **When** parameters are outside physically reasonable ranges, **Then** the system provides clear warnings with physically meaningful explanations (e.g., "Grid size insufficient for ponderomotive energy" rather than "nrmax too small").

---

### User Story 2 - High-Harmonic Generation Calculation (Priority: P2)

A researcher needs to compute HHG spectra for different atoms and laser configurations to understand harmonic cutoff behavior.

**Why this priority**: HHG calculations are a distinct use case requiring different observables (dipole moment/acceleration) and post-processing (FFT). This extends the solver's utility to a second major research domain while reusing the same underlying TDSE infrastructure.

**Independent Test**: Can be fully tested by running HHG calculations with known laser-atom combinations and verifying that harmonic spectra match published results, including correct cutoff position (Ip + 3.17Up) and odd-harmonic structure.

**Acceptance Scenarios**:

1. **Given** a researcher wants HHG spectra, **When** they specify calculation type as "HHG" with laser and atom parameters, **Then** the system computes time-dependent dipole moment and produces HHG power spectrum without user needing to understand dipole operator construction or windowing functions for FFT.

2. **Given** a user requests HHG calculation, **When** computation completes, **Then** the system provides harmonic spectrum with photon energy axis, identifies harmonic orders, and indicates cutoff position with physical interpretation.

3. **Given** a researcher compares different laser intensities, **When** multiple HHG calculations are requested with varying intensity, **Then** the system produces consistent results showing expected intensity scaling of cutoff (Ip + 3.17Up where Up ∝ intensity).

---

### User Story 3 - Custom Potential Models (Priority: P3)

A researcher studying exotic atoms or molecules needs to define custom interaction potentials beyond standard hydrogen/helium models.

**Why this priority**: This provides research flexibility for advanced users while maintaining the simple interface for standard cases. Custom potentials are important for specialized research but not required for the majority of users studying standard atomic systems.

**Independent Test**: Can be tested independently by providing a custom potential function (e.g., soft-core potential) and verifying that ground state energies and ionization dynamics match analytical or published results for that specific potential.

**Acceptance Scenarios**:

1. **Given** a researcher has a custom potential model, **When** they provide a potential via mathematical expression in configuration file (e.g., `potential = "-1/sqrt(r^2 + 0.5)"` for soft-core potential) or as Julia function in separate module, **Then** the system validates physical properties (correct asymptotic behavior, binding energy) and integrates it into calculations without user modifying core solver code.

2. **Given** a user defines custom potential, **When** initial state calculation is performed, **Then** the system computes ground state wavefunction and energy for verification, reporting convergence diagnostics in physically interpretable terms.

---

### User Story 4 - Parameter Exploration and Optimization (Priority: P4)

A researcher needs to run systematic parameter scans (intensity, wavelength, pulse duration) to explore ionization yield trends or optimize HHG efficiency.

**Why this priority**: Once basic calculations work (P1-P3), researchers need tools for systematic exploration. This is lower priority as users could manually script parameter loops using the P1 API, but a native scan capability significantly improves usability.

**Independent Test**: Can be tested by defining a parameter sweep (e.g., intensity from 10¹⁴ to 10¹⁵ W/cm² in 10 steps) and verifying that results show expected physical trends with proper parallel execution and result aggregation.

**Acceptance Scenarios**:

1. **Given** a researcher wants to scan laser intensity, **When** they specify parameter range in configuration file's `[scan]` section (e.g., `parameter = "intensity"`, `range = "1e14:1e14:1e15"` for 10 points), **Then** the system executes calculations in parallel when possible and aggregates results into structured output (tables, arrays) for analysis.

2. **Given** a user requests parameter scan, **When** calculations are running, **Then** the system provides progress monitoring via terminal output (default) and optional structured log files, showing completed/remaining calculations and estimated time remaining in human-readable format, with configurable verbosity level.

3. **Given** a parameter scan encounters physical issues (e.g., grid inadequate at highest intensity), **When** individual calculations fail validation checks, **Then** the system continues other calculations and reports which parameter combinations failed with physical reasoning.

---

### User Story 5 - Result Validation and Benchmarking (Priority: P5)

A new user wants to verify that the solver produces physically correct results before trusting it for novel research.

**Why this priority**: This builds user confidence and is essential for research-grade software, but is lower priority than core functionality. Benchmarks can initially be documented separately from the API itself.

**Independent Test**: Can be tested by running a suite of standard benchmark problems (hydrogen ionization, helium HHG) and automatically comparing results against reference data with quantitative similarity metrics.

**Acceptance Scenarios**:

1. **Given** a new user wants to validate solver accuracy, **When** they request benchmark calculations, **Then** the system runs standard test cases (hydrogen ground state, helium ionization at known intensity) and reports agreement with published results using quantitative metrics.

2. **Given** a user completes benchmark validation, **When** results are within expected tolerances, **Then** the system provides confidence report indicating solver is working correctly and ready for research use.

---

### Edge Cases

- **Insufficient grid resolution**: What happens when user-specified laser intensity produces ponderomotive energy requiring larger spatial grid than provided? System should detect this condition from Keldysh parameter and ponderomotive energy, warn user with physical explanation, and suggest grid parameters.

- **Unphysical parameters**: How does system handle negative intensities, zero wavelength, or other mathematically invalid inputs? System validates all physical parameters at configuration load time (`load_config()`) with clear error messages referencing physical constraints, failing immediately before any computation.

- **Numerical instability**: What happens if time step is too large for stable propagation? System should validate time step against stability criterion (Δt < π/Emax) and warn user, providing recommended value.

- **Memory limitations**: How does system handle cases where requested calculation exceeds available RAM (e.g., lmax=200, nr=1000)? System should estimate memory requirements before starting calculation and either warn user or automatically adjust parameters with user confirmation.

- **Contradictory parameters**: How does system handle physically inconsistent specifications (e.g., requesting ionization calculation but specifying field too weak to ionize)? System should detect when ionization probability will be negligible (Keldysh γ >> 1) and warn user that results may be dominated by numerical noise.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST accept simulation parameters via structured configuration files (TOML or YAML format) that specify atom type, laser parameters (wavelength, intensity, pulse duration, polarization), calculation type (ionization/HHG), and optional numerical overrides

- **FR-002**: System MUST provide high-level interface where users specify only physical parameters without needing to specify numerical details (grid sizes, time steps, angular momentum cutoffs), with automatic determination of these from physical parameters

- **FR-003**: System MUST support standard atomic systems (hydrogen, helium, argon, neon, xenon) with built-in potential models, requiring only atom name as specification

- **FR-004**: System MUST compute and return photoelectron momentum distributions in multiple representations: 1D radial P(p), 2D planar slices P(px,py), P(px,pz), and optionally full 3D P(px,py,pz)

- **FR-005**: System MUST compute and return HHG power spectra S(ω) with photon energy axis in electron-volts, identifying harmonic orders and cutoff position

- **FR-006**: System MUST validate all input parameters at configuration load time (in `load_config()`) against physical constraints and numerical stability criteria, failing immediately with clear error messages and physical interpretation when validation fails, before any expensive computation begins

- **FR-007**: System MUST provide comprehensive documentation including: quick-start guide, parameter reference with physical meanings, example calculations reproducing published results, and troubleshooting guide for common issues

- **FR-008**: System MUST report calculation progress via terminal output (default) and optional structured log files, with physically meaningful indicators (current simulation time, ionization probability, norm conservation) rather than only numerical details (iteration count), and configurable verbosity level

- **FR-009**: System MUST preserve numerical accuracy of original Fortran implementation by following the algorithms in Fortran blueprints (`D_inner_out_volkov_3d_with_prob.f90`, `rescatteing+hhg-he.f90`) as authoritative references, validating against benchmark calculations with maximum acceptable deviations: bound state energies within 10⁻¹⁰ Ha, momentum distributions with correlation > 0.999. Implementation must consult scientific papers (Tong & Chu 1997, PRA articles) frequently to verify algorithm correctness.

- **FR-010**: System MUST provide result objects with clear attribute names (e.g., result.momentum_distribution, result.radial_grid, result.kinetic_energy_spectrum) rather than requiring users to parse output files

- **FR-011**: System MUST include validation diagnostics in results: norm conservation history, energy conservation (field-free periods), symmetry checks (if applicable), with warnings if validation thresholds are exceeded

- **FR-012**: System MUST support parameter scans specified via `[scan]` section in configuration files using range notation (start:step:stop), automatically distributing independent calculations across available computational resources and aggregating results

- **FR-013**: System MUST allow users to specify custom potential functions via (a) mathematical expressions in Julia syntax within configuration files, or (b) direct Julia function/module references for advanced users, with automatic validation of physical properties (asymptotic behavior, boundedness)

- **FR-014**: System MUST provide example scripts demonstrating common use cases: basic ionization, HHG calculation, parameter scan, custom potential, result plotting

### Key Entities

- **Simulation Configuration**: Represents complete specification of a TDSE calculation loaded from structured configuration file (TOML/YAML), including atom type, laser parameters (wavelength, intensity, pulse shape, duration, polarization), calculation type (ionization/HHG), and optional numerical overrides. Configuration files enable parameter versioning, reproducibility, and batch processing.

- **Physical System**: Represents the atom-field system including atomic potential model, ground state properties (energy, quantum numbers), field properties (ponderomotive energy, Keldysh parameter)

- **Computational Grid**: Represents spatial and temporal discretization including radial grid points, angular grids (theta, phi), time step, mapping parameters - mostly transparent to users but accessible for advanced diagnostics

- **Wavefunction State**: Represents quantum state at specific time including radial components, norm, populations in bound states, ionization probability - accessible for monitoring but not requiring user manipulation

- **Calculation Results**: Represents final output of calculation including momentum distributions (1D/2D/3D), HHG spectrum (if applicable), validation diagnostics, metadata (parameters used, computation time), and methods for visualization and export

- **Validation Report**: Represents quality checks on calculation including norm conservation, energy conservation, benchmark comparisons (if applicable), warnings about potential issues, and confidence indicators

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Users with basic laser-atom physics knowledge can set up and run their first TDSE calculation in under 10 minutes from installation to viewing results

- **SC-002**: A standard ionization calculation (helium, 800 nm, 10 optical cycles) completes in under 30 minutes on typical research workstation (8-core desktop)

- **SC-003**: Calculated momentum distributions agree with published benchmark results with quantitative correlation > 0.995 for standard test cases (hydrogen ionization at multiple intensities)

- **SC-004**: At least 90% of common use cases (standard atoms, typical laser parameters, standard observables) require zero manual numerical parameter tuning by users

- **SC-005**: Documentation completeness allows users to successfully complete example calculations and understand physical meaning of all parameters without consulting external papers or asking questions in 80% of cases

- **SC-006**: System automatically detects and warns users about 95% of common pitfalls (insufficient grid, too-large time step, unphysical parameters) at configuration load time, failing immediately with actionable error messages before starting expensive calculations

- **SC-007**: Parameter scans complete with near-linear speedup on multi-core systems (efficiency > 85% on 8-core workstation)

- **SC-008**: Users can integrate solver into research workflows and reproduce published results from literature with total effort (setup + computation + validation) reduced by 70% compared to using original Fortran codes directly

- **SC-009**: HHG calculations produce harmonic spectra with correct cutoff position within 5% of theoretical prediction (Ip + 3.17Up) for all tested intensities

- **SC-010**: System handles edge cases gracefully with informative error messages, maintaining user experience quality where 90% of users can resolve common issues without external support

### Validation Metrics

- Numerical accuracy: Maximum deviation in bound state energies < 10⁻¹⁰ Ha from Fortran reference
- Correlation with Fortran results: Momentum distributions correlation > 0.999
- Performance: < 20% runtime increase compared to optimized Fortran for equivalent calculations
- API usability: Users can complete basic calculation with < 20 lines of code
- Documentation coverage: 100% of public API functions documented with examples

## Assumptions

- Users have basic understanding of atomic physics and laser-atom interactions (know terms like ionization potential, wavelength, intensity, ponderomotive energy)
- Users have access to computational resources suitable for TDSE calculations (8+ core workstation, 16+ GB RAM) for typical problems
- Standard calculation types (ionization, HHG) for standard atoms (H, He, Ar, Ne, Xe) cover 90% of user needs
- Users understand that TDSE calculations are computationally intensive and may take minutes to hours depending on parameters
- Users working with custom potentials have sufficient expertise to specify potential functions correctly and validate results
- Default numerical parameters (automatic determination) will be conservative, favoring accuracy over speed for initial implementation
- Users can install Julia and required packages following standard Julia package manager workflow
- Visualization of results (plotting) can use standard Julia plotting packages rather than being built into solver itself

## Dependencies

- **Fortran reference implementations**: Authoritative blueprint programs (`D_inner_out_volkov_3d_with_prob.f90` for ionization/ATI, `rescatteing+hhg-he.f90` for HHG) located in `explore/` directory - implementation must follow their algorithms, not create alternatives
- **Scientific papers**: Algorithm references that must be consulted frequently - Tong & Chu (1997) for split-operator method, PRA 74.031405(R)(2006) for region splitting, available in `explore/` directory
- Access to reference data from original Fortran codes for validation
- Availability of published benchmark results for standard test cases
- Julia numerical computing packages (FFTW.jl for HHG, SpecialFunctions.jl for special functions, LinearAlgebra.jl for matrix operations)
- Sufficient computational resources to run validation benchmarks during development
- Expert validation of automatically determined numerical parameters against established convergence criteria

## Out of Scope

- Multi-electron systems beyond single-active-electron approximation (e.g., full two-electron helium calculation)
- Molecular systems with multiple nuclei or non-spherical potentials
- Relativistic effects or magnetic field interactions
- Time-dependent potential models or two-color laser fields (initially - may be added in future versions)
- Interactive GUI for parameter specification (command-line/script interface is sufficient for research users)
- Built-in plotting and visualization (users can use standard Julia plotting packages with result data)
- Automatic publication-quality figure generation
- Database or cloud integration for result storage
- Comparison with other theoretical methods (perturbation theory, strong-field approximation) - this is user's responsibility

## Design Decisions

### Documentation Format

**Decision**: In-source docstrings accessible via Julia's help system (`? function_name`), with separate tutorial notebooks showing example workflows

**Rationale**:
- Most familiar to Julia users and integrates with standard Julia workflow
- Documentation lives with code, ensuring consistency and ease of maintenance
- Julia's built-in help system (`?`) is the standard way Julia users discover API details
- Tutorial notebooks (Pluto.jl or Jupyter) provide interactive learning for complex workflows
- Lower maintenance burden compared to separate documentation site hosting

**Implementation**: All public API functions will have comprehensive docstrings following Julia documentation standards, with tutorial notebooks in `examples/` directory demonstrating common research workflows.

---

### Result Format

**Decision**: Automatically saved to HDF5 files with metadata, with helper functions to load and query results

**Rationale**:
- Better suited for large batch calculations typical in research workflows
- Persistent by default - results survive Julia session crashes or interruptions
- HDF5 format is standard in scientific computing, accessible from multiple languages
- Efficient storage for large multi-dimensional arrays (momentum distributions)
- Metadata preservation ensures reproducibility (parameters, timestamps, versions)
- Helper functions can provide convenient data access patterns without requiring users to understand HDF5 structure

**Implementation**: Results automatically saved to `[calculation_name]_results.h5` with structured groups for momentum distributions, validation diagnostics, and metadata. Provide convenience functions like `load_result(filename)` and `query_result(filename, :momentum_px)` for common access patterns.

---

### Parameter Input Mechanism

**Decision**: Structured configuration files (TOML or YAML format) that users create and edit, loaded via API functions

**Rationale**:
- Enables parameter grouping and organization as requested by user
- Supports reproducibility - configuration files can be version-controlled alongside results
- Allows batch processing and parameter scans without code modification
- Standard practice in computational physics software (LAMMPS, Quantum ESPRESSO, VASP)
- Non-programmers can modify parameters without editing Julia code
- Configuration files serve as self-documenting parameter records for publications

**Implementation**: Users create configuration files (e.g., `helium_800nm.toml`) specifying physical parameters. API provides `load_config(filename)` to parse and validate configuration, and `run_tdse(config)` to execute calculation. Example configuration structure:
```toml
[atom]
type = "helium"

[laser]
wavelength = 800.0  # nm
intensity = 5e14     # W/cm²
pulse_duration = 10.0  # fs
polarization = [1.0, 0.0, 0.0]  # x, y, z components

[calculation]
type = "ionization"  # or "HHG"
output_file = "helium_800nm_results.h5"
```

---

### Custom Potential Specification

**Decision**: Two-level approach: (a) mathematical expressions in Julia syntax within configuration files, or (b) direct Julia function/module references for advanced users

**Rationale**:
- Mathematical expressions handle most use cases (soft-core, Yukawa, model potentials) without code
- Expressions in configuration files maintain parameter-code separation and reproducibility
- Direct Julia functions provide full control for complex potentials (coordinate-dependent, time-dependent, multi-parameter)
- Advanced users can leverage Julia's full language features (conditionals, special functions, optimization)
- Consistent with Julia ecosystem philosophy: simple for common cases, powerful for complex ones

**Implementation**:
- **Simple case** (configuration file): Parse Julia expressions like `potential = "-1/sqrt(r^2 + a^2)"` where `a` is a parameter. System evaluates expression at each grid point.
- **Advanced case** (module reference): User creates `MyPotentials.jl` with function `V(r; params...)` and specifies `potential_module = "MyPotentials"` and `potential_function = "V"` in config.
- System validates both: checks asymptotic behavior (r→∞ should approach 0 or -Z/r), boundedness, and computes test ground state.

Example configuration with custom potential:
```toml
[atom]
type = "custom"
potential = "-1/sqrt(r^2 + 0.5)"  # Soft-core potential, a=0.5
charge = 1  # Effective nuclear charge

# Alternative for advanced users:
# potential_module = "MyPotentials"
# potential_function = "custom_V"
```

---

### Progress Monitoring Interface

**Decision**: Dual-channel progress reporting: terminal output (default) with optional structured log files, configurable verbosity level

**Rationale**:
- Terminal output provides immediate feedback for interactive users (most common use case)
- Structured log files enable automated monitoring for batch jobs and long-running parameter scans
- Configurable verbosity allows users to adjust detail level (minimal/normal/verbose)
- Standard practice in scientific computing packages (GROMACS, LAMMPS, Quantum ESPRESSO)
- Log files can be parsed by monitoring scripts for distributed computing environments
- Terminal output suitable for Jupyter notebooks with live updates

**Implementation**:
- Default behavior: print progress to terminal/stdout with physically meaningful indicators (simulation time, ionization probability, norm conservation)
- Configuration option: `log_file = "helium_800nm_progress.log"` enables structured logging with timestamps
- Verbosity levels in config: `verbosity = "normal"` (default), `"minimal"` (warnings only), `"verbose"` (detailed diagnostics)
- Progress format includes estimated time remaining, percent complete, and key physical quantities
- Log files use structured format (JSON lines or CSV) for easy parsing by monitoring tools

---

### Parameter Scan Specification

**Decision**: Separate `[scan]` section in configuration files with range notation (start:step:stop)

**Rationale**:
- Clear separation between single-calculation and scan configurations
- Range notation (start:step:stop) is concise and familiar from Julia/MATLAB
- Allows specifying multiple scan dimensions (2D parameter spaces: intensity × wavelength)
- Keeps main parameter sections clean for single calculations
- Standard practice in computational chemistry (Gaussian scan jobs, ORCA surface scans)
- Configuration file clearly documents scan methodology for reproducibility

**Implementation**:
- Single calculations: omit `[scan]` section, all parameters are scalar values
- Scan calculations: add `[scan]` section specifying parameter(s) to vary and ranges
- Range syntax: `"start:step:stop"` generates array of values (Julia-style)
- Multi-dimensional scans: specify multiple parameters, system generates Cartesian product
- System automatically parallelizes independent parameter combinations
- Aggregated results saved to single HDF5 file with scan dimensions as dataset axes

Example scan configuration:
```toml
[atom]
type = "helium"

[laser]
wavelength = 800.0  # nm - fixed
intensity = 5e14     # W/cm² - will be overridden by scan
pulse_duration = 10.0  # fs - fixed
polarization = [1.0, 0.0, 0.0]

[calculation]
type = "ionization"
output_file = "helium_intensity_scan_results.h5"

[scan]
parameter = "laser.intensity"
range = "1e14:1e14:1e15"  # 10 points from 1e14 to 1e15 W/cm²
```

For 2D scans:
```toml
[scan]
parameters = ["laser.intensity", "laser.wavelength"]
ranges = ["1e14:5e13:5e14", "400:200:1200"]  # 10×5 grid
```

---

### Algorithm Fidelity to Fortran Blueprints and Scientific Papers

**Decision**: Fortran programs are authoritative blueprints - implementation follows their algorithms exactly, scientific papers are consulted frequently for verification

**Rationale**:
- Fortran programs (`D_inner_out_volkov_3d_with_prob.f90`, `rescatteing+hhg-he.f90`) represent validated, production-tested implementations
- These codes have been used in published research and represent years of domain expertise
- Creating alternative algorithms from scratch risks introducing physics errors or numerical instabilities
- Scientific papers (Tong & Chu 1997, PRA articles) provide mathematical foundations and validation criteria
- Constitution mandates "Fortran Blueprint Fidelity" as Core Principle III
- Project goal is modernization (better API, readability) NOT re-invention of physics algorithms

**Implementation Workflow**:
1. **Algorithm Extraction**: Study Fortran code section-by-section (e.g., split-operator propagation in lines 606-736 of `D_inner_out_volkov_3d_with_prob.f90`)
2. **Paper Consultation**: Cross-reference with scientific papers to understand mathematical foundations (e.g., Tong & Chu 1997 Eq. 15 for split-operator formula)
3. **Julia Translation**: Translate algorithm to Julia preserving numerical operations, using physically meaningful variable names
4. **Validation**: Compare outputs against Fortran benchmark results (correlation > 0.999 for momentum distributions)
5. **Documentation**: Document algorithm source (Fortran line numbers, paper equations) in code comments

**What to preserve from Fortran**:
- Split-operator time propagation sequence (3-step method)
- GPS grid construction and mapping formulas
- S-matrix computation and application
- Region splitting transition function
- Volkov state projection formulas
- Absorbing boundary implementation

**What to modernize**:
- Variable naming (Fortran `g` → Julia `radial_wavefunction`)
- Configuration input (hardcoded constants → TOML files)
- Result output (text files → HDF5 with metadata)
- Code organization (monolithic → modular)
- Parallelization (OpenMP → Julia native threading)

**Key References**:
- `explore/D_inner_out_volkov_3d_with_prob.f90` (2278 lines) - ionization/ATI blueprint
- `explore/rescatteing+hhg-he.f90` - HHG blueprint
- `explore/Tong and Chu - 1997.pdf` - split-operator method in energy representation
- `explore/PhysRevA.74.031405.pdf` - region splitting with smooth transition function
- `.specify/memory/constitution.md` - Core Principle III: Fortran Blueprint Fidelity
