# Implementation Plan: Research-Grade TDSE Solver with User-Friendly Julia API

**Branch**: `001-julia-api` | **Date**: 2025-01-21 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/001-julia-api/spec.md`

**Note**: This plan follows the constitution principles for refactoring Fortran TDSE solvers to modern Julia while preserving numerical accuracy and algorithm fidelity.

## Summary

Implement a research-grade Time-Dependent Schrödinger Equation (TDSE) solver in Julia with a user-friendly API that allows researchers to compute photoelectron momentum distributions and HHG spectra with minimal setup. The implementation follows algorithms from authoritative Fortran blueprints (`D_inner_out_volkov_3d_with_prob.f90`, `rescatteing+hhg-he.f90`) and scientific papers (Tong & Chu 1997, PRA articles), translating them to Julia with modern naming conventions, modular architecture, and TOML configuration input. Core technical approach: Extract algorithms from Fortran section-by-section → Consult scientific papers for mathematical foundations → Translate to Julia preserving numerical operations → Validate against Fortran benchmarks (correlation > 0.999).

## Technical Context

**Language/Version**: Julia 1.10+ (LTS recommended for stability)
**Primary Dependencies**:
- LinearAlgebra.jl (BLAS/LAPACK - standard library)
- SpecialFunctions.jl (Bessel, Legendre functions)
- FFTW.jl (HHG spectrum calculations)
- FastGaussQuadrature.jl (Gauss-Legendre grids)
- HDF5.jl (result storage)
- TOML.jl (configuration parsing - standard library)

**Storage**: HDF5 files for results (momentum distributions, validation diagnostics, metadata)
**Testing**: Julia Test framework (standard library) with benchmark suite comparing to Fortran outputs
**Target Platform**: Linux/macOS/Windows workstations (8+ core, 16+ GB RAM)
**Project Type**: Single Julia package (scientific computing library)
**Performance Goals**: < 20% runtime increase vs. optimized Fortran, > 85% parallel efficiency on 8-core systems
**Constraints**:
- Numerical accuracy: < 1e-12 error for bound states, < 1e-10 for continuum states
- Correlation with Fortran: > 0.999 for momentum distributions
- API simplicity: < 20 lines of code for basic calculation

**Scale/Scope**:
- 10 core modules (Grid, Hamiltonian, Eigenstate, Field, Propagator, RegionSplit, Continuum, Observable, IO, Orchestrator)
- ~2000-3000 lines refactored from Fortran
- 5 user stories (P1-P5), 14 functional requirements

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Core Principles Compliance

**I. Numerical Accuracy First** ✓
- Validation: Every module benchmarked against Fortran (< 1e-12 bound state error, < 1e-10 continuum error)
- Precision: Float64 default, Float128 for critical phase calculations
- Summation: Kahan summation for accumulation operations
- Documentation: All approximations and convergence criteria documented

**II. Modular Architecture** ✓
- Structure: 10 self-contained modules per refactoring_plan.tex
- Testing: Each module independently testable
- Interfaces: Public APIs documented with mathematical foundations
- Dependencies: Acyclic dependency graph enforced

**III. Fortran Blueprint Fidelity** ✓
- References: `D_inner_out_volkov_3d_with_prob.f90` (ionization), `rescatteing+hhg-he.f90` (HHG)
- Algorithms: Preserve split-operator, GPS, S-matrix, region splitting, Volkov projection
- Verification: Cross-reference with Tong & Chu 1997, PRA 74.031405(R)(2006)
- Approach: Translate, do NOT reinvent physics algorithms

**IV. Modern Libraries Over Reinvention** ✓
- LinearAlgebra.jl for matrix operations (replaces hardcoded BLAS calls)
- SpecialFunctions.jl for Bessel/Legendre (replaces hardcoded polynomial evaluations)
- FFTW.jl for HHG spectra (replaces custom FFT)
- FastGaussQuadrature.jl for Gauss-Legendre grids (replaces hardcoded zeros/weights)
- HDF5.jl for I/O (replaces text files)

**V. Performance Preservation** ✓
- Parallelization: Multi-threading over angular momentum blocks (l values)
- Vectorization: Array operations where possible
- Pre-allocation: Minimize allocations in time evolution loop
- Target: < 20% degradation vs. Fortran acceptable

**Naming Conventions** ✓
- Descriptive names: `radial_wavefunction` not `g`, `propagator_matrix` not `s`
- Physics clarity: `electric_field` not `et`, `laser_frequency` not `wmga`
- Consistency: Same names across all modules

### Constitution Verification Result

**STATUS**: ✅ **PASS** - All core principles addressed in design

No violations requiring justification. Implementation follows established refactoring blueprint with modern Julia ecosystem libraries.

## Project Structure

### Documentation (this feature)

```text
specs/001-julia-api/
├── plan.md              # This file (/speckit.plan command output)
├── research.md          # Phase 0 output (/speckit.plan command)
├── data-model.md        # Phase 1 output (/speckit.plan command)
├── quickstart.md        # Phase 1 output (/speckit.plan command)
├── contracts/           # Phase 1 output (/speckit.plan command)
│   ├── configuration-schema.toml  # Example config structure
│   └── api-reference.md           # Public API contracts
├── spec.md              # Feature specification
└── checklists/          # Validation checklists
    └── requirements.md
```

### Source Code (repository root)

```text
src/
├── TDSESolver.jl        # Main module (orchestrator)
├── grid/
│   ├── GPSGrid.jl       # GPS method with algebraic mapping
│   └── AngularGrid.jl   # Theta/phi grids for 3D
├── hamiltonian/
│   ├── Potential.jl     # Atomic potentials (H, He, Ar, Ne, Xe, custom)
│   └── FieldInteraction.jl  # Atom-field coupling (length gauge)
├── eigenstate/
│   ├── BoundStates.jl   # Ground state calculation
│   └── SMatrix.jl       # Time evolution operator in energy basis
├── field/
│   └── LaserPulse.jl    # Electric field, vector potential, pulse envelopes
├── propagator/
│   └── SplitOperator.jl # 3-step split-operator propagation
├── region_split/
│   └── RegionSplit.jl   # Inner/outer decomposition
├── continuum/
│   └── VolkovProjection.jl  # Momentum space projection
├── observable/
│   ├── MomentumDistribution.jl  # P(p), P(px,py), P(px,py,pz)
│   └── HHGSpectrum.jl   # Dipole moment, FFT, harmonic spectrum
├── io/
│   ├── ConfigParser.jl  # TOML configuration loading & validation
│   └── ResultsIO.jl     # HDF5 output with metadata
└── utils/
    ├── PhysicalUnits.jl # Atomic units conversions
    └── Validation.jl    # Parameter validation, benchmark comparison

tests/
├── benchmarks/          # Fortran reference outputs
│   ├── hydrogen_ground_state.h5
│   ├── helium_ionization_800nm.h5
│   └── helium_hhg_800nm.h5
├── integration/
│   ├── test_ionization_workflow.jl
│   └── test_hhg_workflow.jl
└── unit/
    ├── test_gps_grid.jl
    ├── test_split_operator.jl
    ├── test_volkov_projection.jl
    └── ... (one per module)

examples/
├── helium_800nm.toml    # Example configuration
├── hydrogen_scan.toml   # Parameter scan example
└── custom_potential.jl  # Advanced custom potential module

docs/
├── tutorial_notebooks/
│   ├── 01_basic_ionization.jl  # Pluto notebook
│   └── 02_hhg_calculation.jl
└── algorithm_references.md  # Fortran line numbers, paper equations
```

**Structure Decision**: Single Julia package structure selected. Scientific computing packages typically follow this pattern (e.g., DifferentialEquations.jl, Quantum.jl). Modular `src/` organization mirrors the 10-module architecture from constitution. `tests/benchmarks/` contains Fortran validation data. `examples/` provides TOML configuration templates. `docs/tutorial_notebooks/` uses Pluto.jl for interactive documentation.

## Complexity Tracking

> No constitution violations requiring justification.

All design choices align with established principles:
- Single package structure is standard for scientific Julia libraries
- 10-module architecture follows refactoring_plan.tex blueprint
- External dependencies (FFTW.jl, HDF5.jl, etc.) replace hardcoded Fortran implementations per Principle IV
- No additional abstraction layers introduced beyond necessary modularization
