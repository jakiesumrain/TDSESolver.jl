# Tasks: Research-Grade TDSE Solver with User-Friendly Julia API

**Input**: Design documents from `/specs/001-julia-api/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/, quickstart.md

**Tests**: Tests are NOT explicitly requested in the feature specification. This implementation uses benchmark validation against Fortran reference outputs instead of traditional unit tests.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story. Each story can be validated independently against Fortran benchmarks.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and Julia package structure

- [X] T001 Create Julia package structure with src/, tests/, examples/, docs/ directories per plan.md
- [X] T002 Initialize Project.toml with package metadata and dependencies (LinearAlgebra, TOML, LegendrePolynomials, Bessels, FFTW, FastGaussQuadrature, HDF5, JLD2, SpecialFunctions)
- [X] T003 [P] Create .gitignore for Julia projects (Manifest.toml, .DS_Store, *.jl.cov, etc.)
- [X] T004 [P] Create README.md with installation instructions and quick start link
- [X] T005 [P] Setup tests/benchmarks/ directory structure and copy Fortran reference outputs (hydrogen_ground_state.h5, helium_ionization_800nm.h5, helium_hhg_800nm.h5)

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core utilities and infrastructure that ALL user stories depend on

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [X] T006 Implement PhysicalUnits module in src/utils/PhysicalUnits.jl with unit conversion functions (au_to_eV, intensity_SI_to_au, wavelength_nm_to_frequency_au) per contracts/api-reference.md
- [X] T007 [P] Implement ConfigParser module in src/io/ConfigParser.jl with load_config() function that parses TOML and validates parameters per contracts/configuration-schema.toml and FR-006
- [X] T008 [P] Implement Validation module in src/utils/Validation.jl with parameter validation functions (validate_physical_parameters, validate_numerical_stability) per FR-006
- [X] T009 Implement Potential module in src/hamiltonian/Potential.jl with get_potential() for built-in atoms (H, He, Ar, Ne, Xe) and parse_custom_potential() for Julia expressions per data-model.md PhysicalSystem entity
- [X] T010 [P] Implement GPSGrid module in src/grid/GPSGrid.jl with create_gps_grid() function using FastGaussQuadrature.jl per contracts/api-reference.md and research.md section 3
- [X] T011 [P] Implement AngularGrid module in src/grid/AngularGrid.jl with create_angular_grid() function per contracts/api-reference.md
- [X] T012 Implement ResultsIO module in src/io/ResultsIO.jl with save_results() and load_result() functions using HDF5.jl per data-model.md CalculationResults entity and research.md section 6

**Checkpoint**: Foundation ready - user story implementation can now begin in parallel

---

## Phase 3: User Story 1 - Basic Ionization Calculation (Priority: P1) 🎯 MVP

**Goal**: Researchers can compute photoelectron momentum distributions for helium atoms with minimal setup (< 10 minutes from installation to results)

**Independent Test**: Run helium ionization at 800 nm, 5×10¹⁴ W/cm², 10 fs pulse. Verify momentum distribution correlates > 0.999 with Fortran benchmark tests/benchmarks/helium_ionization_800nm.h5 per SC-003

### Implementation for User Story 1

- [ ] T013 [P] [US1] Implement BoundStates module in src/eigenstate/BoundStates.jl with compute_ground_state() function using imaginary time propagation or diagonalization per contracts/api-reference.md, following Fortran blueprint lines in D_inner_out_volkov_3d_with_prob.f90
- [ ] T014 [P] [US1] Implement SMatrix module in src/eigenstate/SMatrix.jl with compute_smatrix() function per contracts/api-reference.md formula from Fortran lines 255-293, consulting Tong & Chu 1997
- [ ] T015 [P] [US1] Implement LaserPulse module in src/field/LaserPulse.jl with create_laser_pulse() function including sin² envelope per contracts/api-reference.md
- [ ] T016 [US1] Implement SplitOperator module in src/propagator/SplitOperator.jl with propagate_split_operator!() function implementing 3-step method from Fortran lines 606-736, following algorithm from Tong & Chu 1997 Eq. 15
- [ ] T017 [US1] Implement RegionSplit module in src/region_split/RegionSplit.jl with smooth transition functions from PRA 74.031405(R)(2006) per Fortran blueprint
- [ ] T018 [US1] Implement VolkovProjection module in src/continuum/VolkovProjection.jl with project_volkov() function using Bessels.jl per research.md section 2, following Fortran lines 909-923 and PRA 74.031405 algorithm
- [ ] T019 [US1] Implement MomentumDistribution module in src/observable/MomentumDistribution.jl with compute_momentum_distribution() function computing P(p), P(px,py), P(px,pz), P(px,py,pz) per contracts/api-reference.md
- [ ] T020 [US1] Implement main TDSESolver orchestrator in src/TDSESolver.jl with run_tdse() function that coordinates all modules for ionization calculations per data-model.md state transitions
- [ ] T021 [US1] Create example configuration examples/helium_800nm.toml matching quickstart.md helium ionization example
- [ ] T022 [US1] Validate User Story 1: Run helium_800nm.toml configuration and verify correlation > 0.999 with tests/benchmarks/helium_ionization_800nm.h5 per FR-009 and SC-003

**Checkpoint**: At this point, User Story 1 should be fully functional - users can run basic ionization calculations and get momentum distributions matching Fortran benchmarks

---

## Phase 4: User Story 2 - High-Harmonic Generation Calculation (Priority: P2)

**Goal**: Researchers can compute HHG spectra for different atoms showing correct harmonic cutoff behavior

**Independent Test**: Run helium HHG at 800 nm, verify harmonic spectrum with cutoff position within 5% of Ip + 3.17Up per SC-009, odd-harmonic structure present

### Implementation for User Story 2

- [ ] T023 [P] [US2] Implement HHGSpectrum module in src/observable/HHGSpectrum.jl with compute_hhg_spectrum() function using FFTW.jl per contracts/api-reference.md and research.md section 4, following rescatteing+hhg-he.f90 blueprint
- [ ] T024 [US2] Extend TDSESolver.run_tdse() in src/TDSESolver.jl to support calculation_type = :hhg by computing dipole moment history during propagation per Fortran HHG blueprint
- [ ] T025 [US2] Create example configuration examples/helium_hhg.toml matching quickstart.md HHG example
- [ ] T026 [US2] Validate User Story 2: Run helium_hhg.toml and verify HHG cutoff position within 5% of theory, odd harmonics present, correlation with tests/benchmarks/helium_hhg_800nm.h5 per SC-009

**Checkpoint**: At this point, User Stories 1 AND 2 should both work independently - users can compute both ionization and HHG

---

## Phase 5: User Story 3 - Custom Potential Models (Priority: P3)

**Goal**: Advanced researchers can define custom atomic potentials beyond standard models

**Independent Test**: Define soft-core potential V(r) = -1/sqrt(r² + 0.5), verify ground state energy matches analytical value, ionization dynamics physically reasonable

### Implementation for User Story 3

- [ ] T027 [US3] Extend Potential.parse_custom_potential() in src/hamiltonian/Potential.jl to support Julia expression parsing with Meta.parse() and safe evaluation per spec.md design decision "Custom Potential Specification"
- [ ] T028 [US3] Add custom potential validation in src/hamiltonian/Potential.jl checking asymptotic behavior (r→∞ approaches 0 or -Z/r), boundedness per FR-013
- [ ] T029 [US3] Extend ConfigParser.load_config() in src/io/ConfigParser.jl to handle atom_type = :custom with potential expressions or module references per contracts/configuration-schema.toml
- [ ] T030 [US3] Create example configuration examples/soft_core.toml with custom potential expression matching quickstart.md custom potential example
- [ ] T031 [US3] Create example advanced module examples/custom_potential.jl demonstrating MyPotentials module approach per quickstart.md
- [ ] T032 [US3] Validate User Story 3: Run soft_core.toml and verify ground state energy, momentum distribution physically reasonable for soft-core potential

**Checkpoint**: All three user stories (ionization, HHG, custom potentials) should now be independently functional

---

## Phase 6: User Story 4 - Parameter Exploration and Optimization (Priority: P4)

**Goal**: Researchers can run systematic parameter scans with parallel execution and result aggregation

**Independent Test**: Define intensity scan 1e14 to 1e15 W/cm² in 10 steps, verify ionization probability increases with intensity, parallel execution works, results aggregated in single HDF5 file

### Implementation for User Story 4

- [ ] T033 [US4] Extend ConfigParser.load_config() in src/io/ConfigParser.jl to parse [scan] section with parameter ranges per spec.md design decision "Parameter Scan Specification"
- [ ] T034 [US4] Implement run_scan() function in src/TDSESolver.jl that distributes independent calculations across threads using Threads.@threads per contracts/api-reference.md and FR-012
- [ ] T035 [US4] Extend ResultsIO module in src/io/ResultsIO.jl to aggregate scan results into single HDF5 file with scan dimension as dataset axis per data-model.md
- [ ] T036 [US4] Add progress monitoring for scans in src/TDSESolver.jl showing completed/remaining calculations and ETA per FR-008
- [ ] T037 [US4] Create example configuration examples/intensity_scan.toml matching quickstart.md parameter scan example
- [ ] T038 [US4] Validate User Story 4: Run intensity_scan.toml with ENV["JULIA_NUM_THREADS"]="8", verify parallel speedup > 85% per SC-007, ionization trends physically correct

**Checkpoint**: Parameter scanning capability functional with parallel execution

---

## Phase 7: User Story 5 - Result Validation and Benchmarking (Priority: P5)

**Goal**: New users can verify solver accuracy before trusting for research by running standard benchmarks

**Independent Test**: Run benchmark suite comparing hydrogen ground state, helium ionization, helium HHG against reference data with automatic similarity metrics

### Implementation for User Story 5

- [ ] T039 [P] [US5] Implement ValidationReport generation in src/utils/Validation.jl with validate_calculation() function checking norm conservation, energy conservation, symmetry per contracts/api-reference.md and data-model.md ValidationReport entity
- [ ] T040 [P] [US5] Implement compare_with_fortran() function in src/utils/Validation.jl computing correlation between Julia and Fortran results per FR-009
- [ ] T041 [US5] Create benchmark test suite in tests/integration/test_benchmarks.jl that runs hydrogen ground state, helium ionization, helium HHG and validates against tests/benchmarks/ reference data per quickstart.md validation section
- [ ] T042 [US5] Integrate ValidationReport into CalculationResults in src/TDSESolver.jl so all calculations include validation diagnostics per FR-011 and data-model.md
- [ ] T043 [US5] Validate User Story 5: Run benchmark suite, verify all tests pass with correlations > 0.999 per SC-003, ground state energies within 1e-10 Ha per SC-006

**Checkpoint**: All five user stories complete - full TDSE solver with user-friendly API, validation, and benchmarking

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, examples, and final quality improvements

- [ ] T044 [P] Create tutorial notebook docs/tutorial_notebooks/01_basic_ionization.jl demonstrating User Story 1 workflow using Pluto.jl
- [ ] T045 [P] Create tutorial notebook docs/tutorial_notebooks/02_hhg_calculation.jl demonstrating User Story 2 workflow
- [ ] T046 [P] Create comprehensive docstrings for all public API functions following Julia documentation standards (? function_name should provide complete usage info) per spec.md design decision "Documentation Format"
- [ ] T047 [P] Create docs/algorithm_references.md documenting Fortran line numbers and paper equations for each algorithm per plan.md structure and FR-009
- [ ] T048 Run quickstart.md validation: Execute all four examples (helium ionization, HHG, parameter scan, custom potential) and verify outputs match expected results per SC-001 (< 10 minutes setup)
- [ ] T049 Performance optimization: Profile hot paths, apply @inbounds and @simd where safe, verify < 20% slower than Fortran per constitution Principle V and SC-002
- [ ] T050 [P] Create installation guide in README.md with MKL.jl optional optimization instructions per research.md performance notes
- [ ] T051 Final validation: Run full test suite (benchmarks, integration tests) and verify all success criteria SC-001 through SC-010 are met

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion (T001-T005) - BLOCKS all user stories
- **User Story 1 (Phase 3)**: Depends on Foundational (T006-T012) - MVP delivery
- **User Story 2 (Phase 4)**: Depends on Foundational (T006-T012) and partially on US1 (T020 orchestrator structure) - Can parallelize most tasks
- **User Story 3 (Phase 5)**: Depends on Foundational (T006-T012) and US1 (T020 orchestrator) - Independent of US2
- **User Story 4 (Phase 6)**: Depends on US1 complete (T013-T022) for base calculation capability - Independent of US2, US3
- **User Story 5 (Phase 7)**: Depends on US1, US2 complete (needs both ionization and HHG working) - Independent of US3, US4
- **Polish (Phase 8)**: Depends on all user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Core dependency - ALL other stories need this
  - Provides: run_tdse() orchestrator, ionization calculation capability, momentum distributions
  - Blocking for: US2 (needs orchestrator), US4 (needs base calculation), US5 (needs ionization benchmarks)

- **User Story 2 (P2)**: Independent extension of US1
  - Provides: HHG calculation capability, dipole moment computation, harmonic spectra
  - Depends on: US1 orchestrator structure (T020)
  - Can parallelize: Most HHG implementation (T023) is independent, only T024 modifies orchestrator

- **User Story 3 (P3)**: Independent extension of foundational Potential module
  - Provides: Custom potential capability
  - Depends on: Foundational Potential module (T009), US1 orchestrator (T020)
  - Independent of: US2, US4, US5

- **User Story 4 (P4)**: Builds on US1 for batch processing
  - Provides: Parameter scanning, parallel execution
  - Depends on: US1 complete (T013-T022) - needs working ionization calculations
  - Independent of: US2, US3, US5

- **User Story 5 (P5)**: Cross-cutting validation for US1 + US2
  - Provides: Benchmark validation, confidence reporting
  - Depends on: US1 complete (needs ionization), US2 complete (needs HHG)
  - Independent of: US3, US4

### Within Each User Story

**User Story 1** (Critical path: longest development time):
```
T006-T012 (Foundational) → T013, T014, T015 [P] → T016 → T017 → T018 → T019 → T020 → T021 → T022
```

**User Story 2** (Extends US1):
```
T006-T012 (Foundational), T020 (US1 orchestrator) → T023 [P], T024 → T025 → T026
```

**User Story 3** (Extends Foundational):
```
T009 (Foundational Potential), T020 (US1 orchestrator) → T027 → T028, T029 [P] → T030, T031 [P] → T032
```

**User Story 4** (Builds on US1):
```
T013-T022 (US1 complete) → T033, T034, T035, T036 [P] → T037 → T038
```

**User Story 5** (Requires US1 + US2):
```
T022 (US1 complete), T026 (US2 complete) → T039, T040 [P] → T041 → T042 → T043
```

### Parallel Opportunities

**Phase 1 (Setup)**: T003, T004, T005 can run in parallel after T001, T002

**Phase 2 (Foundational)**: T007, T008, T010, T011 can run in parallel after T006

**Phase 3 (User Story 1)**:
- T013, T014, T015 can run in parallel (different files: BoundStates.jl, SMatrix.jl, LaserPulse.jl)
- T021, T022 can overlap if T022 uses existing configs

**Phase 4 (User Story 2)**:
- T023 (HHGSpectrum.jl) can start immediately in parallel with T024 (orchestrator modification)

**Phase 5 (User Story 3)**:
- T028, T029 can run in parallel (validation vs. config parsing in different files)
- T030, T031 can run in parallel (two different example files)

**Phase 6 (User Story 4)**:
- T033, T034, T035, T036 involve different aspects but same orchestrator file - sequential

**Phase 7 (User Story 5)**:
- T039, T040 can run in parallel (ValidationReport generation vs. Fortran comparison - different functions)

**Phase 8 (Polish)**:
- T044, T045, T046, T047, T050 all work on different files - high parallelism

---

## Parallel Example: User Story 1 Initial Implementation

After Foundational phase completes, three core modules can be developed in parallel:

```bash
# Launch these three tasks together (different files, no dependencies):
Task T013: "Implement BoundStates module in src/eigenstate/BoundStates.jl"
Task T014: "Implement SMatrix module in src/eigenstate/SMatrix.jl"
Task T015: "Implement LaserPulse module in src/field/LaserPulse.jl"

# Once all three complete, proceed to:
Task T016: "Implement SplitOperator module" (needs SMatrix)
```

---

## Implementation Strategy

### MVP First (User Story 1 Only - Fastest Path to Value)

1. **Complete Phase 1: Setup** (T001-T005) - 1 day
2. **Complete Phase 2: Foundational** (T006-T012) - 3-4 days
   - ⚠️ CRITICAL BLOCKER - Nothing else can proceed until done
3. **Complete Phase 3: User Story 1** (T013-T022) - 10-12 days
   - This is the longest critical path
   - Requires careful algorithm translation from Fortran
   - Extensive validation against benchmarks
4. **STOP and VALIDATE**: Run helium_800nm.toml, verify correlation > 0.999 with Fortran
5. **Deploy/Demo MVP**: Users can now compute ionization momentum distributions

**MVP Timeline: ~14-17 days of focused development**

### Incremental Delivery (All User Stories)

1. **Phase 1 + Phase 2** → Foundation ready (4-5 days)
2. **+ Phase 3 (US1)** → Test independently → **MVP Release** (ionization calculations) (10-12 days)
3. **+ Phase 4 (US2)** → Test independently → **Release 2** (adds HHG capability) (3-4 days)
4. **+ Phase 5 (US3)** → Test independently → **Release 3** (adds custom potentials) (2-3 days)
5. **+ Phase 6 (US4)** → Test independently → **Release 4** (adds parameter scans) (3-4 days)
6. **+ Phase 7 (US5)** → Test independently → **Release 5** (adds validation suite) (2-3 days)
7. **+ Phase 8 (Polish)** → **Final Release** (documentation, optimization) (3-4 days)

**Total Timeline: ~27-35 days for complete implementation**

Each release adds value without breaking previous functionality. Users can start using ionization calculations after just ~17 days, while HHG, custom potentials, and scanning are added incrementally.

### Parallel Team Strategy (3 Developers)

With three developers and staggered starts:

**Week 1-2**: All three developers work together on Phase 1 + Phase 2 (Setup + Foundational)
- This is CRITICAL BLOCKING work
- Establishes shared architecture understanding
- Ensures consistent coding style

**Week 3-4**: Phase 3 (User Story 1) - All three developers collaborate
- Developer A: Grid + Eigenstate modules (T010, T011, T013, T014)
- Developer B: Propagation modules (T015, T016, T017)
- Developer C: Observable modules (T018, T019) + Orchestrator (T020)
- **Critical**: Frequent integration due to tight coupling

**Week 5**: After US1 MVP validated, split work:
- Developer A: Phase 4 (User Story 2 - HHG) (T023-T026)
- Developer B: Phase 5 (User Story 3 - Custom Potentials) (T027-T032)
- Developer C: Phase 6 (User Story 4 - Parameter Scans) (T033-T038)

**Week 6**: Convergence:
- Developer A: Phase 7 (User Story 5 - Validation) (T039-T043)
- Developer B: Phase 8 (Polish - Documentation) (T044-T047)
- Developer C: Phase 8 (Polish - Performance) (T048-T051)

**Parallel Team Timeline: ~6 weeks to complete feature**

---

## Notes

- **[P] tasks** indicate different files with no dependencies - safe for parallel execution
- **[Story] labels** (US1, US2, etc.) map tasks to user stories from spec.md for traceability
- Each user story should be independently completable and testable against Fortran benchmarks
- **Fortran Blueprint Fidelity**: T016, T017, T018 are CRITICAL - follow Fortran algorithms exactly per constitution Principle III
- **Scientific Paper Consultation**: T016 (split-operator) MUST reference Tong & Chu 1997, T018 (Volkov) MUST reference PRA 74.031405
- **Validation Gates**: T022, T026, T032, T038, T043 are validation checkpoints - MUST achieve correlation > 0.999 with Fortran per FR-009
- **Performance Critical**: T016 (propagation loop) and T018 (Volkov projection) are hot paths - apply optimizations carefully after correctness validated
- Stop at any checkpoint to validate story independently before proceeding
- Commit after each task or logical group of parallel tasks
- Document all algorithm sources (Fortran line numbers, paper equations) in code comments per T047
