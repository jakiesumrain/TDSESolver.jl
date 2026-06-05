# Benchmark Data

This directory contains reference output data from Fortran TDSE solver implementations for validation.

## Files

- `hydrogen_ground_state.h5` - Hydrogen ground state wavefunction and energy
- `helium_ionization_800nm.h5` - Helium momentum distribution at 800 nm, 5×10¹⁴ W/cm²
- `helium_hhg_800nm.h5` - Helium HHG spectrum at 800 nm

## Usage

These benchmarks are used in validation tests to ensure Julia implementation maintains numerical accuracy (correlation > 0.999) with the original Fortran codes.

## Source

Original Fortran programs:
- `explore/D_inner_out_volkov_3d_with_prob.f90` (ionization/ATI)
- `explore/rescatteing+hhg-he.f90` (HHG)

**Note**: Actual .h5 files need to be generated from Fortran runs and copied here before running validation tests.
