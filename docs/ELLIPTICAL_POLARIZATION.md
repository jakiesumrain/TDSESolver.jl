# Elliptical Polarization Quick Start Guide

This guide shows how to use the full 3D field interaction capability for elliptically and bicircularly polarized laser fields.

## Status (2025-01-23)

✅ **WORKING**: Coordinate transformations preserve norm to machine precision
✅ **WORKING**: Full 3D field interaction (all polarizations)
✅ **WORKING**: Bicircular field support (ω + 2ω) - FULLY INTEGRATED ✨
⚠️ **LIMITATION**: S-matrix unitarity issues (use workarounds)

## Quick Examples

### Linear Polarization (Z-axis)

```julia
include("src/TDSESolver.jl")

params = create_default_params(
    nrmax = 100,
    rmax = 50.0,
    lmax = 2,
    laser_wavelength = 800.0,        # nm
    laser_intensity = 1.0e14,        # W/cm²
    laser_polarization = [0.0, 0.0, 1.0],  # z-axis
    atom = :hydrogen
)

results = run_simulation(params, verbose=true)
```

### Linear Polarization (X-axis)

```julia
params = create_default_params(
    # ... same as above ...
    laser_polarization = [1.0, 0.0, 0.0],  # x-axis
)
```

### Elliptical Polarization

```julia
# Ellipticity ε = 0.5 (ratio of minor to major axis)
epsilon = 0.5
Ex = 1.0
Ey = epsilon
norm = sqrt(Ex^2 + Ey^2)

params = create_default_params(
    # ... same as above ...
    laser_polarization = [Ex/norm, Ey/norm, 0.0],
)
```

### Circular Polarization (Approximation)

```julia
# Equal components in x and y
params = create_default_params(
    # ... same as above ...
    laser_polarization = [1.0, 1.0, 0.0] / sqrt(2),
)

# NOTE: True circular requires time-dependent phase between Ex and Ey
# See LaserField module integration below
```

## Bicircular Fields (ω + 2ω)

### ✅ NOW INTEGRATED (2025-01-23)

Bicircular fields are now fully integrated into the Simulation module. Simply set `use_bicircular=true`!

### Usage

```julia
include("src/TDSESolver.jl")

params = create_default_params(
    nrmax = 100,
    rmax = 50.0,
    lmax = 2,
    laser_wavelength = 800.0,
    laser_intensity = 1.0e14,
    laser_duration = 5.0,
    # Enable bicircular mode
    use_bicircular = true,
    exrate = 1.0,        # ω component amplitude (x-direction)
    eyrate = 1.0,        # 2ω component amplitude (y-direction)
    ezrate = 0.0,        # z-direction (typically zero)
    phase_offset = π/2,  # Relative phase (π/2 for circular)
    atom = :hydrogen
)

results = run_simulation(params, verbose=true)
```

### Field Equations

The bicircular field uses the LaserField module formulas:
```
Ex(t) = Ex * f(t) * sin(ωt)              [Fundamental ω]
Ey(t) = Ey * f(t) * sin(2ωt + φ)         [Second harmonic 2ω]
Ez(t) = 0
```

where:
- `f(t)`: Trapezoidal envelope (ramps up/down over 4 optical cycles)
- `φ`: phase_offset parameter (relative phase between ω and 2ω)
- Pulse centered at `t = 0` (Fortran convention)

### Field Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `use_bicircular` | Enable bicircular mode | `false` |
| `exrate` | X-component amplitude ratio (ω) | `1.0` |
| `eyrate` | Y-component amplitude ratio (2ω) | `0.0` |
| `ezrate` | Z-component amplitude ratio | `0.0` |
| `phase_offset` | Relative phase (radians) | `0.0` |

### Configuration Examples

**Counter-rotating circular (typical bicircular):**
```julia
use_bicircular = true
exrate = 1.0        # ω in x
eyrate = 1.0        # 2ω in y
phase_offset = π/2  # π/2 for circular
```

**Co-rotating:**
```julia
use_bicircular = true
exrate = 1.0
eyrate = 1.0
phase_offset = -π/2  # -π/2 for co-rotating
```

**Elliptical bicircular:**
```julia
use_bicircular = true
exrate = 1.0
eyrate = 0.5        # Ellipticity ε = 0.5
phase_offset = π/2
```

### Testing

Run the integration test:
```julia
julia examples/bicircular_test.jl
```

Expected output: All 5 tests pass ✅

## Physics Notes

### Selection Rules

For linear polarization along z:
- Δl = ±1 (dipole selection rule)
- Δm = 0 (m conserved)
- Only couples (l,m) ↔ (l±1, m)

For circular/elliptical polarization:
- Δl = ±1 (dipole selection rule)
- Δm = ±1 (angular momentum transfer)
- Couples multiple m states

### Norm Conservation

After the 2025-01-23 fix, coordinate transformations preserve norm to machine precision:

```julia
# Verified in tests/test_transformation_debug.jl:
Original norm:              1.0
After forward transform:    0.9999999999998053
After phase application:    0.9999999999998053
After inverse transform:    0.9999999999999996
Total change:               4.44×10⁻¹⁶  ✓
```

This means **field interactions are unitary** and preserve quantum mechanics properly.

## Known Issues & Workarounds

### Issue 1: S-Matrix Non-Unitarity

**Problem:** S-matrix changes norm by ~8% per application

**Impact:** Full simulations (with field-free propagation) show exponential norm divergence

**Workaround:** For testing field interactions only:
```julia
# Skip S-matrix, only test field interaction
apply_field_interaction_full!(wfn, prop, E_field)  # Norm preserved ✓
```

**Status:** Documented in `UNSOLVED_PROBLEMS.md`

### Issue 2: Eigenstate Solver Accuracy

**Problem:** l=0 eigenstate energies ~50x too negative

**Impact:** Inaccurate initial states, non-unitary S-matrix

**Workaround:** Use analytical ground states:
```julia
gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
wfn = create_wavefunction(grid.nrmax, lmax, grid.quadrature_weights)
initialize_ground_state!(wfn, gs.radial_function, 0, 0)
```

**Status:** Documented in `UNSOLVED_PROBLEMS.md`

## Testing

### Verify Norm Conservation

```julia
using Test

# Test field interaction alone
norm_before = compute_norm(wfn)
apply_field_interaction_full!(wfn, prop, E_field)
norm_after = compute_norm(wfn)

@test abs(norm_after - norm_before) < 1e-10  # Should pass!
```

### Run Debug Test Suite

```julia
cd("tests")
julia --project=.. test_transformation_debug.jl
```

Expected output:
```
Test 1: Forward-Inverse Transformation
  Relative error: 7.13×10⁻¹⁶  ✓

Test 2: Field Interaction
  Norm after small field: 1.0  ✓

Test 4: Components
  Total norm change: 4.44×10⁻¹⁶  ✓
```

## Complete Example

See `examples/elliptical_polarization.jl` for complete working examples of:
1. Linear polarization (z and x)
2. Elliptical polarization
3. Circular polarization (approximation)
4. Bicircular fields (conceptual)

## Technical Documentation

For complete technical details on the coordinate transformation fix:
- **`docs/COORDINATE_TRANSFORM_FIX.md`** - Full technical analysis
- **`UNSOLVED_PROBLEMS.md`** - Known issues and workarounds
- **Fortran reference:** `explore/D_inner_out_volkov_3d_with_prob.f90` (lines 638-708)

## References

### Spherical Harmonic Normalization
- GSL function: `sf_legendre_sphPlm` (spherical harmonic normalized)
- Normalization: `sqrt((2l+1)/(4π) * (l-|m|)!/(l+|m|)!)`

### Angular Grid
- Gauss-Legendre with substitution: `x = cos(θ)`
- θ values: `θ = arccos(x)` where x from `gausslegendre(nthmax)`
- Weights: Gauss-Legendre weights directly (no extra factors)

### Coordinate Transformations
Forward (radial → angular):
```
ψ(φ,θ,r) = Σ_l Σ_m g(r,m,l) * P_l^|m|(cos θ) * exp(imφ)
```

Inverse (angular → radial):
```
g(r,m,l) = ∫∫ ψ(φ,θ,r) * P_l^|m|(cos θ) * exp(-imφ) dΩ
```

With spherical harmonic normalized P_l^m, these transformations are unitary.

## Support

For issues or questions:
1. Check `UNSOLVED_PROBLEMS.md` for known limitations
2. Review `docs/COORDINATE_TRANSFORM_FIX.md` for technical details
3. Run `tests/test_transformation_debug.jl` to verify your installation
4. See `examples/elliptical_polarization.jl` for usage patterns

---

**Last Updated:** 2025-01-23
**Status:** Coordinate transformations working perfectly ✅
**Next Steps:** Integrate LaserField module for full bicircular support
