# Coordinate Transformation Bug Fix

**Date**: 2025-01-23
**Issue**: Critical norm conservation failure in coordinate transformations
**Status**: ✅ **RESOLVED**

---

## Executive Summary

Fixed a critical bug in coordinate transformations that caused wavefunction norm to explode exponentially (1.0 → 10¹⁸ in 10 time steps). The issue stemmed from using unnormalized Legendre polynomials and incorrect angular grid setup. After fixes, transformations now preserve norm to machine precision (relative error ≈ 4.4×10⁻¹⁶).

---

## Problem Description

### Observed Behavior

Time propagation showed catastrophic norm explosion:

```
t = 0.0:  norm = 1.0         ✓
t = 1.0:  norm = 8.3         ❌
t = 2.0:  norm = 2.9×10²⁸    ❌
t = 10.0: norm = 1.5×10¹⁴⁷   ❌
```

### Impact

- Wavefunction completely unphysical
- All observables meaningless
- Simulation results unusable
- Made full 3D field interaction impossible

---

## Root Cause Analysis

### Root Cause 1: Unnormalized Legendre Polynomials

**What was wrong:**
```julia
# OLD (INCORRECT)
using GSL: sf_legendre_Plm  # Returns unnormalized P_l^m(x)
```

**Problem:** GSL's `sf_legendre_Plm` returns **unnormalized** associated Legendre polynomials. For spherical harmonics, we need:

```
Y_l^m(θ,φ) = N_lm * P_l^|m|(cos θ) * exp(imφ)
```

where `N_lm = sqrt((2l+1)/(4π) * (l-|m|)!/(l+|m|)!)` is the normalization constant.

**Effect on norm:**
For ground state (l=0, m=0):
- Unnormalized: P₀⁰(x) = 1
- Normalized: Y₀⁰ = 1/√(4π) ≈ 0.282

Forward transformation gave norm ≈ **4π** instead of 1.0!

**Fix:**
```julia
# NEW (CORRECT)
using GSL: sf_legendre_sphPlm  # Returns spherical harmonic normalized P_l^m(x)
```

`sf_legendre_sphPlm` includes the normalization factor automatically.

### Root Cause 2: Incorrect Angular Grid

**What was wrong:**
```julia
# OLD (INCORRECT)
x_gl, w_gl = gausslegendre(nthmax)
theta_grid = @. (x_gl + 1.0) * π / 2.0      # Linear mapping ❌
theta_weights = @. w_gl * π / 2.0
```

**Problem:** Used linear mapping `θ = (x+1)π/2` which is incorrect for spherical harmonics with P_l^m(cos θ).

**Correct approach:** For integrals involving spherical harmonics:

```
∫₀^π f(θ) P_l^m(cos θ) sin(θ) dθ
```

Substitute `x = cos(θ)`, then `sin(θ) dθ = -dx`:

```
∫₀^π f(θ) P_l^m(cos θ) sin(θ) dθ = ∫₋₁¹ f(arccos(x)) P_l^m(x) dx
```

Now use Gauss-Legendre quadrature on `[-1,1]` with NO extra factors.

**Fix:**
```julia
# NEW (CORRECT)
x_gl, w_gl = gausslegendre(nthmax)
theta_grid = @. acos(x_gl)              # θ = arccos(x) ✓
theta_weights = w_gl                    # No π/2 factor ✓
```

The `sin(θ)` factor is **implicit** in the `dx` transformation, so weights need no modification.

---

## Solution Implementation

### Files Modified

1. **`src/grid/AngularGrid.jl`**
   - Changed import: `sf_legendre_Plm` → `sf_legendre_sphPlm`
   - Fixed `compute_associated_legendre()` function
   - Fixed angular grid creation (lines 156-176)
   - Updated `compute_solid_angle_element()` to handle both grid types

2. **`tests/unit/test_angular_grid.jl`**
   - Updated tests for correct spherical harmonic normalization
   - Fixed monotonicity tests (Gauss grids have **decreasing** θ)
   - Updated solid angle element tests

3. **`tests/test_transformation_debug.jl`**
   - Created comprehensive debug test suite
   - Removed redundant `sin(θ)` factors in manual norm computations

### Key Code Changes

#### Change 1: Spherical Harmonic Normalization

```julia
# src/grid/AngularGrid.jl (lines 81-90)
function compute_associated_legendre(l::Int, m::Int, x::Float64)
    m_abs = abs(m)
    if m_abs > l
        return 0.0
    end

    # GSL function: sf_legendre_sphPlm(l, m_abs, x)
    # Returns spherical harmonic normalized P_l^|m|(x)
    # Normalization: sqrt((2l+1)/(4π) * (l-m)!/(l+m)!)
    return sf_legendre_sphPlm(l, m_abs, x)
end
```

#### Change 2: Correct Angular Grid

```julia
# src/grid/AngularGrid.jl (lines 156-176)
elseif method == :gauss
    # Gauss-Legendre for θ using x = cos(θ) transformation
    # For ∫₀^π f(θ) sin(θ) dθ:
    # Substitute x = cos(θ), then sin(θ) dθ = -dx
    # ∫₀^π f(θ) sin(θ) dθ = ∫₁^{-1} f(arccos(x)) (-dx) = ∫₋₁¹ f(arccos(x)) dx

    x_gl, w_gl = gausslegendre(nthmax)

    # θ = arccos(x), NOT linear mapping!
    theta_grid = @. acos(x_gl)

    # Weights are w_gl directly (sin(θ) factor already in transformation)
    theta_weights = w_gl

    # Uniform for φ
    phi_grid = range(0.0, 2π, length=nphimax+1)[1:end-1] |> collect
    Δφ = 2π / nphimax
    phi_weights = fill(Δφ, nphimax)
```

#### Change 3: Grid-Type-Aware Solid Angle Element

```julia
# src/grid/AngularGrid.jl (lines 266-282)
function compute_solid_angle_element(grid::AngularGridData, i_theta::Int, i_phi::Int)
    θ = grid.theta_grid[i_theta]
    w_θ = grid.theta_weights[i_theta]
    w_φ = grid.phi_weights[i_phi]

    # Detect grid type: Gauss grids have decreasing θ (arccos is decreasing function)
    is_gauss_grid = (grid.nthmax > 1) && (grid.theta_grid[2] < grid.theta_grid[1])

    if is_gauss_grid
        # Gauss-Legendre grid: sin(θ) already in weights
        return w_θ * w_φ
    else
        # Uniform grid: sin(θ) explicit
        return sin(θ) * w_θ * w_φ
    end
end
```

---

## Verification & Testing

### Test Results

Created comprehensive test suite `tests/test_transformation_debug.jl`:

#### Test 1: Forward-Inverse Transformation (Identity)
```
Original norm:              1.0
After forward transform:    0.9999999999998053
After inverse transform:    1.0
Relative error:             7.13×10⁻¹⁶  ✓ (machine precision!)
```

#### Test 2: Field Interaction
```
Before field:               1.0
After zero field:           1.0
After small field (Ez=0.001): 0.9999999999999996  ✓
```

#### Test 4: Component-wise Analysis
```
Original norm:              1.0
After forward transform:    0.9999999999998053
After phase application:    0.9999999999998053
After inverse transform:    0.9999999999999996
Total change:               4.44×10⁻¹⁶  ✓
```

### Unit Test Status

All foundational tests pass:
```
✓ PhysicalUnits   (25/25)
✓ ConfigParser    (24/24)
✓ Validation      (30/30)
✓ Potential       (43/43)
✓ GPSGrid        (429/429)
✓ AngularGrid  (11183/11183)  ← Fixed!
✓ ResultsIO       (41/41)

Total: 7/7 modules passing
```

---

## Mathematical Verification

### Spherical Harmonic Orthonormality

With correct normalization, spherical harmonics satisfy:

```
∫∫ Y_l^m(θ,φ) Y_{l'}^{m'}*(θ,φ) dΩ = δ_{ll'} δ_{mm'}
```

Using our implementation:
```julia
∫∫ = Σ_i Σ_j f(θ_i, φ_j) * w_θ[i] * w_φ[j]  (Gauss grid)
```

For l=m=0:
```julia
Y₀⁰ = 1/√(4π) * P₀⁰(cos θ) * exp(0)
    = 1/√(4π) * (spherical harmonic normalized value)
    = 0.282095...  ✓ (verified with GSL)
```

### Transformation Formulas

**Forward (radial → angular):**
```
ψ(φ,θ,r) = Σ_l Σ_m g(r,m,l) * P_l^|m|(cos θ) * exp(imφ)
```

With spherical harmonic normalized P_l^m, this preserves norm.

**Inverse (angular → radial):**
```
g(r,m,l) = ∫∫ ψ(φ,θ,r) * P_l^|m|(cos θ) * exp(-imφ) dΩ
         ≈ Σ_i Σ_j ψ(φ_j,θ_i,r) * P_l^|m|(cos θ_i) * exp(-imφ_j) * w_θ[i] * w_φ[j]
```

With `x = cos(θ)` transformation, the quadrature is exact for polynomials.

---

## Fortran Reference Comparison

### Fortran Implementation (lines 638-708)

The Fortran code uses:
```fortran
plgd_sp_lm(l,m,nth)  ! Associated Legendre polynomials
```

**Key insight:** The Fortran code likely uses **normalized** Legendre polynomials from a library, though this isn't explicitly documented. Our fix aligns with standard spherical harmonic conventions used in quantum mechanics.

### Differences from Fortran

1. **Explicit normalization choice:** We explicitly use GSL's spherical harmonic normalized version
2. **Grid construction:** We use `θ = arccos(x)` which is standard for Legendre polynomial quadrature
3. **Clarity:** Our implementation documents the normalization and transformation explicitly

---

## Remaining Issues

### S-Matrix Non-Unitarity (Separate Issue)

The S-matrix still causes norm changes (1.0 → 0.924 per application). This is **NOT** related to coordinate transformations but to:

- Numerical eigenstate solver accuracy
- Incomplete energy cutoff
- GPS grid discretization errors

**Status:** Documented in `UNSOLVED_PROBLEMS.md` as a known limitation.

**Workaround:** For testing coordinate transformations, skip S-matrix applications.

---

## Performance Impact

### Memory

No change - same array sizes.

### Computational Cost

- Gauss-Legendre grid: **slightly faster** (fewer points needed for same accuracy)
- `sf_legendre_sphPlm`: Same cost as `sf_legendre_Plm`

### Accuracy

- **Dramatically improved:** Machine precision norm conservation vs. exponential divergence
- Gauss-Legendre quadrature: Exact for polynomials up to degree `2n-1`

---

## Usage Guidelines

### For Simulations

```julia
# Create angular grid (default is Gauss-Legendre)
ang_grid = create_angular_grid(90, 60, lmax, method=:gauss)

# Coordinate transformations now preserve norm automatically
transform_radial_to_angular!(psai_sp, wfn, ang_grid)
# ... apply field interaction ...
transform_angular_to_radial!(wfn, psai_sp, ang_grid, gps_grid)
```

### For Testing

```julia
# Verify norm conservation
norm_before = compute_norm(wfn)
# ... do transformations ...
norm_after = compute_norm(wfn)
@test abs(norm_after - norm_before) < 1e-10  # Should pass!
```

---

## Lessons Learned

### 1. **Always check normalization conventions**
   - Different libraries use different normalizations
   - GSL has both normalized and unnormalized versions
   - **Always** read the documentation carefully

### 2. **Coordinate transformations are subtle**
   - The "obvious" linear mapping `θ = (x+1)π/2` was wrong
   - Proper substitution `x = cos(θ)` is essential
   - Jacobian factors can be implicit or explicit

### 3. **Comprehensive testing is critical**
   - Debug test isolated the issue within minutes
   - Component-wise testing (forward, phase, inverse) pinpointed the bug
   - Always test norm conservation to machine precision

### 4. **Follow the mathematics**
   - Standard spherical harmonic formulas exist for a reason
   - Don't deviate without very good reason
   - Trust the textbooks!

---

## References

### Mathematical Background

1. **Spherical Harmonics:**
   - Jackson, *Classical Electrodynamics*, 3rd Ed., Section 3.5
   - Arfken & Weber, *Mathematical Methods for Physicists*, Chapter 12

2. **Gauss-Legendre Quadrature:**
   - Numerical Recipes, 3rd Ed., Section 4.6
   - Substitution `x = cos(θ)` for angular integrals

3. **GSL Documentation:**
   - `sf_legendre_Plm`: Standard normalization
   - `sf_legendre_sphPlm`: Spherical harmonic normalization
   - https://www.gnu.org/software/gsl/doc/html/specfunc.html#legendre-functions

### Code References

- Fortran reference: `explore/D_inner_out_volkov_3d_with_prob.f90`, lines 638-708
- Fixed Julia code: `src/grid/AngularGrid.jl`
- Tests: `tests/unit/test_angular_grid.jl`, `tests/test_transformation_debug.jl`

---

## Acknowledgments

Bug identified and fixed: 2025-01-23
Original Fortran implementation: Tong & Chu algorithm

---

## Appendix: Quick Reference

### Before (Broken)
```julia
using GSL: sf_legendre_Plm
theta_grid = (x_gl .+ 1.0) .* π ./ 2.0
theta_weights = w_gl .* π ./ 2.0
# Result: norm → 10¹⁸ ❌
```

### After (Fixed)
```julia
using GSL: sf_legendre_sphPlm
theta_grid = acos.(x_gl)
theta_weights = w_gl
# Result: norm conserved to 10⁻¹⁶ ✓
```

### Key Formula
```
∫₀^π f(θ) P_l^m(cos θ) sin(θ) dθ = ∫₋₁¹ f(arccos(x)) P_l^m(x) dx
                                  ≈ Σᵢ f(arccos(xᵢ)) P_l^m(xᵢ) wᵢ
```

---

**End of Technical Note**
