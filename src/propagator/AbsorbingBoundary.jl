"""
    AbsorbingBoundary

Module for absorbing boundary conditions in TDSE propagation.

Implements cos^(1/4) absorbing mask following Tong's method to prevent
wavefunction reflection from grid boundaries. Essential for accurate
HHG (High Harmonic Generation) calculations.

# Physics
Without absorbing boundaries, outgoing wavepackets reflect from the grid
boundary and contaminate the inner region, producing spurious features in:
- High harmonic spectra (HHG)
- Dipole acceleration calculations
- Recollision physics

# Fortran Reference
From rescatteing+hhg-he.f90, Lines 766-776:
```fortran
!!!!!!!!!!!!added 20140311 absorber 根据tong's paper 用于处理边界的反射
do l=0,lmax
    do nr=1,nrmax
        if(r(nr).le.r0) then
            g(nr,l) = g(nr,l)
        else
            g(nr,l) = g(nr,l) * dcos(pai*(r(nr)-r0)/2.d0/(rmax-r0))**0.25d0
        endif
    enddo
enddo
```

# Usage
```julia
# Create absorber
absorber = create_absorbing_boundary(grid, r0=100.0)

# Apply after each propagation step
apply_absorbing_boundary!(wfn, absorber)
```
"""
module AbsorbingBoundary

using ..GPSGrid: GPSGridData
using ..Wavefunction: WavefunctionData

export AbsorbingBoundaryData, create_absorbing_boundary, apply_absorbing_boundary!
export compute_absorbed_norm

"""
    AbsorbingBoundaryData

Container for absorbing boundary mask.

# Fields
- `grid::GPSGridData`: Reference to radial grid
- `r0::Float64`: Absorber start radius (a.u.)
- `rmax::Float64`: Maximum grid radius (a.u.)
- `mask::Vector{Float64}`: Pre-computed absorbing mask [nrmax]

# Mask Function
For r > r0:
    mask(r) = cos^(1/4)(π(r - r0) / (2(rmax - r0)))

For r ≤ r0:
    mask(r) = 1.0

The mask smoothly transitions from 1 at r0 to 0 at rmax.
"""
struct AbsorbingBoundaryData
    grid::GPSGridData
    r0::Float64
    rmax::Float64
    mask::Vector{Float64}
end

"""
    create_absorbing_boundary(grid::GPSGridData;
                              r0::Float64=0.0,
                              r0_fraction::Float64=0.8) -> AbsorbingBoundaryData

Create absorbing boundary mask.

# Arguments
- `grid::GPSGridData`: Radial grid
- `r0::Float64=0.0`: Absorber start radius. If 0, uses r0_fraction * rmax
- `r0_fraction::Float64=0.8`: Fraction of rmax to start absorbing (default: 80%)

# Returns
- `AbsorbingBoundaryData`: Absorber with pre-computed mask

# Example
```julia
# Start absorbing at 80% of grid
absorber = create_absorbing_boundary(grid)

# Start absorbing at specific radius
absorber = create_absorbing_boundary(grid, r0=100.0)

# Start absorbing at 70% of grid
absorber = create_absorbing_boundary(grid, r0_fraction=0.7)
```

# Physics Notes
- r0 should be large enough that bound states are unaffected
- Typically r0 > 50 a.u. for hydrogen-like atoms
- For HHG, r0 should be beyond the quiver radius: r0 > E0/ω²
"""
function create_absorbing_boundary(grid::GPSGridData;
                                   r0::Float64=0.0,
                                   r0_fraction::Float64=0.8)
    rmax = maximum(grid.radial_grid)

    # Determine r0
    if r0 <= 0.0
        r0 = r0_fraction * rmax
    end

    # Validate r0
    if r0 >= rmax
        error("r0 ($r0) must be less than rmax ($rmax)")
    end

    # Pre-compute mask
    nrmax = grid.nrmax
    mask = ones(Float64, nrmax)

    for ir in 1:nrmax
        r = grid.radial_grid[ir]

        if r > r0
            # Fortran: dcos(pai*(r(nr)-r0)/2.d0/(rmax-r0))**0.25d0
            arg = π * (r - r0) / (2.0 * (rmax - r0))
            mask[ir] = cos(arg)^0.25
        end
        # else mask[ir] = 1.0 (already initialized)
    end

    @info "Created absorbing boundary" r0=r0 rmax=rmax

    return AbsorbingBoundaryData(grid, r0, rmax, mask)
end

"""
    apply_absorbing_boundary!(wfn::WavefunctionData, absorber::AbsorbingBoundaryData)

Apply absorbing boundary mask to wavefunction.

# Arguments
- `wfn::WavefunctionData`: Wavefunction (modified in-place)
- `absorber::AbsorbingBoundaryData`: Absorber with pre-computed mask

# Algorithm
For each (l, m) channel and radial point:
    g(r, m, l) → g(r, m, l) × mask(r)

# When to Apply
Call this function AFTER each propagation step, typically after the
second S-matrix application in the split-operator scheme:

```julia
# Full propagation step with absorbing boundary
apply_s_matrix!(wfn, prop)
apply_field_interaction_full!(wfn, prop, E_field)
apply_s_matrix!(wfn, prop)
apply_absorbing_boundary!(wfn, absorber)  # Apply here
```

# Performance
Very fast: O(nrmax × lmax × (2lmax+1)) simple multiplications.
"""
function apply_absorbing_boundary!(wfn::WavefunctionData, absorber::AbsorbingBoundaryData)
    nrmax = wfn.nrmax
    lmax = wfn.lmax
    mask = absorber.mask

    # Apply mask to all (l, m) channels
    for l in 0:lmax
        l_index = l + 1
        for m in -l:l
            m_index = m + lmax + 1

            for ir in 1:nrmax
                wfn.g[ir, m_index, l_index] *= mask[ir]
            end
        end
    end
end

"""
    compute_absorbed_norm(wfn::WavefunctionData, absorber::AbsorbingBoundaryData) -> Float64

Compute the norm of wavefunction in the absorbing region (r > r0).

# Arguments
- `wfn::WavefunctionData`: Wavefunction
- `absorber::AbsorbingBoundaryData`: Absorber

# Returns
- `Float64`: Norm in absorbing region

# Usage
Track how much population is being absorbed over time.
Large absorbed population indicates strong ionization.
"""
function compute_absorbed_norm(wfn::WavefunctionData, absorber::AbsorbingBoundaryData)
    nrmax = wfn.nrmax
    lmax = wfn.lmax
    r0 = absorber.r0
    grid = absorber.grid

    absorbed_norm = 0.0

    for l in 0:lmax
        l_index = l + 1
        for m in -l:l
            m_index = m + lmax + 1

            for ir in 1:nrmax
                if grid.radial_grid[ir] > r0
                    absorbed_norm += abs2(wfn.g[ir, m_index, l_index]) * wfn.grid_weights[ir]
                end
            end
        end
    end

    return absorbed_norm
end

"""
    print_absorber_info(absorber::AbsorbingBoundaryData)

Print diagnostic information about absorber.
"""
function print_absorber_info(absorber::AbsorbingBoundaryData)
    println("Absorbing Boundary Information:")
    println("  Start radius r0: $(absorber.r0) a.u.")
    println("  Grid rmax: $(absorber.rmax) a.u.")
    println("  Absorbing width: $(absorber.rmax - absorber.r0) a.u.")

    # Find where mask drops below certain thresholds
    mask = absorber.mask
    grid = absorber.grid

    idx_90 = findfirst(m -> m < 0.9, mask)
    idx_50 = findfirst(m -> m < 0.5, mask)
    idx_10 = findfirst(m -> m < 0.1, mask)

    if idx_90 !== nothing
        println("  Mask = 0.9 at r = $(round(grid.radial_grid[idx_90], digits=2)) a.u.")
    end
    if idx_50 !== nothing
        println("  Mask = 0.5 at r = $(round(grid.radial_grid[idx_50], digits=2)) a.u.")
    end
    if idx_10 !== nothing
        println("  Mask = 0.1 at r = $(round(grid.radial_grid[idx_10], digits=2)) a.u.")
    end
end

end  # module AbsorbingBoundary
