"""
Momentum Distribution Module

Computes momentum distributions from full momentum space wavefunction ψ(φ,θ,p).

**Algorithm follows Fortran D_inner_out_volkov_3d_with_prob.f90, lines 1646-1982:**

1. Transform ψ(φ,θ,p) → |ψ(px,py,pz)|² on 3D Cartesian grid
2. Integrate Cartesian grid for 2D distributions
3. Extract slices at specific planes

**Theory:**

Starting from momentum space wavefunction ψ(φ,θ,p) on spherical grid, transform to
Cartesian representation |ψ(px,py,pz)|², then compute observables:

1. **3D Cartesian Grid:**
   - |ψ(px,py,pz)|² on uniform Cartesian grid
   - Enables all downstream analysis

2. **2D Integrated Distributions (Fortran lines 1944-1949):**
   - P(px,py) = ∫ |ψ(px,py,pz)|² dpz
   - P(px,pz) = ∫ |ψ(px,py,pz)|² dpy
   - P(py,pz) = ∫ |ψ(px,py,pz)|² dpx

3. **2D Slices (Fortran lines 1973-1982):**
   - |ψ(px,py,pz=z₀)|² at specific pz
   - |ψ(px,py=y₀,pz)|² at specific py
   - |ψ(px=x₀,py,pz)|² at specific px

**Coordinate Transformations:**

Spherical (p,θ,φ) ↔ Cartesian (px,py,pz):
- px = p sin(θ) cos(φ)
- py = p sin(θ) sin(φ)
- pz = p cos(θ)

**Usage:**

```julia
using .MomentumDistribution

# Step 1: Transform to 3D Cartesian grid
rate_xyz = transform_to_cartesian_grid(psai_p, momentum_grid, p_max, n_grid)

# Step 2: Compute 2D integrated distributions (like Fortran pxyz_2D-momentum_probe.txt)
pxy_rate, pxz_rate, pyz_rate = compute_integrated_2d_distributions(rate_xyz, p_max, n_grid)

# Step 3: Extract 2D slices (like Fortran plane_xyz_momentum_probe.txt)
slice_xy_at_z0 = extract_2d_slice(rate_xyz, "xy", 0.0, p_max, n_grid)
slice_xz_at_y0 = extract_2d_slice(rate_xyz, "xz", 0.0, p_max, n_grid)
slice_yz_at_x0 = extract_2d_slice(rate_xyz, "yz", 0.0, p_max, n_grid)
```

**References:**
- Fortran: D_inner_out_volkov_3d_with_prob.f90
  - Lines 1646-1933: Cartesian grid transformation (subroutine probe_cartersian)
  - Lines 1944-1949: 2D integrated distributions
  - Lines 1973-1982: 2D slices at origin planes
"""
module MomentumDistribution

export transform_to_cartesian_grid
export compute_integrated_2d_distributions
export extract_2d_slice
export get_axis_grids

using LinearAlgebra
using Printf

"""
    spherical_to_cartesian(p, theta, phi) -> (px, py, pz)

Convert spherical momentum coordinates to Cartesian.

**Formula:**
- px = p sin(θ) cos(φ)
- py = p sin(θ) sin(φ)
- pz = p cos(θ)
"""
function spherical_to_cartesian(p::Float64, theta::Float64, phi::Float64)
    px = p * sin(theta) * cos(phi)
    py = p * sin(theta) * sin(phi)
    pz = p * cos(theta)
    return (px, py, pz)
end

"""
    cartesian_to_spherical(px, py, pz) -> (p, theta, phi)

Convert Cartesian momentum coordinates to spherical.

**Formula:**
- p = √(px² + py² + pz²)
- θ = arccos(pz/p)
- φ = arctan(py/px)
"""
function cartesian_to_spherical(px::Float64, py::Float64, pz::Float64)
    p = sqrt(px^2 + py^2 + pz^2)

    if p < 1e-10
        return (0.0, 0.0, 0.0)
    end

    theta = acos(clamp(pz/p, -1.0, 1.0))
    phi = atan(py, px)

    return (p, theta, phi)
end

"""
    transform_to_cartesian_grid(psai_p, momentum_grid, p_max, n_grid) -> rate_xyz

Transform ψ(φ,θ,p) on spherical grid to |ψ(px,py,pz)|² on Cartesian grid.

**Algorithm (following Fortran lines 1646-1933):**

1. Create 3D Cartesian grid: (px, py, pz) ∈ [-p_max, p_max]³
2. For each Cartesian grid point:
   - Convert (px, py, pz) → (p, θ, φ)
   - Interpolate ψ(φ, θ, p) from spherical grid
   - Store |ψ|² at Cartesian point
3. Apply spherical cutoff: only keep points with p < p_max

**Parameters:**
- `psai_p`: Momentum space wavefunction [iphi, itheta, ip] (from VolkovProjection)
- `momentum_grid`: MomentumGrid structure with p, theta, phi grids
- `p_max`: Maximum momentum for Cartesian grid
- `n_grid`: Number of grid points per axis (Fortran uses nxmax, nymax, nzmax)

**Returns:**
- `rate_xyz`: 3D array [ix, iy, iz] with |ψ(px,py,pz)|²
- Grid corresponds to px, py, pz ∈ [-p_max, p_max] with spacing pgrid = 2*p_max/n_grid

**Fortran equivalent:**
```fortran
rate_xyz(nx,ny,nz) = |zpsai|²  ! Lines 1916-1924
```
"""
function transform_to_cartesian_grid(
    psai_p::Array{ComplexF64,3},
    momentum_grid,
    p_max::Float64,
    n_grid::Int
)
    # Cartesian grid: -n_grid:n_grid indices
    rate_xyz = zeros(Float64, 2*n_grid+1, 2*n_grid+1, 2*n_grid+1)

    # Grid spacing
    pgrid = p_max / n_grid

    # Extract spherical grids
    p_sph = momentum_grid.p
    theta_sph = momentum_grid.theta
    phi_sph = momentum_grid.phi

    println("Transforming to Cartesian grid...")
    println("  Grid: $(2*n_grid+1)³ points")
    println("  p_max: $p_max a.u.")
    println("  Grid spacing: $(@sprintf("%.4f", pgrid)) a.u.")

    # Loop over Cartesian grid
    for iz in 1:(2*n_grid+1)
        nz = iz - n_grid - 1  # -n_grid:n_grid
        pz = nz * pgrid

        for iy in 1:(2*n_grid+1)
            ny = iy - n_grid - 1
            py = ny * pgrid

            for ix in 1:(2*n_grid+1)
                nx = ix - n_grid - 1
                px = nx * pgrid

                # Convert to spherical
                p, theta, phi = cartesian_to_spherical(px, py, pz)

                # Apply spherical cutoff (like Fortran line 1735, 1941)
                if p > p_max || p < 1e-10
                    rate_xyz[ix, iy, iz] = 0.0
                    continue
                end

                # Interpolate ψ from spherical grid (nearest neighbor for now)
                # TODO: Could use linear interpolation for better accuracy
                ip = argmin(abs.(p_sph .- p))
                itheta = argmin(abs.(theta_sph .- theta))

                # Handle φ periodicity
                phi_wrapped = mod(phi, 2π)
                iphi = argmin(abs.(phi_sph .- phi_wrapped))

                # Get wavefunction value
                psi_val = psai_p[iphi, itheta, ip]

                # Store probability density
                rate_xyz[ix, iy, iz] = abs2(psi_val)
            end
        end

        # Progress indicator
        if nz % 10 == 0
            percent = (iz - 1) / (2*n_grid+1) * 100
            println("  Progress: $(@sprintf("%.1f", percent))% (nz = $nz)")
        end
    end

    println("✓ Cartesian grid transformation complete")
    return rate_xyz
end

"""
    compute_integrated_2d_distributions(rate_xyz, p_max, n_grid) -> (pxy_rate, pxz_rate, pyz_rate)

Compute 2D integrated momentum distributions by summing over perpendicular direction.

**Algorithm (following Fortran lines 1938-1952):**

For xy-plane: P(px,py) = Σ_pz |ψ(px,py,pz)|² × pgrid
For xz-plane: P(px,pz) = Σ_py |ψ(px,py,pz)|² × pgrid
For yz-plane: P(py,pz) = Σ_px |ψ(px,py,pz)|² × pgrid

**Parameters:**
- `rate_xyz`: 3D Cartesian grid from `transform_to_cartesian_grid()`
- `p_max`: Maximum momentum
- `n_grid`: Number of grid points per axis

**Returns:**
- `pxy_rate`: 2D array [ix, iy] for P(px, py)
- `pxz_rate`: 2D array [ix, iz] for P(px, pz)
- `pyz_rate`: 2D array [iy, iz] for P(py, pz)

**Fortran equivalent:**
```fortran
pxyrate(nx,ny) = pxyrate(nx,ny) + rate*pgrid  ! Line 1944
pxzrate(nx,nz) = pxzrate(nx,nz) + rate*pgrid  ! Line 1945
pyzrate(ny,nz) = pyzrate(ny,nz) + rate*pgrid  ! Line 1946
```

**Output format matches Fortran file pxyz_2D-momentum_probe.txt**
"""
function compute_integrated_2d_distributions(
    rate_xyz::Array{Float64,3},
    p_max::Float64,
    n_grid::Int
)
    pgrid = p_max / n_grid

    # Initialize 2D distributions
    pxy_rate = zeros(Float64, 2*n_grid+1, 2*n_grid+1)
    pxz_rate = zeros(Float64, 2*n_grid+1, 2*n_grid+1)
    pyz_rate = zeros(Float64, 2*n_grid+1, 2*n_grid+1)

    println("Computing 2D integrated distributions...")

    # Integrate over all three dimensions
    for iz in 1:(2*n_grid+1)
        for iy in 1:(2*n_grid+1)
            for ix in 1:(2*n_grid+1)
                rate = rate_xyz[ix, iy, iz]

                # Integrate for each 2D plane
                pxy_rate[ix, iy] += rate * pgrid  # Sum over pz
                pxz_rate[ix, iz] += rate * pgrid  # Sum over py
                pyz_rate[iy, iz] += rate * pgrid  # Sum over px
            end
        end
    end

    println("✓ 2D integrated distributions complete")
    println("  P(px,py) sum: $(@sprintf("%.6e", sum(pxy_rate)))")
    println("  P(px,pz) sum: $(@sprintf("%.6e", sum(pxz_rate)))")
    println("  P(py,pz) sum: $(@sprintf("%.6e", sum(pyz_rate)))")

    return (pxy_rate, pxz_rate, pyz_rate)
end

"""
    extract_2d_slice(rate_xyz, plane, offset, p_max, n_grid) -> slice_2d

Extract 2D slice from 3D Cartesian grid at specified offset.

**Algorithm (following Fortran lines 1973-1982):**

Extract |ψ(px,py,pz=z₀)|² or similar slices at specific coordinate values.

**Parameters:**
- `rate_xyz`: 3D Cartesian grid from `transform_to_cartesian_grid()`
- `plane`: "xy", "xz", or "yz"
- `offset`: Value of perpendicular coordinate (e.g., pz=0.0 for xy-plane)
- `p_max`: Maximum momentum
- `n_grid`: Number of grid points per axis

**Returns:**
- `slice_2d`: 2D array with probability density at the slice

**Examples:**
```julia
# Origin plane slices (like Fortran plane_xyz_momentum_probe.txt)
slice_xy = extract_2d_slice(rate_xyz, "xy", 0.0, p_max, n_grid)  # pz = 0
slice_xz = extract_2d_slice(rate_xyz, "xz", 0.0, p_max, n_grid)  # py = 0
slice_yz = extract_2d_slice(rate_xyz, "yz", 0.0, p_max, n_grid)  # px = 0

# Off-origin slices
slice_xy_pz1 = extract_2d_slice(rate_xyz, "xy", 1.0, p_max, n_grid)  # pz = 1 a.u.
```

**Fortran equivalent:**
```fortran
rate_xyz(i,j,0)    ! xy-plane at pz=0, Line 1979
rate_xyz(0,j,i)    ! yz-plane at px=0, Line 1979
rate_xyz(j,0,i)    ! xz-plane at py=0, Line 1979
```
"""
function extract_2d_slice(
    rate_xyz::Array{Float64,3},
    plane::String,
    offset::Float64,
    p_max::Float64,
    n_grid::Int
)
    pgrid = p_max / n_grid

    # Find index corresponding to offset
    offset_index = round(Int, offset / pgrid) + n_grid + 1

    # Check bounds
    if offset_index < 1 || offset_index > 2*n_grid+1
        error("Offset $offset is outside grid range [-$p_max, $p_max]")
    end

    # Extract slice based on plane
    if plane == "xy"
        # pz = offset
        slice_2d = rate_xyz[:, :, offset_index]
        println("Extracted xy-plane slice at pz = $(@sprintf("%.3f", offset)) a.u.")
    elseif plane == "xz"
        # py = offset
        slice_2d = rate_xyz[:, offset_index, :]
        println("Extracted xz-plane slice at py = $(@sprintf("%.3f", offset)) a.u.")
    elseif plane == "yz"
        # px = offset
        slice_2d = rate_xyz[offset_index, :, :]
        println("Extracted yz-plane slice at px = $(@sprintf("%.3f", offset)) a.u.")
    else
        error("Invalid plane: $plane. Must be 'xy', 'xz', or 'yz'")
    end

    return slice_2d
end

"""
    get_axis_grids(p_max, n_grid) -> (px_axis, py_axis, pz_axis)

Get Cartesian axis grids matching the 3D rate_xyz array indices.

**Returns:**
- Tuple of (px_axis, py_axis, pz_axis) vectors with momentum values

**Usage:**
```julia
px_axis, py_axis, pz_axis = get_axis_grids(p_max, n_grid)
# Then px_axis[ix] gives the px value for index ix in rate_xyz or pxy_rate
```
"""
function get_axis_grids(p_max::Float64, n_grid::Int)
    pgrid = p_max / n_grid

    # Create axis vectors
    px_axis = [(i - n_grid - 1) * pgrid for i in 1:(2*n_grid+1)]
    py_axis = [(i - n_grid - 1) * pgrid for i in 1:(2*n_grid+1)]
    pz_axis = [(i - n_grid - 1) * pgrid for i in 1:(2*n_grid+1)]

    return (px_axis, py_axis, pz_axis)
end

end  # module MomentumDistribution
