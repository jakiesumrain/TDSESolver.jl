"""
Region Splitting Module

Implements smooth wavefunction decomposition into inner (bound) and outer (continuum) regions
for photoelectron momentum distribution calculations.

Based on algorithm from PRA 74, 031405(R) (2006)

**Theory:**

Split wavefunction into inner and outer regions:
    ψ(r,t) = ψ_inner(r,t) + ψ_outer(r,t)

Using smooth transition function:
    f_split(r) = 1 / (1 + exp(-(r - R_c)/Δ))

Where:
- R_c: Critical radius (typically 80-120 a.u.)
- Δ: Smoothness parameter (typically 3-10 a.u.)

**Properties:**
- f_split(r << R_c) ≈ 0  (inner region)
- f_split(r >> R_c) ≈ 1  (outer region)
- Smooth transition over ~4Δ width

**Splitting:**
    ψ_inner(r) = [1 - f_split(r)] · ψ(r)
    ψ_outer(r) = f_split(r) · ψ(r)

**Usage:**

```julia
using .RegionSplit

# Create splitter with R_c=100 a.u., Δ=5 a.u.
splitter = create_region_splitter(radial_grid, R_c=100.0, delta=5.0)

# Split wavefunction (modifies wfn_outer in-place)
split_wavefunction!(wfn_outer, wfn_total, splitter)

# Check populations
pop_inner, pop_outer = compute_split_populations(wfn_total, wfn_outer, radial_grid)
```

**References:**
- PRA 74, 031405(R) (2006): Region splitting method
- Fortran: D_inner_out_volkov_3d_with_prob.f90, lines 467-474, 786-794
"""
module RegionSplit

export RegionSplitterData
export create_region_splitter
export split_wavefunction!
export compute_split_populations
export get_split_function

using LinearAlgebra

"""
    RegionSplitterData

Data structure for region splitting.

**Fields:**
- `R_c::Float64`: Critical radius (transition center)
- `delta::Float64`: Smoothness parameter (transition width ~4Δ)
- `split_function::Vector{Float64}`: f_split(r) values on radial grid
- `nrmax::Int`: Number of radial grid points
"""
struct RegionSplitterData
    R_c::Float64
    delta::Float64
    split_function::Vector{Float64}
    nrmax::Int
end

"""
    create_region_splitter(radial_grid, R_c=100.0, delta=5.0) -> RegionSplitterData

Create region splitter with smooth transition function.

**Algorithm (PRA 74, 031405(R) Eq. 5):**

```
f_split(r) = 1 / (1 + exp(-(r - R_c)/Δ))
```

**Parameters:**
- `radial_grid`: GPSGridData with radial grid points
- `R_c`: Critical radius (default 100.0 a.u.)
  - Should be in outer region but before absorber
  - Typical: 0.5 × rmax to 0.7 × rmax
- `delta`: Smoothness parameter (default 5.0 a.u.)
  - Larger Δ → smoother transition, wider boundary
  - Smaller Δ → sharper transition, more localized
  - Typical: 3-10 a.u.

**Returns:**
- `RegionSplitterData`: Splitter with precomputed f_split(r)

**Example:**
```julia
# Standard splitting for rmax=150 a.u.
splitter = create_region_splitter(radial_grid, R_c=100.0, delta=5.0)

# Sharp splitting for precise localization
splitter = create_region_splitter(radial_grid, R_c=80.0, delta=3.0)
```

**Fortran Reference:**
- D_inner_out_volkov_3d_with_prob.f90, lines 469-474:
```fortran
r_split = 100d0
del_r_split = 5d0
do i=1,nrmax
    split_zone(i) = 1.d0/(1.d0 + dexp(-(r(i)-r_split)/del_r_split))
enddo
```
"""
function create_region_splitter(radial_grid, R_c::Float64=100.0, delta::Float64=5.0)
    nrmax = radial_grid.nrmax
    r = radial_grid.radial_grid

    # Compute smooth transition function f_split(r)
    # f_split(r) = 1 / (1 + exp(-(r - R_c)/delta))
    split_function = zeros(Float64, nrmax)

    for i in 1:nrmax
        # Smooth sigmoid function
        argument = -(r[i] - R_c) / delta

        # Avoid overflow for large negative arguments
        if argument < -50.0
            split_function[i] = 1.0  # exp(-50) ≈ 2e-22 ≈ 0
        elseif argument > 50.0
            split_function[i] = 0.0  # exp(50) → ∞, 1/(1+∞) → 0
        else
            split_function[i] = 1.0 / (1.0 + exp(argument))
        end
    end

    @info "Created region splitter" R_c delta nrmax transition_width=4*delta

    return RegionSplitterData(R_c, delta, split_function, nrmax)
end

"""
    split_wavefunction!(wfn_outer, wfn_total, splitter)

Split total wavefunction into outer region component.

**Algorithm:**

```
ψ_outer(r,m,l) = f_split(r) · ψ_total(r,m,l)
```

Implicitly:
```
ψ_inner(r,m,l) = [1 - f_split(r)] · ψ_total(r,m,l)
            = ψ_total(r,m,l) - ψ_outer(r,m,l)
```

**Parameters:**
- `wfn_outer`: WavefunctionData to store outer region (modified in-place)
- `wfn_total`: WavefunctionData with full wavefunction (unchanged)
- `splitter`: RegionSplitterData with f_split(r)

**Modifies:**
- `wfn_outer.g[ir, m_idx, l+1]` ← f_split(r) × wfn_total.g[ir, m_idx, l+1]

**Example:**
```julia
# Initialize outer wavefunction (zero)
wfn_outer = Wavefunction.create_wavefunction(nrmax, lmax, weights)

# Extract outer region
split_wavefunction!(wfn_outer, wfn_total, splitter)

# Inner region (if needed): wfn_inner = wfn_total - wfn_outer
```

**Fortran Reference:**
- D_inner_out_volkov_3d_with_prob.f90, lines 788-794:
```fortran
do l=0,lmax
    do m=-l,l
        do i=1,nrmax
            psai_p_temp(i,m,l) = g(i,m,l) * split_zone(i)  ! Outer
            g(i,m,l) = g(i,m,l) - psai_p_temp(i,m,l)        ! Inner
        enddo
    enddo
enddo
```
"""
function split_wavefunction!(wfn_outer, wfn_total, splitter::RegionSplitterData)
    nrmax = splitter.nrmax
    lmax = wfn_total.lmax

    # Extract split function
    f_split = splitter.split_function

    # Apply splitting: ψ_outer(r) = f_split(r) · ψ_total(r)
    for l in 0:lmax
        for m in -l:l
            m_idx = m + lmax + 1
            for ir in 1:nrmax
                wfn_outer.g[ir, m_idx, l+1] = f_split[ir] * wfn_total.g[ir, m_idx, l+1]
            end
        end
    end

    return nothing
end

"""
Compute populations in inner and outer regions.

Returns tuple (pop_inner, pop_outer)
"""
function compute_split_populations(wfn_total, wfn_outer, radial_grid)
    # Total norm
    norm_total_sq = 0.0
    norm_outer_sq = 0.0

    lmax = wfn_total.lmax
    nrmax = radial_grid.nrmax
    weights = radial_grid.quadrature_weights

    for l in 0:lmax
        for m in -l:l
            m_idx = m + lmax + 1
            for ir in 1:nrmax
                # Total
                norm_total_sq += abs2(wfn_total.g[ir, m_idx, l+1]) * weights[ir]

                # Outer
                norm_outer_sq += abs2(wfn_outer.g[ir, m_idx, l+1]) * weights[ir]
            end
        end
    end

    pop_total = norm_total_sq
    pop_outer = norm_outer_sq
    pop_inner = pop_total - pop_outer

    return (pop_inner, pop_outer)
end

"""
    get_split_function(splitter) -> Vector{Float64}

Get the splitting function f_split(r) values.

**Returns:**
- `Vector{Float64}`: f_split(r_i) for i=1:nrmax

**Example:**
```julia
f_split = get_split_function(splitter)

# Plot transition region
using Plots
plot(radial_grid.radial_grid, f_split,
     xlabel="r (a.u.)", ylabel="f_split(r)",
     title="Region Splitting Function")
```
"""
function get_split_function(splitter::RegionSplitterData)
    return splitter.split_function
end

end  # module RegionSplit
