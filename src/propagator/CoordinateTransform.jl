"""
    CoordinateTransform

Module for transforming wavefunction between radial and angular coordinate representations.

Implements transformations matching Fortran D_inner_out_volkov_3d_with_prob.f90:
- Lines 638-655: Radial → Angular (g → psai_sp)
- Lines 688-708: Angular → Radial (psai_sp → g)

# Physics
The wavefunction can be represented in two ways:

1. Radial representation: g(r, m, l)
   - Decomposed into angular momentum channels
   - Good for field-free evolution (S-matrix diagonal in l)

2. Angular coordinate representation: ψ(φ, θ, r)
   - Wavefunction in real space angular coordinates
   - Good for applying dipole interaction V = -r⃗·E⃗

# Transformation Formulas
Forward (radial → angular):
    ψ(φ, θ, r) = Σ_l Σ_m g(r, m, l) Y_l^m(θ, φ)
               = Σ_m [Σ_l g(r, m, l) P_l^|m|(cos θ)] exp(imφ)

Inverse (angular → radial):
    g(r, m, l) = ∫∫ ψ(φ, θ, r) Y_l^m*(θ, φ) sin(θ) dθ dφ
               ≈ Σ_θ Σ_φ ψ(φ, θ, r) P_l^|m|(cos θ) exp(-imφ) w_θ w_φ

# Fortran Reference
D_inner_out_volkov_3d_with_prob.f90, lines 635-708
"""
module CoordinateTransform

using LinearAlgebra
using Printf

# Required types from other modules
using ..Wavefunction: WavefunctionData
using ..AngularGrid: AngularGridData
using ..GPSGrid: GPSGridData

export transform_radial_to_angular!, transform_angular_to_radial!

"""
    transform_radial_to_angular!(psai_sp::Array{ComplexF64,3},
                                  wfn::WavefunctionData,
                                  grid::AngularGridData)

Transform wavefunction from radial to angular representation.

# Arguments
- `psai_sp::Array{ComplexF64,3}`: Output array [nfi, nth, nr] (modified in-place)
- `wfn::WavefunctionData`: Input wavefunction in radial representation g(nr, m, l)
- `grid::AngularGridData`: Angular grid with precomputed transformation arrays

# Algorithm (Fortran lines 638-655)
```fortran
do nr=1,nrmax
  do nth=1,nthmax_sp
    Teemp(m) = Σ_l g(r,m,l) * P_l^|m|(cos θ)
    do nfi=1,nfimax_sp
      psai_sp(nfi,nth,nr) = Σ_m Teemp(m) * exp(imφ)
    enddo
  enddo
enddo
```

This computes: ψ(φ,θ,r) = Σ_l Σ_m g(r,m,l) P_l^|m|(cos θ) exp(imφ)

# Performance
- Complexity: O(nrmax × nθ × nφ × lmax²)
- Can parallelize over nr (OpenMP in Fortran)
"""
function transform_radial_to_angular!(psai_sp::Array{ComplexF64,3},
                                      wfn::WavefunctionData,
                                      grid::AngularGridData)
    nrmax = wfn.nrmax
    lmax = wfn.lmax
    nthmax = grid.nthmax
    nfimax = grid.nphimax

    # Initialize output
    fill!(psai_sp, 0.0 + 0.0im)

    # Match Fortran indexing: g_rever_lm(l,m,nr) vs g[nr, m_idx, l_idx]
    # We'll access g[nr, m+lmax+1, l+1] directly

    # Loop over radial points (can be parallelized)
    for nr in 1:nrmax
        # Loop over theta points
        for nth in 1:nthmax
            # Temporary array for m summation (Fortran: Teemp(m))
            Teemp = zeros(ComplexF64, 2*lmax+1)

            # Sum over l for each m (Fortran lines 641-646)
            for m in -lmax:lmax
                i = abs(m)  # Minimum l for given |m|

                for l in i:lmax
                    l_idx = l + 1
                    m_idx = m + lmax + 1

                    # Fortran line 644: Teemp(m) += g_rever_lm(l,m,nr) * plgd_sp_lm(l,m,nth)
                    Teemp[m+lmax+1] += wfn.g[nr, m_idx, l_idx] * grid.plgd_sp[l_idx, m_idx, nth]
                end
            end

            # Sum over m with phi phase factors (Fortran lines 648-652)
            for nfi in 1:nfimax
                for m in -lmax:lmax
                    m_idx = m + lmax + 1

                    # Fortran line 650: psai_sp(nfi,nth,nr) += Teemp(m) * expfi_sp(m,nfi)
                    psai_sp[nfi, nth, nr] += Teemp[m_idx] * grid.expfi_sp[m_idx, nfi]
                end
            end
        end
    end
end

"""
    transform_angular_to_radial!(wfn::WavefunctionData,
                                  psai_sp::Array{ComplexF64,3},
                                  grid::AngularGridData,
                                  gps_grid::GPSGridData)

Transform wavefunction from angular to radial representation.

# Arguments
- `wfn::WavefunctionData`: Output wavefunction in radial representation (modified in-place)
- `psai_sp::Array{ComplexF64,3}`: Input array [nfi, nth, nr]
- `grid::AngularGridData`: Angular grid with precomputed transformation arrays
- `gps_grid::GPSGridData`: GPS grid for quadrature weights

# Algorithm (Fortran lines 688-708)
```fortran
do nr=1,nrmax
  ! First sum over phi
  TTemp(nth,m) = Σ_φ psai_sp(nfi,nth,nr) * exp(-imφ) * w_φ

  ! Then sum over theta
  g(m,l,nr) = Σ_θ TTemp(nth,m) * P_l^|m|(cos θ) * w_θ
enddo
```

This computes: g(r,m,l) = ∫∫ ψ(φ,θ,r) P_l^|m|(cos θ) exp(-imφ) dθ dφ

# Performance
- Complexity: O(nrmax × nθ × nφ × lmax²)
- Can parallelize over nr (OpenMP in Fortran)
"""
function transform_angular_to_radial!(wfn::WavefunctionData,
                                      psai_sp::Array{ComplexF64,3},
                                      grid::AngularGridData,
                                      gps_grid::GPSGridData)
    nrmax = wfn.nrmax
    lmax = wfn.lmax
    nthmax = grid.nthmax
    nfimax = grid.nphimax

    # Temporary storage (Fortran: g_split_rever(m,l,nr))
    g_split_rever = zeros(ComplexF64, 2*lmax+1, lmax+1, nrmax)

    # Loop over radial points (can be parallelized)
    for nr in 1:nrmax
        # Temporary array for phi summation (Fortran: TTemp(nth,m))
        TTemp = zeros(ComplexF64, nthmax, 2*lmax+1)

        # Sum over phi for each (theta, m) (Fortran lines 691-697)
        for m in -lmax:lmax
            m_idx = m + lmax + 1

            for nth in 1:nthmax
                for nfi in 1:nfimax
                    # Fortran line 694: TTemp(nth,m) += psai_sp(nfi,nth,nr) * expfi_sp_conj(nfi,m) * wfi_sp(nfi)
                    TTemp[nth, m_idx] += psai_sp[nfi, nth, nr] *
                                          grid.expfi_sp_conj[nfi, m_idx] *
                                          grid.phi_weights[nfi]
                end
            end
        end

        # Sum over theta for each (l, m) (Fortran lines 699-705)
        for l in 0:lmax
            l_idx = l + 1

            for m in -l:l
                m_idx = m + lmax + 1

                for nth in 1:nthmax
                    # Fortran line 702: g_split_rever(m,l,nr) += TTemp(nth,m) * plgd_sp_rever(nth,m,l) * wth_sp(nth)
                    # Note: plgd_sp_rever is same as plgd_sp (symmetric)
                    g_split_rever[m_idx, l_idx, nr] += TTemp[nth, m_idx] *
                                                        grid.plgd_sp[l_idx, m_idx, nth] *
                                                        grid.theta_weights[nth]
                end
            end
        end
    end

    # Copy back to wfn.g (Fortran lines 711-714: reverse array to reduce time consumption)
    for l in 0:lmax
        l_idx = l + 1

        for m in -l:l
            m_idx = m + lmax + 1

            for nr in 1:nrmax
                wfn.g[nr, m_idx, l_idx] = g_split_rever[m_idx, l_idx, nr]
            end
        end
    end
end

"""
    print_transform_info(wfn::WavefunctionData, grid::AngularGrid Data)

Print diagnostic information about coordinate transformations.
"""
function print_transform_info(wfn::WavefunctionData, grid::AngularGridData)
    println("Coordinate Transformation Information:")
    println("  Radial representation: g[$(wfn.nrmax), $(2*wfn.lmax+1), $(wfn.lmax+1)]")
    println("  Angular representation: ψ[$(grid.nphimax), $(grid.nthmax), $(wfn.nrmax)]")

    radial_size = sizeof(wfn.g) / 1024^2
    angular_size = sizeof(ComplexF64) * grid.nphimax * grid.nthmax * wfn.nrmax / 1024^2

    @printf("  Radial array memory: %.2f MB\n", radial_size)
    @printf("  Angular array memory: %.2f MB\n", angular_size)
    @printf("  Total memory: %.2f MB\n", radial_size + angular_size)

    # Estimate computational cost
    nops_forward = wfn.nrmax * grid.nthmax * (wfn.lmax^2 + grid.nphimax * wfn.lmax)
    nops_inverse = wfn.nrmax * (grid.nthmax * grid.nphimax * wfn.lmax + grid.nthmax * wfn.lmax^2)

    @printf("  Forward transform ops: %.2e\n", Float64(nops_forward))
    @printf("  Inverse transform ops: %.2e\n", Float64(nops_inverse))
end

end  # module CoordinateTransform
