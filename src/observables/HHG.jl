"""
High Harmonic Generation (HHG) Module

Calculates HHG spectrum from TDSE wavefunction using two methods:

**Method 1: Eigenstate Expansion** (Active in Fortran lines 891-899)
- Projects wavefunction onto eigenstate basis
- Computes dipole transitions between states
- Requires eigenvalues and eigenvectors

**Method 2: Direct Coordinate Space** (Commented in Fortran lines 905-909)
- Direct calculation from wavefunction
- Acceleration or length gauge
- No eigenstate projection needed

**References:**
- Fortran: rescatteing+hhg-he.f90, lines 891-899 (Method 1), 905-909 (Method 2)
- FFT subroutines: lines 1027-1356

**Usage:**
```julia
# Enable HHG in simulation
params = SimulationParams(
    enable_hhg = true,
    hhg_method = :eigenstate,  # or :acceleration, :length
)

# After simulation
spectrum = compute_hhg_spectrum(results.observables.dipole_moment, dt, ω₀)
```
"""
module HHG

export compute_dipole_eigenstate!
export compute_dipole_acceleration!
export compute_dipole_length!
export compute_hhg_spectrum
export save_dipole_to_file
export save_hhg_spectrum_to_file

using LinearAlgebra
using FFTW
using Printf

"""
    project_to_eigenstates(wfn, ham, grid) -> C_nl

Project wavefunction onto eigenstate basis to get time-dependent coefficients.

**Algorithm:**
    Cₙₗ = ⟨φₙₗ|ψ⟩ = Σᵣ φₙₗ*(r) g(r,m=0,l) w(r)

**Parameters:**
- `wfn`: WavefunctionData with g(r,m,l)
- `ham`: HamiltonianData with eigenvectors
- `grid`: GPSGridData with quadrature weights

**Returns:**
- `C_nl`: Matrix [n, l+1] of eigenstate coefficients
"""
function project_to_eigenstates(wfn, ham, grid)
    n_max = size(ham.eigenvectors, 2)
    lmax = ham.lmax
    nrmax = grid.nrmax

    C_nl = zeros(ComplexF64, n_max, lmax+1)

    # Project onto each eigenstate
    for l in 0:lmax
        for n in 1:n_max
            # Only project m=0 component (z-axis symmetry)
            m_idx = wfn.lmax + 1  # m=0

            for ir in 1:nrmax
                φ_nl = ham.eigenvectors[ir, n, l+1]
                g_val = wfn.g[ir, m_idx, l+1]
                w = grid.quadrature_weights[ir]

                C_nl[n, l+1] += conj(φ_nl) * g_val * w
            end
        end
    end

    return C_nl
end

"""
    compute_dipole_eigenstate!(obs, wfn, ham, grid, t)

Compute dipole moment using eigenstate expansion (Method 1).

**Algorithm (Fortran lines 891-899):**
```fortran
Dipole_1S = (0,0)
do n=1,nrmax
    do nr=1,nrmax
        Dipole_1S = f_r(nr) * conj(vec(nr,1,0)) * exp(i*E_n*t) *
                    vec(nr,n,1) * Cp(n,1) * (1/(2l+1)) * coef2(nr) + Dipole_1S
    enddo
enddo
```

**Formula:**
    d(t) = Σₙ ⟨1s|r|np⟩ Cₙ,ₗ₌₁(t) exp(i(Eₙₗ - E₁ₛ)t)

Where:
- `Cₙₗ(t)`: Time-dependent eigenstate coefficients from projection
- `⟨1s|r|np⟩`: Dipole matrix element (length gauge)
- Phase accounts for energy difference between states

**Parameters:**
- `obs`: ObservablesData (modified in-place)
- `wfn`: WavefunctionData
- `ham`: HamiltonianData with eigenstates
- `grid`: GPSGridData
- `t`: Current time
"""
function compute_dipole_eigenstate!(obs, wfn, ham, grid, t)
    if obs.dipole_moment === nothing
        return  # HHG not enabled
    end

    # Project wavefunction onto eigenstates
    C_nl = project_to_eigenstates(wfn, ham, grid)

    nrmax = grid.nrmax
    n_max = size(ham.eigenvectors, 2)
    r = grid.radial_grid
    weights = grid.quadrature_weights

    # Ground state: n=1, l=0
    φ_1s = ham.eigenvectors[:, 1, 1]  # l+1 = 1
    E_1s = ham.eigenvalues[1, 1]

    dipole = 0.0 + 0.0im

    # Sum over excited states: n, l=1 (p-states)
    for n in 1:n_max
        # Dipole matrix element: ⟨1s|r|np⟩
        dipole_matrix = 0.0 + 0.0im

        φ_np = ham.eigenvectors[:, n, 2]  # l=1, index = l+1 = 2
        E_np = ham.eigenvalues[n, 2]

        for ir in 1:nrmax
            # Length gauge: ⟨1s|r|np⟩
            dipole_matrix += conj(φ_1s[ir]) * r[ir] * φ_np[ir] * weights[ir]
        end

        # Time-dependent coefficient for (n, l=1) state
        C_np = C_nl[n, 2]

        # Phase factor: exp(i(Eₙₗ - E₁ₛ)t)
        phase = exp(im * (E_np - E_1s) * t)

        # Angular momentum factor: 1/(2l+1) = 1/3 for l=1
        angular_factor = 1.0 / 3.0

        # Accumulate dipole
        dipole += dipole_matrix * C_np * phase * angular_factor
    end

    # Additional angular factor from Y₁₀
    dipole *= sqrt(3.0 / (4π))

    push!(obs.dipole_moment, dipole)
end

"""
    compute_dipole_acceleration!(obs, wfn, pot, grid, E_field)

Compute dipole moment using acceleration form (Method 2a).

**Algorithm (Fortran lines 905-909, commented):**
```fortran
! dipolea(mtt) = (-cos(theta)/r² - E(t)) * |ψ(theta,r)|² * weights
```

**Formula:**
    d_acc(t) = ⟨ψ|-∂V/∂z + E_z(t)|ψ⟩

For z-polarized field and l=0 ↔ l=1, m=0 transition:
    d_acc = ∫ g₁₀*(r) × (-Z/r² + E_z) × g₀₀(r) × w(r) dr × √(3/4π)

**Parameters:**
- `obs`: ObservablesData (modified in-place)
- `wfn`: WavefunctionData
- `pot`: PotentialFunction (for nuclear charge Z)
- `grid`: GPSGridData
- `E_field`: Electric field [Ex, Ey, Ez] at current time
"""
function compute_dipole_acceleration!(obs, wfn, pot, grid, E_field)
    if obs.dipole_moment === nothing
        return  # HHG not enabled
    end

    nrmax = grid.nrmax
    r = grid.radial_grid
    weights = grid.quadrature_weights

    # Get nuclear charge
    Z = if pot.atom_type == :hydrogen
        1.0
    elseif pot.atom_type == :helium
        2.0
    else
        1.0  # Default
    end

    # Extract l=0 and l=1, m=0 components
    m_idx = wfn.lmax + 1  # m=0
    g_l0 = wfn.g[:, m_idx, 1]  # l=0
    g_l1 = wfn.g[:, m_idx, 2]  # l=1

    dipole = 0.0 + 0.0im

    # Radial integral: ⟨l=1| -∂V/∂z + E_z |l=0⟩
    for ir in 1:nrmax
        # Force term: -∂V/∂z = -Z cos(θ)/r²
        # For z-axis (θ=0), cos(θ)=1
        # Plus external field E_z
        if r[ir] > 1e-10
            force = Z / r[ir]^2 + E_field[3]
        else
            force = E_field[3]
        end

        # Dipole matrix element
        dipole += conj(g_l1[ir]) * force * g_l0[ir] * weights[ir]
    end

    # Angular factor from Y₁₀
    dipole *= sqrt(3.0 / (4π))

    push!(obs.dipole_moment, dipole)
end

"""
    compute_dipole_length!(obs, wfn, grid)

Compute dipole moment using length form (Method 2b).

**Algorithm (Fortran line 906, commented):**
```fortran
! dipolel(mtt) = (-cos(theta)/r²) * |ψ(theta,r)|² * weights
```

**Formula:**
    d_len(t) = ⟨ψ|r cos(θ)|ψ⟩ = ⟨ψ|z|ψ⟩

For l=0 ↔ l=1, m=0 transition:
    d_len = ∫ g₁₀*(r) × r × g₀₀(r) × w(r) dr × √(3/4π)

**Parameters:**
- `obs`: ObservablesData (modified in-place)
- `wfn`: WavefunctionData
- `grid`: GPSGridData
"""
function compute_dipole_length!(obs, wfn, grid)
    if obs.dipole_moment === nothing
        return  # HHG not enabled
    end

    nrmax = grid.nrmax
    r = grid.radial_grid
    weights = grid.quadrature_weights

    # Extract l=0 and l=1, m=0 components
    m_idx = wfn.lmax + 1  # m=0
    g_l0 = wfn.g[:, m_idx, 1]  # l=0
    g_l1 = wfn.g[:, m_idx, 2]  # l=1

    dipole = 0.0 + 0.0im

    # Radial integral: ⟨l=1| r |l=0⟩
    for ir in 1:nrmax
        dipole += conj(g_l1[ir]) * r[ir] * g_l0[ir] * weights[ir]
    end

    # Angular factor from Y₁₀
    dipole *= sqrt(3.0 / (4π))

    push!(obs.dipole_moment, dipole)
end

"""
    hann_window(n::Int) -> Vector{Float64}

Create Hann window function for FFT.

**Formula:**
    w[i] = 0.5 × (1 - cos(2π×i/(n-1)))

Used to reduce spectral leakage in FFT.
"""
function hann_window(n::Int)
    return [0.5 * (1 - cos(2π * i / (n-1))) for i in 0:n-1]
end

"""
    blackman_window(n::Int) -> Vector{Float64}

Create Blackman window function for FFT.

**Formula:**
    w[i] = 0.42 - 0.5×cos(2π×i/(n-1)) + 0.08×cos(4π×i/(n-1))

Provides better sidelobe suppression than Hann.
"""
function blackman_window(n::Int)
    a0 = 0.42
    a1 = 0.5
    a2 = 0.08
    return [a0 - a1*cos(2π*i/(n-1)) + a2*cos(4π*i/(n-1)) for i in 0:n-1]
end

"""
    compute_hhg_spectrum(dipole, dt, ω₀; window=:hann) -> (harmonic_orders, power_spectrum, frequencies)

Compute HHG power spectrum from dipole time series using FFT.

**Algorithm (following Fortran FFT subroutines lines 1027-1356):**
1. Apply window function to reduce spectral leakage
2. Compute FFT: D(ω) = FFT[d(t)]
3. Compute power spectrum: S(ω) = |D(ω)|²
4. Convert to harmonic orders: n = ω/ω₀

**Parameters:**
- `dipole`: Vector of complex dipole moments d(t)
- `dt`: Time step (a.u.)
- `ω₀`: Fundamental laser frequency (a.u.)
- `window`: Window function (:none, :hann, :blackman)

**Returns:**
- `harmonic_orders`: Harmonic order array n = ω/ω₀
- `power_spectrum`: Power spectral density S(ω) = |D(ω)|²
- `frequencies`: Angular frequency array ω (a.u.)

**Example:**
```julia
# Compute spectrum
h_orders, power, freqs = compute_hhg_spectrum(obs.dipole_moment, dt, ω₀, window=:hann)

# Plot HHG spectrum
using Plots
plot(h_orders[1:100], power[1:100], yscale=:log10,
     xlabel="Harmonic Order", ylabel="Power (arb. units)")
```

**References:**
- Fortran: four1 (line 1097), realft (line 1156), dftcor (line 1271)
"""
function compute_hhg_spectrum(dipole::Vector{ComplexF64},
                              dt::Float64,
                              ω₀::Float64;
                              window::Symbol=:hann)
    n_points = length(dipole)

    @info "Computing HHG spectrum" n_points dt ω₀ window

    # Apply window function
    if window == :hann
        w = hann_window(n_points)
        @info "Applied Hann window"
    elseif window == :blackman
        w = blackman_window(n_points)
        @info "Applied Blackman window"
    elseif window == :none
        w = ones(n_points)
        @info "No window applied"
    else
        @warn "Unknown window type: $window, using Hann"
        w = hann_window(n_points)
    end

    # Apply window to dipole
    dipole_windowed = dipole .* w

    # FFT to frequency domain (using FFTW)
    spectrum_complex = fft(dipole_windowed)

    # Power spectrum: S(ω) = |D(ω)|²
    power_spectrum = abs2.(spectrum_complex)

    # Frequency axis: ω = 2π × k / (n × dt)
    frequencies = fftfreq(n_points, 1.0/dt) .* 2π  # Convert to angular frequency

    # Harmonic order axis: n = ω/ω₀
    harmonic_orders = frequencies ./ ω₀

    @info "Spectrum computed" max_harmonic=maximum(abs.(harmonic_orders)) peak_power=maximum(power_spectrum)

    return (harmonic_orders, power_spectrum, frequencies)
end

"""
    save_dipole_to_file(dipole, times, filename)

Save dipole time series to text file (matching Fortran format).

**Format:** time  real(dipole)  imag(dipole)

**References:**
- Fortran output: write(148,*) t/T_cycle, REAL(Dipole_1S), IMAG(Dipole_1S)
- Files: dipolea.txt, dipolel.txt, dipole_1s.txt
"""
function save_dipole_to_file(dipole::Vector{ComplexF64},
                            times::Vector{Float64},
                            filename::String;
                            time_unit::Float64=1.0)
    open(filename, "w") do io
        for (t, d) in zip(times, dipole)
            @printf(io, "%.8e  %.8e  %.8e\n", t/time_unit, real(d), imag(d))
        end
    end
    @info "Dipole saved to file" filename n_points=length(dipole)
end

"""
    save_hhg_spectrum_to_file(harmonic_orders, power_spectrum, filename)

Save HHG spectrum to text file.

**Format:** harmonic_order  power_spectral_density
"""
function save_hhg_spectrum_to_file(harmonic_orders::Vector{Float64},
                                   power_spectrum::Vector{Float64},
                                   filename::String)
    open(filename, "w") do io
        # Header
        println(io, "# Harmonic_Order  Power_Spectral_Density")

        for (h, p) in zip(harmonic_orders, power_spectrum)
            @printf(io, "%.6f  %.12e\n", h, p)
        end
    end
    @info "HHG spectrum saved" filename
end

end  # module HHG
