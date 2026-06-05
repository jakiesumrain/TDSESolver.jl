"""
    PhysicalUnits

Module for converting between SI units and atomic units (a.u.) used in TDSE calculations.

Atomic units conventions:
- Energy: Hartree (Ha) = 27.2114 eV
- Length: Bohr radius a₀ = 0.529177 Å
- Time: ℏ/Ha = 24.1888 as
- Electric field: Ha/(e·a₀)
- Intensity: (1/2) * ε₀ * c * E₀² in W/cm²
"""
module PhysicalUnits

export au_to_eV, eV_to_au, intensity_SI_to_au, intensity_au_to_SI,
       wavelength_nm_to_frequency_au, frequency_au_to_wavelength_nm,
       time_au_to_fs, time_fs_to_au

# Physical constants in atomic units
const HARTREE_TO_EV = 27.2114  # 1 Ha = 27.2114 eV
const BOHR_TO_ANGSTROM = 0.529177  # 1 a₀ = 0.529177 Å
const AU_TIME_TO_FS = 0.024188843  # 1 a.u. time = 24.1888 as = 0.024188843 fs
const SPEED_OF_LIGHT_AU = 137.036  # c in atomic units (fine structure constant⁻¹)
const INTENSITY_AU_TO_SI = 3.5094452e16  # 1 a.u. intensity = 3.5095×10¹⁶ W/cm²

"""
    au_to_eV(energy_au::Float64) -> Float64

Convert energy from atomic units (Hartree) to electron-volts.

# Arguments
- `energy_au::Float64`: Energy in atomic units (Ha)

# Returns
- Energy in electron-volts (eV)

# Example
```julia
julia> au_to_eV(0.5)  # Convert 0.5 Ha to eV
13.6057
```
"""
au_to_eV(energy_au::Float64) = energy_au * HARTREE_TO_EV

"""
    eV_to_au(energy_eV::Float64) -> Float64

Convert energy from electron-volts to atomic units (Hartree).

# Arguments
- `energy_eV::Float64`: Energy in electron-volts (eV)

# Returns
- Energy in atomic units (Ha)
"""
eV_to_au(energy_eV::Float64) = energy_eV / HARTREE_TO_EV

"""
    intensity_SI_to_au(I_W_cm2::Float64) -> Float64

Convert laser intensity from SI units (W/cm²) to atomic units.

# Arguments
- `I_W_cm2::Float64`: Intensity in W/cm²

# Returns
- Intensity in atomic units

# Example
```julia
julia> intensity_SI_to_au(5e14)  # 5×10¹⁴ W/cm² to a.u.
0.01425
```
"""
intensity_SI_to_au(I_W_cm2::Float64) = I_W_cm2 / INTENSITY_AU_TO_SI

"""
    intensity_au_to_SI(I_au::Float64) -> Float64

Convert laser intensity from atomic units to SI units (W/cm²).

# Arguments
- `I_au::Float64`: Intensity in atomic units

# Returns
- Intensity in W/cm²
"""
intensity_au_to_SI(I_au::Float64) = I_au * INTENSITY_AU_TO_SI

"""
    wavelength_nm_to_frequency_au(λ_nm::Float64) -> Float64

Convert wavelength in nanometers to angular frequency in atomic units.

# Arguments
- `λ_nm::Float64`: Wavelength in nanometers (nm)

# Returns
- Angular frequency ω in atomic units (Ha/ℏ)

# Example
```julia
julia> wavelength_nm_to_frequency_au(800.0)  # 800 nm laser
0.0569903
```

# Formula
ω [a.u.] = (2πc [a.u.]) / (λ [nm] / 0.529177 Å/a₀ × 10⁻¹)
"""
function wavelength_nm_to_frequency_au(λ_nm::Float64)
    # Convert nm to a.u. (Bohr radii)
    λ_au = λ_nm / BOHR_TO_ANGSTROM * 10.0  # nm → Å → a₀
    # ω = 2πc/λ in atomic units
    ω_au = 2π * SPEED_OF_LIGHT_AU / λ_au
    return ω_au
end

"""
    frequency_au_to_wavelength_nm(ω_au::Float64) -> Float64

Convert angular frequency in atomic units to wavelength in nanometers.

# Arguments
- `ω_au::Float64`: Angular frequency in atomic units (Ha/ℏ)

# Returns
- Wavelength in nanometers (nm)
"""
function frequency_au_to_wavelength_nm(ω_au::Float64)
    # λ = 2πc/ω
    λ_au = 2π * SPEED_OF_LIGHT_AU / ω_au
    # Convert a.u. to nm
    λ_nm = λ_au * BOHR_TO_ANGSTROM / 10.0  # a₀ → Å → nm
    return λ_nm
end

"""
    time_au_to_fs(t_au::Float64) -> Float64

Convert time from atomic units to femtoseconds.

# Arguments
- `t_au::Float64`: Time in atomic units (ℏ/Ha)

# Returns
- Time in femtoseconds (fs)
"""
time_au_to_fs(t_au::Float64) = t_au * AU_TIME_TO_FS

"""
    time_fs_to_au(t_fs::Float64) -> Float64

Convert time from femtoseconds to atomic units.

# Arguments
- `t_fs::Float64`: Time in femtoseconds (fs)

# Returns
- Time in atomic units (ℏ/Ha)
"""
time_fs_to_au(t_fs::Float64) = t_fs / AU_TIME_TO_FS

end  # module PhysicalUnits
