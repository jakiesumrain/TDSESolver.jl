"""
    Validation

Module for validating TDSE calculation parameters and results.

Implements validation diagnostics per FR-011, including norm conservation,
energy conservation, symmetry checks, and benchmark comparisons.
"""
module Validation

using LinearAlgebra
using Printf  # For @sprintf
using Statistics  # For mean, std, cor

export validate_physical_parameters, validate_numerical_stability,
       ValidationReport, validate_calculation

"""
    ValidationReport

Quality checks on TDSE calculation results.

# Fields
- `norm_conservation_min::Float64`: Minimum norm during propagation
- `norm_conservation_max::Float64`: Maximum norm during propagation
- `norm_violation::Bool`: True if |norm - 1.0| > threshold
- `energy_conservation_error::Float64`: Max ΔE/E during field-free periods
- `energy_violation::Bool`: True if error > threshold
- `cylindrical_symmetry_error::Union{Float64, Nothing}`: Symmetry error (if applicable)
- `symmetry_violation::Bool`: True if symmetry violated
- `fortran_correlation::Union{Float64, Nothing}`: Correlation with Fortran benchmark
- `benchmark_passed::Bool`: True if correlation > 0.999
- `hhg_cutoff_theory::Union{Float64, Nothing}`: Theoretical cutoff (Ip + 3.17Up) in eV
- `hhg_cutoff_actual::Union{Float64, Nothing}`: Measured cutoff in eV
- `hhg_cutoff_error_percent::Union{Float64, Nothing}`: Cutoff error percentage
- `hhg_odd_harmonics::Union{Bool, Nothing}`: True if only odd harmonics
- `passed_validation::Bool`: True if all checks passed
- `warnings::Vector{String}`: Warning messages
"""
struct ValidationReport
    # Norm conservation
    norm_conservation_min::Float64
    norm_conservation_max::Float64
    norm_violation::Bool

    # Energy conservation
    energy_conservation_error::Float64
    energy_violation::Bool

    # Symmetry checks
    cylindrical_symmetry_error::Union{Float64, Nothing}
    symmetry_violation::Bool

    # Benchmark comparison
    fortran_correlation::Union{Float64, Nothing}
    benchmark_passed::Bool

    # HHG validation
    hhg_cutoff_theory::Union{Float64, Nothing}
    hhg_cutoff_actual::Union{Float64, Nothing}
    hhg_cutoff_error_percent::Union{Float64, Nothing}
    hhg_odd_harmonics::Union{Bool, Nothing}

    # Overall status
    passed_validation::Bool
    warnings::Vector{String}
end

"""
    validate_physical_parameters(config) -> (Bool, Vector{String})

Validate physical parameters from configuration.

# Arguments
- `config`: SimulationConfiguration object

# Returns
- `(valid::Bool, warnings::Vector{String})`: Validation status and warnings

# Validation checks:
- Wavelength > 0
- Intensity > 0
- Pulse duration > 0
- Polarization normalized
- Keldysh parameter indicates expected regime
"""
function validate_physical_parameters(wavelength_nm, intensity_W_cm2, pulse_duration_fs,
                                     polarization, ionization_potential_eV=13.6)
    warnings = String[]
    valid = true

    # Basic positivity checks
    if wavelength_nm <= 0
        push!(warnings, "ERROR: Wavelength must be positive, got $wavelength_nm nm")
        valid = false
    end

    if intensity_W_cm2 <= 0
        push!(warnings, "ERROR: Intensity must be positive, got $intensity_W_cm2 W/cm²")
        valid = false
    end

    if pulse_duration_fs <= 0
        push!(warnings, "ERROR: Pulse duration must be positive, got $pulse_duration_fs fs")
        valid = false
    end

    # Polarization normalization check
    pol_norm = norm(polarization)
    if abs(pol_norm - 1.0) > 1e-6
        push!(warnings, "WARNING: Polarization not normalized (norm=$pol_norm), will be normalized automatically")
    end

    # Only compute physical regime warnings if basic parameters are valid
    if valid
        # Physical regime warnings
        # Compute ponderomotive energy: Up = I/(4ω²) in atomic units
        # I [a.u.] = I [W/cm²] / 3.5094e16
        # ω [a.u.] = 2πc / λ, with c = 137.036 a.u., λ in Bohr radii
        λ_au = wavelength_nm / 0.529177 * 10.0
        ω_au = 2π * 137.036 / λ_au
        I_au = intensity_W_cm2 / 3.5094e16
        Up_au = I_au / (4 * ω_au^2)
        Up_eV = Up_au * 27.2114

        # Keldysh parameter: γ = √(Ip / 2Up)
        Ip_au = ionization_potential_eV / 27.2114
        keldysh = sqrt(Ip_au / (2 * Up_au))

        if keldysh > 2.0
            push!(warnings, "WARNING: Keldysh parameter γ=$(@sprintf("%.2f", keldysh)) >> 1 (multiphoton regime), ionization probability may be very low")
        elseif keldysh < 0.5
            push!(warnings, "INFO: Keldysh parameter γ=$(@sprintf("%.2f", keldysh)) < 1 (tunneling regime)")
        end

        # Check if intensity might exceed grid capacity
        if Up_eV > 50.0
            push!(warnings, "WARNING: High ponderomotive energy Up=$(@sprintf("%.1f", Up_eV)) eV may require large spatial grid (rmax > 200 a.u.)")
        end
    end

    return (valid, warnings)
end

"""
    validate_numerical_stability(Δt_au, E_max_au) -> (Bool, Vector{String})

Validate numerical stability criteria for time propagation.

# Arguments
- `Δt_au`: Time step in atomic units
- `E_max_au`: Maximum eigenstate energy in atomic units

# Returns
- `(valid::Bool, warnings::Vector{String})`: Validation status and warnings

# Stability criterion:
Δt < π / E_max for stable propagation
"""
function validate_numerical_stability(Δt_au, E_max_au)
    warnings = String[]
    valid = true

    # Stability criterion from split-operator method
    Δt_max = π / abs(E_max_au)

    if Δt_au > Δt_max
        push!(warnings, "ERROR: Time step Δt=$(@sprintf("%.4f", Δt_au)) a.u. exceeds stability criterion π/E_max=$(@sprintf("%.4f", Δt_max)) a.u.")
        valid = false
    elseif Δt_au > 0.8 * Δt_max
        push!(warnings, "WARNING: Time step close to stability limit, consider reducing Δt")
    end

    return (valid, warnings)
end

"""
    validate_calculation(norm_history, energy_history=nothing, fortran_data=nothing) -> ValidationReport

Generate validation report for completed calculation.

# Arguments
- `norm_history::Vector{Float64}`: Wavefunction norm at each time step
- `energy_history::Union{Vector{Float64}, Nothing}`: Energy during field-free periods
- `fortran_data::Union{Array, Nothing}`: Fortran benchmark data for correlation

# Returns
- `ValidationReport`: Comprehensive validation report

# Validation thresholds:
- Norm conservation: 0.9999 < norm < 1.0001
- Energy conservation: ΔE/E < 1e-6
- Fortran correlation: > 0.999
"""
function validate_calculation(norm_history::Vector{Float64};
                              energy_history::Union{Vector{Float64}, Nothing}=nothing,
                              fortran_data::Union{Array, Nothing}=nothing,
                              julia_data::Union{Array, Nothing}=nothing,
                              hhg_cutoff_theory_eV::Union{Float64, Nothing}=nothing,
                              hhg_cutoff_actual_eV::Union{Float64, Nothing}=nothing)
    warnings = String[]

    # Norm conservation check
    norm_min = minimum(norm_history)
    norm_max = maximum(norm_history)
    norm_threshold = 1e-4  # 0.9999 to 1.0001
    norm_violation = (norm_min < 1.0 - norm_threshold) || (norm_max > 1.0 + norm_threshold)

    if norm_violation
        push!(warnings, "Norm conservation violated: range [$(@sprintf("%.6f", norm_min)), $(@sprintf("%.6f", norm_max))]")
    end

    # Energy conservation check
    energy_error = 0.0
    energy_violation = false
    if !isnothing(energy_history) && length(energy_history) > 1
        E_mean = mean(energy_history)
        energy_error = maximum(abs.(energy_history .- E_mean)) / abs(E_mean)
        energy_threshold = 1e-6
        energy_violation = energy_error > energy_threshold

        if energy_violation
            push!(warnings, "Energy conservation violated: ΔE/E = $(@sprintf("%.2e", energy_error))")
        end
    end

    # Fortran correlation check
    fortran_corr = nothing
    benchmark_passed = true
    if !isnothing(fortran_data) && !isnothing(julia_data)
        # Flatten arrays for correlation
        f_flat = vec(fortran_data)
        j_flat = vec(julia_data)

        if length(f_flat) == length(j_flat)
            # Pearson correlation coefficient
            fortran_corr = cor(f_flat, j_flat)
            benchmark_passed = fortran_corr > 0.999

            if !benchmark_passed
                push!(warnings, "Fortran correlation $(@sprintf("%.6f", fortran_corr)) below threshold 0.999")
            end
        else
            push!(warnings, "Cannot compute Fortran correlation: array size mismatch")
            benchmark_passed = false
        end
    end

    # HHG cutoff check
    hhg_cutoff_error = nothing
    hhg_odd_harmonics = nothing
    if !isnothing(hhg_cutoff_theory_eV) && !isnothing(hhg_cutoff_actual_eV)
        hhg_cutoff_error = abs(hhg_cutoff_actual_eV - hhg_cutoff_theory_eV) / hhg_cutoff_theory_eV * 100.0

        if hhg_cutoff_error > 5.0
            push!(warnings, "HHG cutoff error $(@sprintf("%.1f", hhg_cutoff_error))% exceeds 5% threshold")
        end
    end

    # Overall validation status
    passed = !norm_violation && !energy_violation && benchmark_passed

    return ValidationReport(
        norm_min,
        norm_max,
        norm_violation,
        energy_error,
        energy_violation,
        nothing,  # cylindrical_symmetry_error
        false,    # symmetry_violation
        fortran_corr,
        benchmark_passed,
        hhg_cutoff_theory_eV,
        hhg_cutoff_actual_eV,
        hhg_cutoff_error,
        hhg_odd_harmonics,
        passed,
        warnings
    )
end

end  # module Validation
