"""
    ConfigParser

Module for parsing and validating TOML configuration files for TDSE calculations.

Implements fail-fast validation at load time per FR-006, converting physical parameters
from SI units to atomic units for internal calculations.
"""
module ConfigParser

using TOML
using LinearAlgebra
include("../utils/PhysicalUnits.jl")
using .PhysicalUnits

export SimulationConfiguration, ParameterScan, load_config

"""
    ParameterScan

Specifies parameter sweep for batch calculations.

# Fields
- `parameter::Union{Symbol, Nothing}`: Single parameter to scan (e.g., :intensity)
- `range::Union{StepRangeLen, Nothing}`: Range for single parameter
- `parameters::Union{Vector{Symbol}, Nothing}`: Multiple parameters (2D scan)
- `ranges::Union{Vector{StepRangeLen}, Nothing}`: Ranges for multiple parameters
"""
struct ParameterScan
    parameter::Union{Symbol, Nothing}
    range::Union{StepRangeLen, Nothing}
    parameters::Union{Vector{Symbol}, Nothing}
    ranges::Union{Vector{StepRangeLen}, Nothing}

    function ParameterScan(param, rng, params, rngs)
        # Validate: either single or multi-parameter, not both
        if !isnothing(param) && !isnothing(params)
            error("Cannot specify both 'parameter' and 'parameters' in scan configuration")
        end
        if isnothing(param) && isnothing(params)
            error("Must specify either 'parameter' or 'parameters' in scan configuration")
        end
        new(param, rng, params, rngs)
    end
end

"""
    SimulationConfiguration

Complete specification of a TDSE calculation loaded from TOML file.

# Fields
- `atom_type::Symbol`: Atom type (:hydrogen, :helium, :argon, :neon, :xenon, :custom)
- `custom_potential::Union{String, Function, Nothing}`: Custom potential specification
- `charge::Float64`: Effective nuclear charge (for custom potentials)
- `wavelength_nm::Float64`: Laser wavelength (nm)
- `wavelength_au::Float64`: Laser wavelength (atomic units)
- `intensity_W_cm2::Float64`: Laser intensity (W/cm²)
- `intensity_au::Float64`: Laser intensity (atomic units)
- `pulse_duration_fs::Float64`: Pulse duration FWHM (fs)
- `pulse_duration_au::Float64`: Pulse duration (atomic units)
- `polarization::Vector{Float64}`: Normalized polarization vector [ex, ey, ez]
- `calculation_type::Symbol`: :ionization or :hhg
- `output_file::String`: HDF5 output file path
- `log_file::Union{String, Nothing}`: Optional log file path
- `verbosity::Symbol`: :minimal, :normal, or :verbose
- `numerical_overrides::Dict{Symbol, Any}`: Optional numerical parameter overrides
- `scan::Union{ParameterScan, Nothing}`: Optional parameter scan specification
"""
struct SimulationConfiguration
    # Atomic system
    atom_type::Symbol
    custom_potential::Union{String, Function, Nothing}
    charge::Float64

    # Laser parameters (SI units for user reference)
    wavelength_nm::Float64
    intensity_W_cm2::Float64
    pulse_duration_fs::Float64
    polarization::Vector{Float64}

    # Laser parameters (atomic units for calculation)
    wavelength_au::Float64
    frequency_au::Float64
    intensity_au::Float64
    pulse_duration_au::Float64

    # Calculation settings
    calculation_type::Symbol
    output_file::String
    log_file::Union{String, Nothing}
    verbosity::Symbol

    # Optional overrides
    numerical_overrides::Dict{Symbol, Any}
    scan::Union{ParameterScan, Nothing}
end

"""
    load_config(filepath::String) -> SimulationConfiguration

Load and validate TOML configuration file.

Implements fail-fast validation per FR-006: all parameters are validated at load time,
failing immediately with descriptive error messages if any constraint is violated.

# Arguments
- `filepath::String`: Path to TOML configuration file

# Returns
- `SimulationConfiguration`: Validated configuration object

# Throws
- `ArgumentError`: If any parameter is invalid or missing

# Example
```julia
config = load_config("helium_800nm.toml")
```
"""
function load_config(filepath::String)
    # Parse TOML file
    if !isfile(filepath)
        throw(ArgumentError("Configuration file not found: $filepath"))
    end

    toml_data = TOML.parsefile(filepath)

    # Extract atom configuration
    if !haskey(toml_data, "atom")
        throw(ArgumentError("Missing [atom] section in configuration"))
    end
    atom_section = toml_data["atom"]

    atom_type_str = get(atom_section, "type", nothing)
    if isnothing(atom_type_str)
        throw(ArgumentError("Missing 'type' in [atom] section"))
    end
    atom_type = Symbol(lowercase(atom_type_str))

    # Validate atom type
    valid_atoms = [:hydrogen, :helium, :argon, :neon, :xenon, :custom]
    if !(atom_type in valid_atoms)
        throw(ArgumentError("Invalid atom type '$atom_type'. Must be one of: $(join(valid_atoms, ", "))"))
    end

    # Handle custom potential
    custom_potential = get(atom_section, "potential", nothing)
    charge = get(atom_section, "charge", 1.0)

    if atom_type == :custom && isnothing(custom_potential)
        throw(ArgumentError("Custom atom type requires 'potential' specification"))
    end

    # Extract laser parameters
    if !haskey(toml_data, "laser")
        throw(ArgumentError("Missing [laser] section in configuration"))
    end
    laser_section = toml_data["laser"]

    wavelength_nm = get(laser_section, "wavelength", nothing)
    intensity_W_cm2 = get(laser_section, "intensity", nothing)
    pulse_duration_fs = get(laser_section, "pulse_duration", nothing)
    polarization = get(laser_section, "polarization", [1.0, 0.0, 0.0])

    # Validate laser parameters
    if isnothing(wavelength_nm)
        throw(ArgumentError("Missing 'wavelength' in [laser] section"))
    end
    if wavelength_nm <= 0
        throw(ArgumentError("Wavelength must be positive, got $wavelength_nm nm"))
    end

    if isnothing(intensity_W_cm2)
        throw(ArgumentError("Missing 'intensity' in [laser] section"))
    end
    if intensity_W_cm2 <= 0
        throw(ArgumentError("Intensity must be positive, got $intensity_W_cm2 W/cm²"))
    end

    if isnothing(pulse_duration_fs)
        throw(ArgumentError("Missing 'pulse_duration' in [laser] section"))
    end
    if pulse_duration_fs <= 0
        throw(ArgumentError("Pulse duration must be positive, got $pulse_duration_fs fs"))
    end

    # Validate and normalize polarization
    if length(polarization) != 3
        throw(ArgumentError("Polarization must be 3-element vector [ex, ey, ez], got $(length(polarization)) elements"))
    end
    pol_norm = norm(polarization)
    if pol_norm < 1e-10
        throw(ArgumentError("Polarization vector cannot be zero"))
    end
    polarization = polarization ./ pol_norm  # Normalize

    # Convert to atomic units
    frequency_au = wavelength_nm_to_frequency_au(wavelength_nm)
    intensity_au = intensity_SI_to_au(intensity_W_cm2)
    pulse_duration_au = time_fs_to_au(pulse_duration_fs)
    wavelength_au = wavelength_nm / 0.529177 * 10.0  # nm → a₀

    # Extract calculation settings
    if !haskey(toml_data, "calculation")
        throw(ArgumentError("Missing [calculation] section in configuration"))
    end
    calc_section = toml_data["calculation"]

    calc_type_str = get(calc_section, "type", "ionization")
    calc_type = Symbol(lowercase(calc_type_str))
    if !(calc_type in [:ionization, :hhg])
        throw(ArgumentError("Invalid calculation type '$calc_type'. Must be 'ionization' or 'hhg'"))
    end

    output_file = get(calc_section, "output_file", "results.h5")
    log_file = get(calc_section, "log_file", nothing)

    verbosity_str = get(calc_section, "verbosity", "normal")
    verbosity = Symbol(lowercase(verbosity_str))
    if !(verbosity in [:minimal, :normal, :verbose])
        throw(ArgumentError("Invalid verbosity '$verbosity'. Must be 'minimal', 'normal', or 'verbose'"))
    end

    # Extract numerical overrides
    numerical_overrides = Dict{Symbol, Any}()
    if haskey(toml_data, "numerical")
        for (key, val) in toml_data["numerical"]
            numerical_overrides[Symbol(key)] = val
        end
    end

    # Extract scan configuration
    scan = nothing
    if haskey(toml_data, "scan")
        scan_section = toml_data["scan"]

        if haskey(scan_section, "parameter")
            # Single parameter scan
            param_str = scan_section["parameter"]
            range_str = scan_section["range"]

            # Parse range string "start:step:stop"
            range_parts = split(range_str, ":")
            if length(range_parts) != 3
                throw(ArgumentError("Range must be 'start:step:stop' format, got '$range_str'"))
            end
            start = parse(Float64, range_parts[1])
            step = parse(Float64, range_parts[2])
            stop = parse(Float64, range_parts[3])

            param_symbol = Symbol(replace(param_str, "." => "_"))
            scan = ParameterScan(param_symbol, range(start, step=step, stop=stop), nothing, nothing)

        elseif haskey(scan_section, "parameters")
            # Multi-parameter scan
            param_strs = scan_section["parameters"]
            range_strs = scan_section["ranges"]

            if length(param_strs) != length(range_strs)
                throw(ArgumentError("Number of parameters must match number of ranges"))
            end

            params = [Symbol(replace(p, "." => "_")) for p in param_strs]
            ranges = []
            for range_str in range_strs
                range_parts = split(range_str, ":")
                if length(range_parts) != 3
                    throw(ArgumentError("Range must be 'start:step:stop' format"))
                end
                start = parse(Float64, range_parts[1])
                step = parse(Float64, range_parts[2])
                stop = parse(Float64, range_parts[3])
                push!(ranges, range(start, step=step, stop=stop))
            end

            scan = ParameterScan(nothing, nothing, params, ranges)
        end
    end

    # Construct configuration object
    return SimulationConfiguration(
        atom_type,
        custom_potential,
        charge,
        wavelength_nm,
        intensity_W_cm2,
        pulse_duration_fs,
        polarization,
        wavelength_au,
        frequency_au,
        intensity_au,
        pulse_duration_au,
        calc_type,
        output_file,
        log_file,
        verbosity,
        numerical_overrides,
        scan
    )
end

end  # module ConfigParser
