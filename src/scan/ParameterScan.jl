"""
    ParameterScan

Module for parallel parameter scanning across multiple TDSE simulations.

Enables systematic exploration of parameter space (intensity, wavelength, etc.)
with automatic result aggregation and progress monitoring.

# Usage

```julia
# Load configuration with [scan] section
config = load_config("intensity_scan.toml")

# Run scan with parallel execution
results = run_parameter_scan(config, verbose=true)

# Access individual results
for (i, result) in enumerate(results.simulation_results)
    println("Result \$i: Ionization = \$(result.observables.ionization_probs[end])")
end
```

# Implementation Notes

Based on Fortran approach where multiple simulations are run with varying parameters,
then results are aggregated for analysis (plotting ionization vs intensity, etc.).
"""
module ParameterScan

using Distributed
using Printf
using ..ConfigParser: SimulationConfiguration, ParameterScan as ScanConfig
using ..Simulation: SimulationParams, SimulationResults, run_simulation, create_default_params

export run_parameter_scan, ScanResults, print_scan_summary

"""
    ScanResults

Container for parameter scan results.

# Fields
- `scan_config::ScanConfig`: Scan configuration
- `parameter_values::Vector{Float64}`: Parameter values used in scan
- `simulation_results::Vector{SimulationResults}`: Results for each parameter value
- `total_time::Float64`: Total wall time for scan (seconds)
"""
struct ScanResults
    scan_config::ScanConfig
    parameter_values::Vector{Float64}
    simulation_results::Vector{SimulationResults}
    total_time::Float64
end

"""
    run_parameter_scan(config::SimulationConfiguration;
                       verbose::Bool=true,
                       parallel::Bool=false) -> ScanResults

Execute parameter scan over specified range.

# Arguments
- `config::SimulationConfiguration`: Configuration with scan section
- `verbose::Bool=true`: Print progress information
- `parallel::Bool=false`: Use parallel execution (requires Distributed.jl setup)

# Returns
- `ScanResults`: Container with all simulation results

# Example
```julia
config = load_config("intensity_scan.toml")
results = run_parameter_scan(config, verbose=true)
```

# Performance
- Serial: Simulations run one at a time
- Parallel: Simulations distributed across available workers (use `addprocs(n)`)

# Supported Parameters
Currently supports scanning over:
- `laser_intensity`: Laser intensity (W/cm²)
- `laser_wavelength`: Laser wavelength (nm)
- `laser_duration`: Pulse duration (fs)
- `laser_cep`: Carrier-envelope phase (radians)
"""
function run_parameter_scan(config::SimulationConfiguration;
                            verbose::Bool=true,
                            parallel::Bool=false)
    if isnothing(config.scan)
        error("Configuration does not contain scan section")
    end

    scan = config.scan
    start_time = time()

    # Determine parameter and range
    if !isnothing(scan.parameter)
        # Single parameter scan
        param_name = scan.parameter
        param_range = scan.range
    else
        error("Multi-parameter scans not yet implemented. Use single parameter scan.")
    end

    # Convert range to vector
    parameter_values = collect(param_range)
    n_sims = length(parameter_values)

    if verbose
        println("="^80)
        println("PARAMETER SCAN")
        println("="^80)
        println("Parameter: $(String(param_name))")
        println("Range: $(parameter_values[1]) to $(parameter_values[end])")
        println("Number of simulations: $n_sims")
        println("Parallel execution: $(parallel ? "enabled" : "disabled")")
        println("="^80)
        println()
    end

    # Create base parameters from configuration
    base_params = config_to_params(config)

    # Run simulations
    simulation_results = if parallel
        # Parallel execution
        if verbose
            println("Running simulations in parallel across $(nworkers()) workers...")
        end
        run_scan_parallel(base_params, param_name, parameter_values, verbose)
    else
        # Serial execution
        run_scan_serial(base_params, param_name, parameter_values, verbose)
    end

    elapsed_time = time() - start_time

    if verbose
        println()
        println("="^80)
        println("SCAN COMPLETE")
        println("="^80)
        println("Total time: $(@sprintf("%.2f", elapsed_time)) seconds")
        println("Average time per simulation: $(@sprintf("%.2f", elapsed_time/n_sims)) seconds")
        println("="^80)
    end

    return ScanResults(scan, parameter_values, simulation_results, elapsed_time)
end

"""
    run_scan_serial(base_params, param_name, values, verbose) -> Vector{SimulationResults}

Run parameter scan serially (one simulation at a time).
"""
function run_scan_serial(base_params::SimulationParams,
                        param_name::Symbol,
                        values::Vector{Float64},
                        verbose::Bool)
    results = Vector{SimulationResults}(undef, length(values))

    for (i, value) in enumerate(values)
        if verbose
            println("Simulation $i/$(length(values)): $(String(param_name)) = $value")
        end

        # Create modified parameters
        params = modify_param(base_params, param_name, value)

        # Run simulation
        result = run_simulation(params, verbose=false)
        results[i] = result

        if verbose
            final_ionization = result.observables.ionization_probs[end]
            println("  → Ionization probability: $(@sprintf("%.6f", final_ionization))")
            println()
        end
    end

    return results
end

"""
    run_scan_parallel(base_params, param_name, values, verbose) -> Vector{SimulationResults}

Run parameter scan in parallel across available workers.

# Requirements
Must call `addprocs(n)` before running to add workers.
"""
function run_scan_parallel(base_params::SimulationParams,
                          param_name::Symbol,
                          values::Vector{Float64},
                          verbose::Bool)
    # Parallel map over parameter values
    results = pmap(values) do value
        params = modify_param(base_params, param_name, value)
        run_simulation(params, verbose=false)
    end

    return results
end

"""
    modify_param(params, name, value) -> SimulationParams

Create new SimulationParams with specified parameter modified.
"""
function modify_param(params::SimulationParams, name::Symbol, value::Float64)
    # Convert param struct to named tuple, modify value, create new struct
    param_dict = Dict{Symbol, Any}()

    for field in fieldnames(SimulationParams)
        if field == name
            param_dict[field] = value
        else
            param_dict[field] = getfield(params, field)
        end
    end

    return SimulationParams(; param_dict...)
end

"""
    config_to_params(config::SimulationConfiguration) -> SimulationParams

Convert configuration to simulation parameters.
"""
function config_to_params(config::SimulationConfiguration)
    # Build keyword arguments dictionary
    kwargs = Dict{Symbol, Any}(
        :atom => config.atom_type,
        :custom_potential => config.custom_potential,
        :laser_wavelength => config.wavelength_nm,
        :laser_intensity => config.intensity_W_cm2,
        :laser_duration => config.pulse_duration_fs,
        :laser_polarization => config.polarization
    )

    # Add numerical overrides if present
    if !isempty(config.numerical_overrides)
        merge!(kwargs, config.numerical_overrides)
    end

    # Create parameters with all keyword arguments
    return create_default_params(; kwargs...)
end

"""
    print_scan_summary(results::ScanResults)

Print summary of scan results showing key observables vs parameter.
"""
function print_scan_summary(results::ScanResults)
    println("="^80)
    println("SCAN RESULTS SUMMARY")
    println("="^80)
    println()

    param_name = !isnothing(results.scan_config.parameter) ?
                 String(results.scan_config.parameter) : "parameters"

    println("Parameter: $param_name")
    println("Number of simulations: $(length(results.parameter_values))")
    println()

    # Print table header
    println(@sprintf("%-20s  %-15s  %-15s", param_name, "P_ionization", "Wall time (s)"))
    println("-"^80)

    # Print each result
    for (i, (value, result)) in enumerate(zip(results.parameter_values, results.simulation_results))
        ionization = result.observables.ionization_probs[end]
        wall_time = result.wall_time

        println(@sprintf("%-20.6e  %-15.6f  %-15.2f", value, ionization, wall_time))
    end

    println("-"^80)
    println("Total scan time: $(@sprintf("%.2f", results.total_time)) seconds")
    println("="^80)
end

end  # module ParameterScan
