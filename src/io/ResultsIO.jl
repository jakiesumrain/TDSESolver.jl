"""
    ResultsIO

Module for saving and loading TDSE calculation results using HDF5 format.

Implements result persistence per data-model.md CalculationResults entity,
with metadata and structured organization for reproducibility.
"""
module ResultsIO

using HDF5
using Dates

export CalculationResults, save_results, load_result, query_result

"""
    CalculationResults

Complete TDSE calculation results with validation diagnostics and metadata.

# Fields
- `momentum_distribution_radial::Union{Vector{Float64}, Nothing}`: P(p)
- `momentum_distribution_2D_xy::Union{Matrix{Float64}, Nothing}`: P(px, py)
- `momentum_distribution_2D_xz::Union{Matrix{Float64}, Nothing}`: P(px, pz)
- `momentum_distribution_3D::Union{Array{Float64, 3}, Nothing}`: P(px, py, pz)
- `p_grid::Vector{Float64}`: Radial momentum grid
- `px_grid::Vector{Float64}`: Momentum x grid
- `py_grid::Vector{Float64}`: Momentum y grid
- `pz_grid::Vector{Float64}`: Momentum z grid
- `hhg_spectrum::Union{Vector{Float64}, Nothing}`: S(ω)
- `photon_energy_grid::Union{Vector{Float64}, Nothing}`: ω in eV
- `harmonic_orders::Union{Vector{Int}, Nothing}`: Identified harmonic orders
- `cutoff_position::Union{Float64, Nothing}`: Cutoff photon energy (eV)
- `norm_history::Vector{Float64}`: Norm at each time step
- `energy_conservation::Vector{Float64}`: Energy during field-free periods
- `correlation_with_fortran::Union{Float64, Nothing}`: Benchmark correlation
- `computation_time::Float64`: Wall-clock time (seconds)
- `timestamp::String`: ISO 8601 timestamp
"""
struct CalculationResults
    # Momentum distributions (ionization)
    momentum_distribution_radial::Union{Vector{Float64}, Nothing}
    momentum_distribution_2D_xy::Union{Matrix{Float64}, Nothing}
    momentum_distribution_2D_xz::Union{Matrix{Float64}, Nothing}
    momentum_distribution_3D::Union{Array{Float64, 3}, Nothing}

    # Momentum grids
    p_grid::Vector{Float64}
    px_grid::Vector{Float64}
    py_grid::Vector{Float64}
    pz_grid::Vector{Float64}

    # HHG spectrum (HHG calculations)
    hhg_spectrum::Union{Vector{Float64}, Nothing}
    photon_energy_grid::Union{Vector{Float64}, Nothing}
    harmonic_orders::Union{Vector{Int}, Nothing}
    cutoff_position::Union{Float64, Nothing}

    # Validation diagnostics
    norm_history::Vector{Float64}
    energy_conservation::Vector{Float64}
    correlation_with_fortran::Union{Float64, Nothing}

    # Metadata
    computation_time::Float64
    timestamp::String
end

"""
    save_results(results::CalculationResults, filepath::String; metadata::Dict=Dict())

Save calculation results to HDF5 file with metadata.

# Arguments
- `results::CalculationResults`: Results to save
- `filepath::String`: HDF5 output file path
- `metadata::Dict`: Additional metadata (wavelength, intensity, etc.)

# HDF5 Structure:
```
/momentum_distributions/
    P_radial [dataset]
    P_px_py [dataset]
    P_px_pz [dataset]
    P_3D [dataset]
    p_grid [dataset]
    px_grid, py_grid, pz_grid [datasets]
/hhg_spectrum/ (if applicable)
    spectrum [dataset]
    photon_energy [dataset]
    harmonic_orders [dataset]
/validation_diagnostics/
    norm_history [dataset]
    energy_conservation [dataset]
    correlation_with_fortran [attribute]
/metadata/
    wavelength_nm [attribute]
    intensity_W_cm2 [attribute]
    computation_time_s [attribute]
    timestamp [attribute]
```

# Example
```julia
save_results(results, "helium_800nm_results.h5",
            metadata=Dict("wavelength_nm" => 800.0,
                         "intensity_W_cm2" => 5e14))
```
"""
function save_results(results::CalculationResults, filepath::String;
                     metadata::Dict{String, Any}=Dict{String, Any}())
    h5open(filepath, "w") do file
        # Create groups
        mom_group = create_group(file, "momentum_distributions")
        diag_group = create_group(file, "validation_diagnostics")
        meta_group = create_group(file, "metadata")

        # Save momentum distributions
        if !isnothing(results.momentum_distribution_radial)
            mom_group["P_radial"] = results.momentum_distribution_radial
        end

        if !isnothing(results.momentum_distribution_2D_xy)
            mom_group["P_px_py"] = results.momentum_distribution_2D_xy
        end

        if !isnothing(results.momentum_distribution_2D_xz)
            mom_group["P_px_pz"] = results.momentum_distribution_2D_xz
        end

        if !isnothing(results.momentum_distribution_3D)
            mom_group["P_3D"] = results.momentum_distribution_3D
        end

        # Save momentum grids
        mom_group["p_grid"] = results.p_grid
        mom_group["px_grid"] = results.px_grid
        mom_group["py_grid"] = results.py_grid
        mom_group["pz_grid"] = results.pz_grid

        # Save HHG spectrum if present
        if !isnothing(results.hhg_spectrum)
            hhg_group = create_group(file, "hhg_spectrum")
            hhg_group["spectrum"] = results.hhg_spectrum

            if !isnothing(results.photon_energy_grid)
                hhg_group["photon_energy"] = results.photon_energy_grid
            end

            if !isnothing(results.harmonic_orders)
                hhg_group["harmonic_orders"] = results.harmonic_orders
            end

            if !isnothing(results.cutoff_position)
                attributes(hhg_group)["cutoff_eV"] = results.cutoff_position
            end
        end

        # Save validation diagnostics
        diag_group["norm_history"] = results.norm_history

        if !isempty(results.energy_conservation)
            diag_group["energy_conservation"] = results.energy_conservation
        end

        if !isnothing(results.correlation_with_fortran)
            attributes(diag_group)["fortran_correlation"] = results.correlation_with_fortran
        end

        # Save metadata
        attributes(meta_group)["computation_time_s"] = results.computation_time
        attributes(meta_group)["timestamp"] = results.timestamp

        # Save additional metadata from user
        for (key, value) in metadata
            if value isa Number || value isa String
                attributes(meta_group)[key] = value
            elseif value isa Vector
                meta_group[key] = value
            else
                @warn "Skipping metadata key '$key' with unsupported type $(typeof(value))"
            end
        end

        println("Results saved to: $filepath")
    end
end

"""
    load_result(filepath::String) -> CalculationResults

Load calculation results from HDF5 file.

# Arguments
- `filepath::String`: Path to HDF5 results file

# Returns
- `CalculationResults`: Reconstructed results object

# Example
```julia
results = load_result("helium_800nm_results.h5")
P_px_py = results.momentum_distribution_2D_xy
```
"""
function load_result(filepath::String)
    if !isfile(filepath)
        error("Results file not found: $filepath")
    end

    h5open(filepath, "r") do file
        # Load momentum distributions
        mom_group = file["momentum_distributions"]

        P_radial = haskey(mom_group, "P_radial") ? read(mom_group, "P_radial") : nothing
        P_px_py = haskey(mom_group, "P_px_py") ? read(mom_group, "P_px_py") : nothing
        P_px_pz = haskey(mom_group, "P_px_pz") ? read(mom_group, "P_px_pz") : nothing
        P_3D = haskey(mom_group, "P_3D") ? read(mom_group, "P_3D") : nothing

        # Load momentum grids
        p_grid = read(mom_group, "p_grid")
        px_grid = read(mom_group, "px_grid")
        py_grid = read(mom_group, "py_grid")
        pz_grid = read(mom_group, "pz_grid")

        # Load HHG spectrum if present
        hhg_spectrum = nothing
        photon_energy = nothing
        harmonic_orders = nothing
        cutoff = nothing

        if haskey(file, "hhg_spectrum")
            hhg_group = file["hhg_spectrum"]
            hhg_spectrum = read(hhg_group, "spectrum")

            if haskey(hhg_group, "photon_energy")
                photon_energy = read(hhg_group, "photon_energy")
            end

            if haskey(hhg_group, "harmonic_orders")
                harmonic_orders = read(hhg_group, "harmonic_orders")
            end

            if haskey(attributes(hhg_group), "cutoff_eV")
                cutoff = read(attributes(hhg_group), "cutoff_eV")
            end
        end

        # Load validation diagnostics
        diag_group = file["validation_diagnostics"]
        norm_history = read(diag_group, "norm_history")

        energy_cons = haskey(diag_group, "energy_conservation") ?
                     read(diag_group, "energy_conservation") : Float64[]

        fortran_corr = haskey(attributes(diag_group), "fortran_correlation") ?
                      read(attributes(diag_group), "fortran_correlation") : nothing

        # Load metadata
        meta_group = file["metadata"]
        comp_time = read(attributes(meta_group), "computation_time_s")
        timestamp = read(attributes(meta_group), "timestamp")

        return CalculationResults(
            P_radial, P_px_py, P_px_pz, P_3D,
            p_grid, px_grid, py_grid, pz_grid,
            hhg_spectrum, photon_energy, harmonic_orders, cutoff,
            norm_history, energy_cons, fortran_corr,
            comp_time, timestamp
        )
    end
end

"""
    query_result(filepath::String, dataset::Symbol) -> Array

Query specific dataset from HDF5 file without loading entire file.

# Arguments
- `filepath::String`: Path to HDF5 results file
- `dataset::Symbol`: Dataset to query (e.g., :P_px_py, :norm_history)

# Returns
- Dataset array

# Example
```julia
P_xy = query_result("results.h5", :P_px_py)
norm = query_result("results.h5", :norm_history)
```
"""
function query_result(filepath::String, dataset::Symbol)
    dataset_paths = Dict(
        :P_radial => "momentum_distributions/P_radial",
        :P_px_py => "momentum_distributions/P_px_py",
        :P_px_pz => "momentum_distributions/P_px_pz",
        :P_3D => "momentum_distributions/P_3D",
        :p_grid => "momentum_distributions/p_grid",
        :px_grid => "momentum_distributions/px_grid",
        :py_grid => "momentum_distributions/py_grid",
        :pz_grid => "momentum_distributions/pz_grid",
        :hhg_spectrum => "hhg_spectrum/spectrum",
        :photon_energy => "hhg_spectrum/photon_energy",
        :norm_history => "validation_diagnostics/norm_history",
        :energy_conservation => "validation_diagnostics/energy_conservation"
    )

    if !haskey(dataset_paths, dataset)
        error("Unknown dataset: $dataset. Available: $(keys(dataset_paths))")
    end

    path = dataset_paths[dataset]

    h5open(filepath, "r") do file
        if !haskey(file, path)
            error("Dataset '$path' not found in file")
        end
        return read(file, path)
    end
end

end  # module ResultsIO
