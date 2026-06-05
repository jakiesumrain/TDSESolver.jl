"""
Unit tests for ResultsIO module.

Tests HDF5 save/load functionality, metadata handling, and query operations.
"""

using Test

# Include the module
include("../../src/io/ResultsIO.jl")
using .ResultsIO

@testset "ResultsIO Tests" begin

    @testset "CalculationResults Construction" begin
        # Minimal results object
        results = CalculationResults(
            nothing, nothing, nothing, nothing,  # Momentum distributions
            [0.1, 0.2, 0.3], [0.1], [0.1], [0.1],  # Grids
            nothing, nothing, nothing, nothing,  # HHG data
            ones(100), Float64[], nothing,  # Validation
            123.45, "2025-01-15T10:30:00"  # Metadata
        )

        @test results.computation_time == 123.45
        @test results.timestamp == "2025-01-15T10:30:00"
        @test length(results.norm_history) == 100
    end

    @testset "Save and Load Round-Trip (Minimal)" begin
        # Create minimal results
        p_grid = collect(0.0:0.1:2.0)
        norm_history = ones(50) .* (1.0 .+ randn(50) * 1e-5)

        results_original = CalculationResults(
            nothing, nothing, nothing, nothing,
            p_grid, Float64[], Float64[], Float64[],
            nothing, nothing, nothing, nothing,
            norm_history, Float64[], nothing,
            10.5, "2025-01-15T12:00:00"
        )

        # Save to temporary file
        temp_file = tempname() * ".h5"
        try
            save_results(results_original, temp_file)

            # Load back
            results_loaded = load_result(temp_file)

            # Verify fields match
            @test results_loaded.p_grid ≈ results_original.p_grid
            @test results_loaded.norm_history ≈ results_original.norm_history
            @test results_loaded.computation_time ≈ results_original.computation_time
            @test results_loaded.timestamp == results_original.timestamp

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Save and Load with Momentum Distributions" begin
        # Create results with momentum distributions
        p_grid = collect(0.0:0.1:2.0)
        px_grid = collect(-1.0:0.1:1.0)
        py_grid = collect(-1.0:0.1:1.0)
        pz_grid = collect(-0.5:0.1:0.5)

        P_radial = exp.(-(p_grid .^ 2))
        P_2D = exp.(-(px_grid .^ 2 .+ py_grid' .^ 2))

        results_original = CalculationResults(
            P_radial, P_2D, nothing, nothing,
            p_grid, px_grid, py_grid, pz_grid,
            nothing, nothing, nothing, nothing,
            ones(100), Float64[], nothing,
            25.3, "2025-01-15T14:30:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results_original, temp_file)
            results_loaded = load_result(temp_file)

            @test results_loaded.momentum_distribution_radial ≈ results_original.momentum_distribution_radial
            @test results_loaded.momentum_distribution_2D_xy ≈ results_original.momentum_distribution_2D_xy
            @test results_loaded.px_grid ≈ results_original.px_grid
            @test results_loaded.py_grid ≈ results_original.py_grid

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Save and Load with HHG Spectrum" begin
        # Create results with HHG data
        photon_energy = collect(0.0:0.1:100.0)
        hhg_spectrum = exp.(-(photon_energy .- 50.0) .^ 2 / 100.0)
        harmonic_orders = [1, 3, 5, 7, 9, 11, 13, 15]
        cutoff = 75.5

        results_original = CalculationResults(
            nothing, nothing, nothing, nothing,
            Float64[], Float64[], Float64[], Float64[],
            hhg_spectrum, photon_energy, harmonic_orders, cutoff,
            ones(200), Float64[], nothing,
            50.2, "2025-01-15T16:00:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results_original, temp_file)
            results_loaded = load_result(temp_file)

            @test results_loaded.hhg_spectrum ≈ results_original.hhg_spectrum
            @test results_loaded.photon_energy_grid ≈ results_original.photon_energy_grid
            @test results_loaded.harmonic_orders == results_original.harmonic_orders
            @test results_loaded.cutoff_position ≈ results_original.cutoff_position

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Save with Metadata" begin
        results = CalculationResults(
            nothing, nothing, nothing, nothing,
            Float64[], Float64[], Float64[], Float64[],
            nothing, nothing, nothing, nothing,
            ones(100), Float64[], 0.998,
            15.7, "2025-01-15T18:00:00"
        )

        metadata = Dict{String, Any}(
            "wavelength_nm" => 800.0,
            "intensity_W_cm2" => 5e14,
            "atom_type" => "helium",
            "pulse_duration_fs" => 10.0
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results, temp_file, metadata=metadata)

            # Verify metadata saved (HDF5 read)
            using HDF5
            h5open(temp_file, "r") do file
                meta_group = file["metadata"]
                @test read(attributes(meta_group), "wavelength_nm") ≈ 800.0
                @test read(attributes(meta_group), "intensity_W_cm2") ≈ 5e14
                @test read(attributes(meta_group), "atom_type") == "helium"
            end

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Query Specific Datasets" begin
        # Create results with multiple datasets
        p_grid = collect(0.0:0.1:1.0)
        P_radial = exp.(-(p_grid .^ 2))
        norm_history = ones(50)

        results = CalculationResults(
            P_radial, nothing, nothing, nothing,
            p_grid, Float64[], Float64[], Float64[],
            nothing, nothing, nothing, nothing,
            norm_history, Float64[], nothing,
            10.0, "2025-01-15T20:00:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results, temp_file)

            # Query individual datasets
            P_queried = query_result(temp_file, :P_radial)
            @test P_queried ≈ P_radial

            p_queried = query_result(temp_file, :p_grid)
            @test p_queried ≈ p_grid

            norm_queried = query_result(temp_file, :norm_history)
            @test norm_queried ≈ norm_history

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Query Nonexistent Dataset" begin
        results = CalculationResults(
            nothing, nothing, nothing, nothing,
            Float64[], Float64[], Float64[], Float64[],
            nothing, nothing, nothing, nothing,
            ones(10), Float64[], nothing,
            5.0, "2025-01-15T21:00:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results, temp_file)

            # Should error when querying dataset that wasn't saved
            @test_throws Exception query_result(temp_file, :P_px_py)

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Query Invalid Dataset Name" begin
        results = CalculationResults(
            nothing, nothing, nothing, nothing,
            Float64[], Float64[], Float64[], Float64[],
            nothing, nothing, nothing, nothing,
            ones(10), Float64[], nothing,
            5.0, "2025-01-15T22:00:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results, temp_file)

            # Should error for unknown dataset symbol
            @test_throws Exception query_result(temp_file, :invalid_dataset)

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Energy Conservation Data" begin
        energy_history = -0.5 .+ randn(100) * 1e-7

        results = CalculationResults(
            nothing, nothing, nothing, nothing,
            Float64[], Float64[], Float64[], Float64[],
            nothing, nothing, nothing, nothing,
            ones(100), energy_history, nothing,
            12.0, "2025-01-15T23:00:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results, temp_file)
            results_loaded = load_result(temp_file)

            @test results_loaded.energy_conservation ≈ results.energy_conservation

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Fortran Correlation Data" begin
        results = CalculationResults(
            nothing, nothing, nothing, nothing,
            Float64[], Float64[], Float64[], Float64[],
            nothing, nothing, nothing, nothing,
            ones(100), Float64[], 0.9995,
            10.0, "2025-01-16T00:00:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results, temp_file)
            results_loaded = load_result(temp_file)

            @test !isnothing(results_loaded.correlation_with_fortran)
            @test results_loaded.correlation_with_fortran ≈ 0.9995

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "File Not Found" begin
        @test_throws Exception load_result("nonexistent_file.h5")
        @test_throws Exception query_result("nonexistent_file.h5", :P_radial)
    end

    @testset "3D Momentum Distribution" begin
        # Create 3D momentum distribution
        px_grid = collect(-1.0:0.2:1.0)
        py_grid = collect(-1.0:0.2:1.0)
        pz_grid = collect(-0.5:0.2:0.5)

        P_3D = zeros(length(px_grid), length(py_grid), length(pz_grid))
        for i in 1:length(px_grid)
            for j in 1:length(py_grid)
                for k in 1:length(pz_grid)
                    p_squared = px_grid[i]^2 + py_grid[j]^2 + pz_grid[k]^2
                    P_3D[i, j, k] = exp(-p_squared)
                end
            end
        end

        results = CalculationResults(
            nothing, nothing, nothing, P_3D,
            Float64[], px_grid, py_grid, pz_grid,
            nothing, nothing, nothing, nothing,
            ones(100), Float64[], nothing,
            30.0, "2025-01-16T01:00:00"
        )

        temp_file = tempname() * ".h5"
        try
            save_results(results, temp_file)
            results_loaded = load_result(temp_file)

            @test results_loaded.momentum_distribution_3D ≈ results.momentum_distribution_3D
            @test size(results_loaded.momentum_distribution_3D) == size(P_3D)

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Complete Results Object" begin
        # Create fully populated results
        p_grid = collect(0.0:0.1:2.0)
        px_grid = collect(-1.0:0.1:1.0)
        py_grid = collect(-1.0:0.1:1.0)
        pz_grid = collect(-0.5:0.1:0.5)

        P_radial = exp.(-(p_grid .^ 2))
        P_2D_xy = exp.(-(px_grid .^ 2 .+ py_grid' .^ 2))
        P_2D_xz = exp.(-(px_grid .^ 2 .+ pz_grid' .^ 2))

        photon_energy = collect(0.0:0.5:100.0)
        hhg_spectrum = exp.(-(photon_energy .- 50.0) .^ 2 / 100.0)
        harmonic_orders = [1, 3, 5, 7, 9]
        cutoff = 80.0

        norm_history = ones(500) .* (1.0 .+ randn(500) * 1e-5)
        energy_history = -0.9 .+ randn(500) * 1e-7

        results = CalculationResults(
            P_radial, P_2D_xy, P_2D_xz, nothing,
            p_grid, px_grid, py_grid, pz_grid,
            hhg_spectrum, photon_energy, harmonic_orders, cutoff,
            norm_history, energy_history, 0.9998,
            125.7, "2025-01-16T02:30:00"
        )

        temp_file = tempname() * ".h5"
        try
            metadata = Dict{String, Any}(
                "wavelength_nm" => 800.0,
                "intensity_W_cm2" => 5e14,
                "atom_type" => "helium"
            )

            save_results(results, temp_file, metadata=metadata)
            results_loaded = load_result(temp_file)

            # Verify all fields
            @test results_loaded.momentum_distribution_radial ≈ results.momentum_distribution_radial
            @test results_loaded.momentum_distribution_2D_xy ≈ results.momentum_distribution_2D_xy
            @test results_loaded.momentum_distribution_2D_xz ≈ results.momentum_distribution_2D_xz
            @test results_loaded.hhg_spectrum ≈ results.hhg_spectrum
            @test results_loaded.harmonic_orders == results.harmonic_orders
            @test results_loaded.cutoff_position ≈ results.cutoff_position
            @test results_loaded.norm_history ≈ results.norm_history
            @test results_loaded.energy_conservation ≈ results.energy_conservation
            @test results_loaded.correlation_with_fortran ≈ results.correlation_with_fortran
            @test results_loaded.computation_time ≈ results.computation_time
            @test results_loaded.timestamp == results.timestamp

        finally
            rm(temp_file, force=true)
        end
    end

end

println("\n✓ ResultsIO tests passed")
