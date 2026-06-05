"""
Unit tests for Validation module.

Tests parameter validation and result diagnostics.
"""

using Test

# Include the module
include("../../src/utils/Validation.jl")
using .Validation

@testset "Validation Tests" begin

    @testset "Physical Parameter Validation" begin
        # Valid parameters
        valid, warnings = validate_physical_parameters(800.0, 5e14, 10.0, [1.0, 0.0, 0.0], 24.6)
        @test valid == true
        @test length(warnings) >= 0  # May have info messages

        # Negative wavelength
        valid, warnings = validate_physical_parameters(-800.0, 5e14, 10.0, [1.0, 0.0, 0.0])
        @test valid == false
        @test any(occursin("ERROR", w) for w in warnings)

        # Negative intensity
        valid, warnings = validate_physical_parameters(800.0, -1e14, 10.0, [1.0, 0.0, 0.0])
        @test valid == false

        # Zero pulse duration
        valid, warnings = validate_physical_parameters(800.0, 5e14, 0.0, [1.0, 0.0, 0.0])
        @test valid == false

        # Unnormalized polarization (should warn but still valid)
        valid, warnings = validate_physical_parameters(800.0, 5e14, 10.0, [2.0, 0.0, 0.0])
        @test any(occursin("WARNING", w) for w in warnings)
    end

    @testset "Keldysh Parameter Warnings" begin
        # Multiphoton regime (γ >> 1) - should warn about low ionization
        # Use low intensity for high Keldysh
        valid, warnings = validate_physical_parameters(800.0, 1e12, 10.0, [1.0, 0.0, 0.0], 24.6)
        @test any(occursin("Keldysh", w) && occursin("multiphoton", w) for w in warnings)

        # Tunneling regime (γ < 1) - should give info
        valid, warnings = validate_physical_parameters(800.0, 1e15, 10.0, [1.0, 0.0, 0.0], 24.6)
        # May have tunneling regime info
    end

    @testset "Numerical Stability Validation" begin
        # Stable time step
        Δt = 0.05  # a.u.
        E_max = 10.0  # a.u.
        valid, warnings = validate_numerical_stability(Δt, E_max)
        @test valid == true

        # Unstable time step (too large)
        Δt_large = 1.0  # Much larger than π/E_max ≈ 0.314
        valid, warnings = validate_numerical_stability(Δt_large, E_max)
        @test valid == false
        @test any(occursin("ERROR", w) && occursin("stability", w) for w in warnings)

        # Just at the limit
        Δt_limit = π / E_max
        valid, warnings = validate_numerical_stability(Δt_limit * 0.99, E_max)
        @test valid == true
    end

    @testset "Calculation Validation Report" begin
        # Good norm conservation
        norm_history = ones(1000) .* (1.0 .+ randn(1000) * 1e-5)  # Small fluctuations
        report = validate_calculation(norm_history)

        @test report.norm_conservation_min > 0.9999
        @test report.norm_conservation_max < 1.0001
        @test report.norm_violation == false

        # Poor norm conservation
        norm_history_bad = ones(1000) .* 0.95  # Significant loss
        report_bad = validate_calculation(norm_history_bad)
        @test report_bad.norm_violation == true
        @test length(report_bad.warnings) > 0
    end

    @testset "Energy Conservation Check" begin
        # Constant energy (good conservation)
        energy_history = ones(100) * (-0.5)  # Hydrogen ground state
        report = validate_calculation(ones(100),
                                     energy_history=energy_history)

        @test report.energy_conservation_error < 1e-10
        @test report.energy_violation == false

        # Varying energy (poor conservation)
        energy_history_bad = -0.5 .+ 0.1 * sin.(range(0, 2π, length=100))
        report_bad = validate_calculation(ones(100),
                                         energy_history=energy_history_bad)
        @test report_bad.energy_conservation_error > 1e-6
        @test report_bad.energy_violation == true
    end

    @testset "Fortran Correlation Check" begin
        # Perfect correlation
        fortran_data = rand(100, 100)
        julia_data = copy(fortran_data)

        report = validate_calculation(ones(10),
                                     fortran_data=fortran_data,
                                     julia_data=julia_data)

        @test !isnothing(report.fortran_correlation)
        @test report.fortran_correlation ≈ 1.0 atol=1e-10
        @test report.benchmark_passed == true

        # Poor correlation
        julia_data_bad = rand(100, 100)
        report_bad = validate_calculation(ones(10),
                                         fortran_data=fortran_data,
                                         julia_data=julia_data_bad)

        @test report_bad.fortran_correlation < 0.999
        @test report_bad.benchmark_passed == false
    end

    @testset "HHG Cutoff Validation" begin
        # Good cutoff agreement
        theory_cutoff = 50.0  # eV
        actual_cutoff = 49.5  # eV (within 1%)

        report = validate_calculation(ones(10),
                                     hhg_cutoff_theory_eV=theory_cutoff,
                                     hhg_cutoff_actual_eV=actual_cutoff)

        @test !isnothing(report.hhg_cutoff_error_percent)
        @test report.hhg_cutoff_error_percent < 5.0

        # Poor cutoff agreement
        actual_cutoff_bad = 55.0  # eV (>10% error)
        report_bad = validate_calculation(ones(10),
                                         hhg_cutoff_theory_eV=theory_cutoff,
                                         hhg_cutoff_actual_eV=actual_cutoff_bad)

        @test report_bad.hhg_cutoff_error_percent > 5.0
        @test any(occursin("cutoff", w) for w in report_bad.warnings)
    end

end

println("\n✓ Validation tests passed")
