"""
Unit tests for PhysicalUnits module.

Tests unit conversion functions and physical constants.
"""

using Test

# Include the module
include("../../src/utils/PhysicalUnits.jl")
using .PhysicalUnits

@testset "PhysicalUnits Tests" begin

    @testset "Energy Conversions" begin
        # Test Ha ↔ eV conversions
        @test au_to_eV(1.0) ≈ 27.2114 atol=1e-4
        @test eV_to_au(27.2114) ≈ 1.0 atol=1e-6

        # Round-trip conversion
        E_au = 0.5
        @test eV_to_au(au_to_eV(E_au)) ≈ E_au atol=1e-10

        # Known values
        @test au_to_eV(0.5) ≈ 13.6057 atol=1e-4  # Hydrogen ionization potential
        @test au_to_eV(0.9034) ≈ 24.5828 atol=1e-3  # Helium ionization potential (24.5874 eV)
    end

    @testset "Intensity Conversions" begin
        # Test W/cm² ↔ a.u. conversions
        I_SI = 5e14  # W/cm²
        I_au = intensity_SI_to_au(I_SI)
        @test I_au > 0
        @test intensity_au_to_SI(I_au) ≈ I_SI rtol=1e-10

        # Known conversion factor
        @test intensity_SI_to_au(3.5094452e16) ≈ 1.0 atol=1e-3

        # Round-trip
        I_test = 1e15
        @test intensity_au_to_SI(intensity_SI_to_au(I_test)) ≈ I_test rtol=1e-10
    end

    @testset "Wavelength/Frequency Conversions" begin
        # 800 nm laser (common Ti:Sapphire wavelength)
        λ_nm = 800.0
        ω_au = wavelength_nm_to_frequency_au(λ_nm)

        # Check positive frequency
        @test ω_au > 0

        # Check order of magnitude (800nm ≈ 1.55 eV ≈ 0.057 Ha)
        ω_eV = au_to_eV(ω_au)
        @test ω_eV ≈ 1.55 atol=0.01

        # Round-trip conversion
        λ_roundtrip = frequency_au_to_wavelength_nm(ω_au)
        @test λ_roundtrip ≈ λ_nm rtol=1e-6

        # Common laser wavelengths
        @test wavelength_nm_to_frequency_au(400.0) ≈ 2 * wavelength_nm_to_frequency_au(800.0) rtol=0.01
    end

    @testset "Time Conversions" begin
        # Test fs ↔ a.u. conversions
        t_fs = 10.0  # femtoseconds
        t_au = time_fs_to_au(t_fs)

        @test t_au > 0
        @test time_au_to_fs(t_au) ≈ t_fs rtol=1e-10

        # Known conversion (1 a.u. time ≈ 24.2 as ≈ 0.0242 fs)
        @test time_au_to_fs(1.0) ≈ 0.024188843 rtol=1e-6

        # Round-trip
        t_test = 100.0
        @test time_au_to_fs(time_fs_to_au(t_test)) ≈ t_test rtol=1e-10
    end

    @testset "Physical Consistency" begin
        # Check that energy-time uncertainty is consistent
        # ΔE·Δt ≥ ℏ/2, in atomic units ℏ=1, so ΔE·Δt ≥ 0.5

        # 1 Ha × 1 a.u. time = 1 (in atomic units)
        E_au = 1.0
        t_au = 1.0
        product_au = E_au * t_au
        @test product_au ≈ 1.0

        # Same product in SI units should be ℏ
        E_eV = au_to_eV(E_au)
        t_fs = time_au_to_fs(t_au)
        # ℏ ≈ 0.658 eV·fs
        product_SI = E_eV * t_fs
        @test product_SI ≈ 0.658 atol=0.01
    end

    @testset "Edge Cases" begin
        # Very small values
        @test au_to_eV(1e-10) ≈ 1e-10 * 27.2114 rtol=1e-6
        @test time_fs_to_au(1e-6) > 0

        # Very large values
        @test intensity_SI_to_au(1e20) > 0
        @test wavelength_nm_to_frequency_au(10000.0) > 0

        # Zero should work for some functions
        @test au_to_eV(0.0) == 0.0
        @test eV_to_au(0.0) == 0.0
    end

end

println("\n✓ PhysicalUnits tests passed")
