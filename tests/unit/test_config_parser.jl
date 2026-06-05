"""
Unit tests for ConfigParser module.

Tests TOML parsing, validation, and configuration object creation.
"""

using Test
using TOML

# Include the module
include("../../src/io/ConfigParser.jl")
using .ConfigParser

@testset "ConfigParser Tests" begin

    @testset "Valid Configuration Loading" begin
        # Create a temporary valid config file
        config_content = """
        [atom]
        type = "helium"

        [laser]
        wavelength = 800.0
        intensity = 5e14
        pulse_duration = 10.0
        polarization = [1.0, 0.0, 0.0]

        [calculation]
        type = "ionization"
        output_file = "test_results.h5"
        verbosity = "normal"
        """

        temp_file = tempname() * ".toml"
        write(temp_file, config_content)

        try
            config = load_config(temp_file)

            @test config.atom_type == :helium
            @test config.wavelength_nm == 800.0
            @test config.intensity_W_cm2 == 5e14
            @test config.pulse_duration_fs == 10.0
            @test config.polarization ≈ [1.0, 0.0, 0.0]
            @test config.calculation_type == :ionization
            @test config.output_file == "test_results.h5"
            @test config.verbosity == :normal

            # Check atomic unit conversions happened
            @test config.frequency_au > 0
            @test config.intensity_au > 0
            @test config.pulse_duration_au > 0

        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Polarization Normalization" begin
        config_content = """
        [atom]
        type = "hydrogen"

        [laser]
        wavelength = 800.0
        intensity = 1e14
        pulse_duration = 10.0
        polarization = [2.0, 0.0, 0.0]

        [calculation]
        type = "ionization"
        output_file = "test.h5"
        """

        temp_file = tempname() * ".toml"
        write(temp_file, config_content)

        try
            config = load_config(temp_file)
            # Should be normalized
            @test sum(config.polarization.^2) ≈ 1.0 atol=1e-10
            @test config.polarization[1] ≈ 1.0
        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Invalid Configurations" begin
        # Missing required section
        @test_throws ArgumentError begin
            temp_file = tempname() * ".toml"
            write(temp_file, "[laser]\nwavelength = 800.0")
            try
                load_config(temp_file)
            finally
                rm(temp_file, force=true)
            end
        end

        # Negative wavelength
        @test_throws ArgumentError begin
            temp_file = tempname() * ".toml"
            config_content = """
            [atom]
            type = "helium"
            [laser]
            wavelength = -800.0
            intensity = 1e14
            pulse_duration = 10.0
            [calculation]
            type = "ionization"
            output_file = "test.h5"
            """
            write(temp_file, config_content)
            try
                load_config(temp_file)
            finally
                rm(temp_file, force=true)
            end
        end

        # Invalid atom type
        @test_throws ArgumentError begin
            temp_file = tempname() * ".toml"
            config_content = """
            [atom]
            type = "krypton"
            [laser]
            wavelength = 800.0
            intensity = 1e14
            pulse_duration = 10.0
            [calculation]
            type = "ionization"
            output_file = "test.h5"
            """
            write(temp_file, config_content)
            try
                load_config(temp_file)
            finally
                rm(temp_file, force=true)
            end
        end

        # Invalid calculation type
        @test_throws ArgumentError begin
            temp_file = tempname() * ".toml"
            config_content = """
            [atom]
            type = "helium"
            [laser]
            wavelength = 800.0
            intensity = 1e14
            pulse_duration = 10.0
            [calculation]
            type = "scattering"
            output_file = "test.h5"
            """
            write(temp_file, config_content)
            try
                load_config(temp_file)
            finally
                rm(temp_file, force=true)
            end
        end
    end

    @testset "Custom Potential Configuration" begin
        config_content = """
        [atom]
        type = "custom"
        potential = "-1/sqrt(r^2 + 0.5)"
        charge = 1.0

        [laser]
        wavelength = 800.0
        intensity = 5e14
        pulse_duration = 10.0
        polarization = [1.0, 0.0, 0.0]

        [calculation]
        type = "ionization"
        output_file = "custom.h5"
        """

        temp_file = tempname() * ".toml"
        write(temp_file, config_content)

        try
            config = load_config(temp_file)
            @test config.atom_type == :custom
            @test config.custom_potential == "-1/sqrt(r^2 + 0.5)"
            @test config.charge == 1.0
        finally
            rm(temp_file, force=true)
        end
    end

    @testset "Parameter Scan Configuration" begin
        config_content = """
        [atom]
        type = "helium"

        [laser]
        wavelength = 800.0
        intensity = 5e14
        pulse_duration = 10.0
        polarization = [1.0, 0.0, 0.0]

        [calculation]
        type = "ionization"
        output_file = "scan.h5"

        [scan]
        parameter = "laser.intensity"
        range = "1e14:1e14:5e14"
        """

        temp_file = tempname() * ".toml"
        write(temp_file, config_content)

        try
            config = load_config(temp_file)
            @test !isnothing(config.scan)
            @test config.scan.parameter == :laser_intensity
            @test length(config.scan.range) == 5  # 1e14, 2e14, 3e14, 4e14, 5e14
        finally
            rm(temp_file, force=true)
        end
    end

    @testset "File Not Found" begin
        @test_throws ArgumentError load_config("nonexistent_file.toml")
    end

end

println("\n✓ ConfigParser tests passed")
