"""
Unit tests for Potential module.

Tests built-in atomic potentials and custom potential parsing.
"""

using Test

# Include the module
include("../../src/hamiltonian/Potential.jl")
using .Potential

@testset "Potential Tests" begin

    @testset "Hydrogen Potential" begin
        pot = get_potential(:hydrogen)

        @test pot.atom_type == :hydrogen
        @test occursin("Coulomb", pot.description)

        # Test V(r) = -1/r
        @test pot.V(1.0) ≈ -1.0
        @test pot.V(2.0) ≈ -0.5
        @test pot.V(10.0) ≈ -0.1

        # Test derivative dV/dr = 1/r²
        @test pot.dVdr(1.0) ≈ 1.0
        @test pot.dVdr(2.0) ≈ 0.25
        @test pot.dVdr(5.0) ≈ 0.04

        # Asymptotic behavior
        @test abs(pot.V(100.0)) < 0.011  # Allow small numerical tolerance
    end

    @testset "Helium Potential" begin
        pot = get_potential(:helium)

        @test pot.atom_type == :helium

        # Should give reasonable ionization potential
        # V(r) at large r should approach 0
        @test abs(pot.V(100.0)) < 0.02

        # Should be more attractive than hydrogen at same r
        @test pot.V(1.0) < -1.0
    end

    @testset "Noble Gas Potentials" begin
        for atom in [:argon, :neon, :xenon]
            pot = get_potential(atom)
            @test pot.atom_type == atom

            # Basic sanity checks
            @test pot.V(1.0) < 0  # Attractive
            @test pot.V(10.0) > pot.V(1.0)  # Less attractive at larger r
            @test abs(pot.V(100.0)) < 0.05  # Approaches zero
        end
    end

    @testset "Custom Potential Parsing" begin
        # Soft-core potential
        expr = "-1/sqrt(r^2 + 0.5)"
        pot = parse_custom_potential(expr)

        @test pot.atom_type == :custom
        @test occursin(expr, pot.description)

        # Test evaluation
        V_r1 = pot.V(1.0)
        @test isfinite(V_r1)
        @test V_r1 < 0  # Should be attractive

        # Test at origin (should be finite for soft-core)
        V_r0 = pot.V(0.01)
        @test isfinite(V_r0)

        # Test derivative
        dV = pot.dVdr(1.0)
        @test isfinite(dV)
    end

    @testset "Custom Potential Expressions" begin
        # Yukawa potential
        pot_yukawa = parse_custom_potential("-2 * exp(-r/2.0) / r")
        @test isfinite(pot_yukawa.V(1.0))

        # Gaussian potential
        pot_gauss = parse_custom_potential("-exp(-r^2)")
        @test pot_gauss.V(0.0) ≈ -1.0
        @test abs(pot_gauss.V(5.0)) < 0.01

        # Simple harmonic oscillator
        pot_ho = parse_custom_potential("0.5 * r^2")
        @test pot_ho.V(0.0) ≈ 0.0
        @test pot_ho.V(2.0) ≈ 2.0
    end

    @testset "Invalid Expressions" begin
        # Invalid syntax
        @test_throws Exception parse_custom_potential("1/0")

        # Non-mathematical expression (should fail safely)
        @test_throws Exception parse_custom_potential("system('ls')")
    end

    @testset "Asymptotic Validation" begin
        # Test hydrogen
        pot_h = get_potential(:hydrogen)
        @test validate_potential_asymptotic(pot_h, r_test=100.0)

        # Test custom potential
        pot_custom = parse_custom_potential("-1/r")
        @test validate_potential_asymptotic(pot_custom, r_test=100.0)
    end

    @testset "Potential Energy Values" begin
        # Hydrogen at specific radii
        pot_h = get_potential(:hydrogen)

        # At Bohr radius (r=1 a.u.), V = -1 Ha
        @test pot_h.V(1.0) ≈ -1.0

        # At r=0.5 a.u., V = -2 Ha
        @test pot_h.V(0.5) ≈ -2.0 atol=1e-10

        # Helium should have effective Z > 1
        pot_he = get_potential(:helium)
        @test abs(pot_he.V(1.0)) > 1.0  # More attractive than hydrogen
    end

    @testset "Derivative Accuracy" begin
        pot = get_potential(:hydrogen)

        # Numerical check of derivative
        r = 2.0
        h = 1e-6
        numerical_deriv = (pot.V(r + h) - pot.V(r - h)) / (2h)

        @test pot.dVdr(r) ≈ numerical_deriv rtol=1e-4
    end

end

println("\n✓ Potential tests passed")
