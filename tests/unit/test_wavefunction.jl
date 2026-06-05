"""
Unit tests for Wavefunction module.

Tests wavefunction representation, initialization, and operations.
"""

using Test

# Include the modules
include("../../src/grid/GPSGrid.jl")
using .GPSGrid

include("../../src/wavefunction/Wavefunction.jl")
using .Wavefunction

@testset "Wavefunction Tests" begin

    @testset "Wavefunction Creation" begin
        grid = create_gps_grid(100, 50.0)
        wfn = create_wavefunction(grid.nrmax, 20, grid.quadrature_weights)

        @test wfn.nrmax == 100
        @test wfn.lmax == 20
        @test size(wfn.g) == (100, 41, 21)  # (nrmax, 2*lmax+1, lmax+1)
        @test length(wfn.grid_weights) == 100
    end

    @testset "Wavefunction Indexing" begin
        grid = create_gps_grid(50, 30.0)
        wfn = create_wavefunction(grid.nrmax, 10, grid.quadrature_weights)

        # Test array dimensions match expected indexing
        @test size(wfn.g, 1) == wfn.nrmax
        @test size(wfn.g, 2) == 2*wfn.lmax + 1  # m from -lmax to +lmax
        @test size(wfn.g, 3) == wfn.lmax + 1     # l from 0 to lmax
    end

    @testset "Ground State Initialization" begin
        grid = create_gps_grid(100, 50.0)
        wfn = create_wavefunction(grid.nrmax, 10, grid.quadrature_weights)

        # Create a simple normalized ground state (Gaussian-like)
        ground_state = exp.(-grid.radial_grid ./ 2.0)
        ground_state ./= sqrt(sum(ground_state.^2 .* grid.quadrature_weights))

        initialize_ground_state!(wfn, ground_state, 1, 0)

        # Check normalization
        norm = compute_norm(wfn)
        @test norm ≈ 1.0 atol=1e-10

        # Check only (l=0, m=0) component is populated
        g_00 = get_component(wfn, 0, 0)
        @test maximum(abs.(g_00 .- ground_state)) < 1e-10

        # Check other components are zero
        g_10 = get_component(wfn, 1, 0)
        @test maximum(abs.(g_10)) < 1e-14
    end

    @testset "Norm Computation" begin
        grid = create_gps_grid(100, 50.0)
        wfn = create_wavefunction(grid.nrmax, 5, grid.quadrature_weights)

        # Initially zero
        @test compute_norm(wfn) ≈ 0.0

        # Set some components
        radial = exp.(-grid.radial_grid ./ 3.0)
        set_component!(wfn, 0, 0, radial .+ 0.0im)

        norm_before = compute_norm(wfn)
        @test norm_before > 0.0

        # Normalize
        wfn.g ./= norm_before
        norm_after = compute_norm(wfn)
        @test norm_after ≈ 1.0 atol=1e-10
    end

    @testset "Radial Density" begin
        grid = create_gps_grid(100, 50.0)
        wfn = create_wavefunction(grid.nrmax, 5, grid.quadrature_weights)

        # Set ground state
        ground_state = exp.(-grid.radial_grid ./ 2.0)
        initialize_ground_state!(wfn, ground_state, 1, 0)

        density = compute_radial_density(wfn)

        @test length(density) == wfn.nrmax
        @test all(density .>= 0.0)  # Density is non-negative

        # For normalized ground state, integral of density with weights should be 1
        integral = sum(density .* grid.quadrature_weights)
        @test integral ≈ 1.0 atol=1e-10
    end

    @testset "Component Get/Set" begin
        grid = create_gps_grid(50, 30.0)
        wfn = create_wavefunction(grid.nrmax, 10, grid.quadrature_weights)

        # Set a component
        test_radial = rand(ComplexF64, grid.nrmax)
        set_component!(wfn, 3, 2, test_radial)

        # Get it back
        retrieved = get_component(wfn, 3, 2)
        @test retrieved ≈ test_radial

        # Check other components are still zero
        other = get_component(wfn, 3, 1)
        @test all(abs.(other) .< 1e-14)
    end

    @testset "Component Bounds Checking" begin
        grid = create_gps_grid(50, 30.0)
        wfn = create_wavefunction(grid.nrmax, 5, grid.quadrature_weights)

        # l out of range
        @test_throws Exception get_component(wfn, 10, 0)
        @test_throws Exception get_component(wfn, -1, 0)

        # m out of range for given l
        @test_throws Exception get_component(wfn, 2, 3)  # |m| > l
        @test_throws Exception get_component(wfn, 2, -3)

        # Valid ranges should work
        @test isa(get_component(wfn, 3, -2), Vector{ComplexF64})
        @test isa(get_component(wfn, 3, 2), Vector{ComplexF64})
    end

    @testset "Copy Wavefunction" begin
        grid = create_gps_grid(50, 30.0)
        wfn1 = create_wavefunction(grid.nrmax, 5, grid.quadrature_weights)

        # Set some values
        ground_state = exp.(-grid.radial_grid)
        initialize_ground_state!(wfn1, ground_state)

        # Copy
        wfn2 = copy_wavefunction(wfn1)

        # Check they're equal
        @test wfn2.g ≈ wfn1.g
        @test wfn2.nrmax == wfn1.nrmax
        @test wfn2.lmax == wfn1.lmax

        # Modify copy
        wfn2.g .*= 2.0

        # Original should be unchanged
        @test !(wfn2.g ≈ wfn1.g)
        @test compute_norm(wfn2) ≈ 2.0 * compute_norm(wfn1)
    end

    @testset "Multiple Components" begin
        grid = create_gps_grid(100, 50.0)
        wfn = create_wavefunction(grid.nrmax, 10, grid.quadrature_weights)

        # Set multiple components
        r1 = exp.(-grid.radial_grid ./ 2.0) .+ 0.0im
        r2 = exp.(-grid.radial_grid ./ 3.0) .+ 0.0im
        r3 = exp.(-grid.radial_grid ./ 4.0) .+ 0.0im

        set_component!(wfn, 0, 0, r1)
        set_component!(wfn, 1, 0, r2)
        set_component!(wfn, 2, 1, r3)

        # Check all three are set correctly
        @test get_component(wfn, 0, 0) ≈ r1
        @test get_component(wfn, 1, 0) ≈ r2
        @test get_component(wfn, 2, 1) ≈ r3

        # Norm should be sum of individual norms
        norm_total = compute_norm(wfn)
        norm1 = sqrt(sum(abs2.(r1) .* grid.quadrature_weights))
        norm2 = sqrt(sum(abs2.(r2) .* grid.quadrature_weights))
        norm3 = sqrt(sum(abs2.(r3) .* grid.quadrature_weights))

        @test norm_total ≈ sqrt(norm1^2 + norm2^2 + norm3^2) atol=1e-10
    end

    @testset "Input Validation" begin
        grid = create_gps_grid(50, 30.0)

        # nrmax too small
        @test_throws Exception create_wavefunction(5, 10, grid.quadrature_weights)

        # Negative lmax
        @test_throws Exception create_wavefunction(50, -1, grid.quadrature_weights)

        # Mismatched grid weights
        wrong_weights = zeros(Float64, 30)
        @test_throws Exception create_wavefunction(50, 10, wrong_weights)
    end

    @testset "Physical Reasonableness" begin
        grid = create_gps_grid(200, 100.0)
        wfn = create_wavefunction(grid.nrmax, 20, grid.quadrature_weights)

        # Hydrogen ground state: ψ(r) ∝ exp(-r)
        psi_1s = @. 2.0 * exp(-grid.radial_grid)  # Normalized for ∫|ψ|²r²dr=1
        # For GPS grid, need to normalize with GPS weights
        psi_1s ./= sqrt(sum(psi_1s.^2 .* grid.quadrature_weights))

        initialize_ground_state!(wfn, psi_1s, 1, 0)

        # Should be normalized
        @test compute_norm(wfn) ≈ 1.0 atol=1e-10

        # Density should decay exponentially
        density = compute_radial_density(wfn)
        @test density[1] > density[end]  # Decays outward

        # Peak should be near origin
        max_idx = argmax(density)
        @test max_idx < wfn.nrmax ÷ 4  # Peak in first quarter of grid
    end

    @testset "Memory Usage" begin
        grid = create_gps_grid(400, 150.0)
        lmax_large = 50
        wfn = create_wavefunction(grid.nrmax, lmax_large, grid.quadrature_weights)

        # Check memory size is reasonable
        mem_mb = sizeof(wfn.g) / 1024^2
        expected_mb = (400 * (2*50+1) * (50+1) * 16) / 1024^2  # ComplexF64 = 16 bytes

        @test mem_mb ≈ expected_mb atol=0.1
        @test mem_mb < 100.0  # Should be less than 100 MB for these parameters
    end

end

println("\n✓ Wavefunction tests passed")
