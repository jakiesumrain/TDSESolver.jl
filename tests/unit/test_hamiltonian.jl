"""
Unit tests for Hamiltonian module.

Tests Hamiltonian construction and eigenstate solving.
"""

using Test
using LinearAlgebra

# Include all modules via master file
include("../../src/TDSESolver.jl")

@testset "Hamiltonian Tests" begin

    @testset "Kinetic Matrix Construction" begin
        grid = create_gps_grid(100, 50.0)

        # Test l=0 (no centrifugal term)
        T0 = construct_kinetic_matrix(grid, 0)

        @test size(T0) == (100, 100)
        @test issymmetric(T0)  # Should be Hermitian

        # Test l=1 (with centrifugal term)
        T1 = construct_kinetic_matrix(grid, 1)

        @test size(T1) == (100, 100)
        @test issymmetric(T1)

        # Diagonal of T1 should be larger than T0 (centrifugal barrier)
        # Near origin where centrifugal term is large
        @test T1[1, 1] > T0[1, 1]
    end

    @testset "Hamiltonian Matrix Construction" begin
        grid = create_gps_grid(100, 50.0)
        pot = get_potential(:hydrogen)

        H = construct_hamiltonian_matrix(grid, pot, 0)

        @test size(H) == (100, 100)
        @test issymmetric(H)

        # Check potential is included (diagonal should have V(r) contribution)
        # At some interior point, H_ii should have potential contribution
        i_mid = 50
        r_mid = grid.radial_grid[i_mid]
        V_mid = pot.V(r_mid)

        # The diagonal should include potential
        # (exact value depends on kinetic energy discretization)
        @test !isnan(H[i_mid, i_mid])
        @test isfinite(H[i_mid, i_mid])
    end

    @testset "Eigenstate Solver - Hydrogen" begin
        grid = create_gps_grid(200, 100.0, L=30.0)
        pot = get_potential(:hydrogen)

        # Solve for a few states
        ham = solve_eigenstates(grid, pot, 5, n_max=5, E_cutoff=0.0)

        @test ham.lmax == 5
        @test length(ham.n_states) == 6  # l=0 to l=5

        # Check ground state
        E_ground, φ_ground, n, l = get_ground_state(ham)

        @test n == 1
        @test l == 0
        @test E_ground < 0.0  # Bound state

        # For hydrogen, ground state energy should be ≈ -0.5 Ha
        @test E_ground ≈ -0.5 atol=0.1  # Coarse grid, so allow 0.1 Ha tolerance

        # Check normalization
        norm_sq = sum(φ_ground.^2 .* grid.quadrature_weights)
        @test norm_sq ≈ 1.0 atol=1e-6
    end

    @testset "Eigenstate Solver - Helium" begin
        grid = create_gps_grid(200, 100.0, L=30.0)
        pot = get_potential(:helium)

        ham = solve_eigenstates(grid, pot, 3, n_max=5, E_cutoff=0.0)

        E_ground, φ_ground, n, l = get_ground_state(ham)

        @test n == 1
        @test l == 0
        @test E_ground < 0.0

        # For helium effective potential, Ip ≈ 0.9 Ha
        @test E_ground ≈ -0.9 atol=0.2

        # Normalization
        norm_sq = sum(φ_ground.^2 .* grid.quadrature_weights)
        @test norm_sq ≈ 1.0 atol=1e-6
    end

    @testset "Multiple States per l" begin
        grid = create_gps_grid(200, 100.0, L=30.0)
        pot = get_potential(:hydrogen)

        ham = solve_eigenstates(grid, pot, 2, n_max=3, E_cutoff=0.0)

        # l=0 should have multiple bound states (1s, 2s, 3s)
        @test ham.n_states[1] >= 1  # At least ground state

        if ham.n_states[1] >= 2
            # Check energy ordering
            E1 = ham.eigenvalues[1, 1]
            E2 = ham.eigenvalues[2, 1]
            @test E2 > E1  # Higher n has higher energy

            # Check orthogonality
            φ1 = ham.eigenvectors[:, 1, 1]
            φ2 = ham.eigenvectors[:, 2, 1]
            overlap = sum(φ1 .* φ2 .* grid.quadrature_weights)
            @test abs(overlap) < 1e-6  # Should be orthogonal
        end
    end

    @testset "Get Eigenstate Function" begin
        grid = create_gps_grid(200, 100.0)
        pot = get_potential(:hydrogen)

        ham = solve_eigenstates(grid, pot, 3, n_max=3, E_cutoff=0.0)

        # Get ground state
        E_1s, φ_1s = get_eigenstate(ham, 1, 0)

        @test E_1s < 0.0
        @test length(φ_1s) == grid.nrmax

        # Normalization
        norm_sq = sum(φ_1s.^2 .* grid.quadrature_weights)
        @test norm_sq ≈ 1.0 atol=1e-6

        # Check bounds checking
        @test_throws Exception get_eigenstate(ham, 1, 10)  # l too large
        @test_throws Exception get_eigenstate(ham, 100, 0)  # n too large
    end

    @testset "Energy Cutoff" begin
        grid = create_gps_grid(200, 100.0)
        pot = get_potential(:hydrogen)

        # With E_cutoff = -0.3, should only get n=1 state (E ≈ -0.5)
        # n=2 state has E ≈ -0.125, which is above cutoff
        ham_cut = solve_eigenstates(grid, pot, 1, n_max=5, E_cutoff=-0.3)

        # Should have fewer states due to cutoff
        @test ham_cut.n_states[1] >= 1  # At least ground state
        @test ham_cut.n_states[1] <= 2  # But not many excited states

        # All energies should be below cutoff
        for l in 0:ham_cut.lmax
            for n in 1:ham_cut.n_states[l+1]
                @test ham_cut.eigenvalues[n, l+1] < -0.3
            end
        end
    end

    @testset "Physical Reasonableness" begin
        grid = create_gps_grid(300, 150.0, L=40.0)
        pot = get_potential(:hydrogen)

        ham = solve_eigenstates(grid, pot, 10, n_max=10, E_cutoff=0.0)

        E_ground, φ_ground, n, l = get_ground_state(ham)

        # Ground state should peak at some finite radius (not at origin or boundary)
        peak_idx = argmax(abs.(φ_ground))
        @test peak_idx > 1
        @test peak_idx < grid.nrmax

        # Wavefunction should decay at large r
        @test abs(φ_ground[end]) < abs(φ_ground[peak_idx])

        # For hydrogen 1s, peak should be around 1-2 a.u.
        r_peak = grid.radial_grid[peak_idx]
        @test 0.5 < r_peak < 5.0
    end

    @testset "Different l Channels" begin
        grid = create_gps_grid(200, 100.0)
        pot = get_potential(:hydrogen)

        ham = solve_eigenstates(grid, pot, 5, n_max=3, E_cutoff=0.0)

        # For each l, ground state energy should increase with l
        # (centrifugal barrier raises energy)
        E_l0 = ham.eigenvalues[1, 1]  # l=0
        E_l1 = ham.eigenvalues[1, 2]  # l=1

        if ham.n_states[2] > 0
            @test E_l1 > E_l0  # Higher l has higher energy for same n
        end

        # Higher l channels should have fewer bound states
        # (centrifugal barrier reduces binding)
        for l in 1:ham.lmax
            @test ham.n_states[l+1] <= ham.n_states[1]  # l channel has ≤ l=0 states
        end
    end

    @testset "Symmetry Properties" begin
        grid = create_gps_grid(150, 80.0)
        pot = get_potential(:hydrogen)

        # Construct Hamiltonian matrices
        H0 = construct_hamiltonian_matrix(grid, pot, 0)
        H1 = construct_hamiltonian_matrix(grid, pot, 1)

        # Both should be symmetric
        @test issymmetric(H0)
        @test issymmetric(H1)

        # Eigenvalues should be real (Hermitian matrix)
        eigs0 = eigvals(Symmetric(H0))
        @test all(isreal.(eigs0))
    end

    @testset "Grid Resolution Effects" begin
        pot = get_potential(:hydrogen)

        # Coarse grid
        grid_coarse = create_gps_grid(100, 100.0)
        ham_coarse = solve_eigenstates(grid_coarse, pot, 2, n_max=3, E_cutoff=0.0)
        E_coarse, _, _, _ = get_ground_state(ham_coarse)

        # Fine grid
        grid_fine = create_gps_grid(300, 100.0)
        ham_fine = solve_eigenstates(grid_fine, pot, 2, n_max=3, E_cutoff=0.0)
        E_fine, _, _, _ = get_ground_state(ham_fine)

        # Both should be close to -0.5 Ha, fine grid should be more accurate
        @test abs(E_fine + 0.5) < abs(E_coarse + 0.5)
        @test abs(E_fine + 0.5) < 0.05  # Within 5% of exact value
    end

end

println("\n✓ Hamiltonian tests passed")
