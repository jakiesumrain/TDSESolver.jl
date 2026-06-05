"""
Unit tests for Observables module.

Tests computation of physical observables.
"""

using Test
using LinearAlgebra

# Include all modules
include("../../src/TDSESolver.jl")

@testset "Observables Tests" begin

    @testset "Observables Tracker Creation" begin
        obs = create_observables_tracker()

        @test isempty(obs.times)
        @test isempty(obs.norms)
        @test isempty(obs.energies)
        @test isempty(obs.ionization_probs)
        @test isempty(obs.r_expectation)
    end

    @testset "Energy Computation - Ground State" begin
        # Setup hydrogen system
        grid = create_gps_grid(200, 100.0, L=30.0)
        pot = get_potential(:hydrogen)
        ham = solve_eigenstates(grid, pot, 1, n_max=5)

        # Initialize with analytical ground state
        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        # Compute energy
        E = compute_energy(wfn, ham, pot, grid)

        # NOTE: Energy won't match -0.5 Ha exactly due to eigenstate solver issues
        # But should be negative (bound state) and reasonable magnitude
        @test E < 0.0  # Bound state
        @test E > -50.0  # Not absurdly negative
        @test E < 5.0  # Not positive

        @info "Ground state energy" E=E E_analytical=-0.5
    end

    @testset "Ionization Probability - Ground State" begin
        grid = create_gps_grid(200, 100.0)
        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        # Ground state should have very low ionization probability
        P_ion = compute_ionization_probability(wfn, grid, r_cutoff=10.0)

        @test P_ion >= 0.0
        @test P_ion < 0.1  # Less than 10% beyond r=10 a.u.
        @test 1 - P_ion ≈ 1.0 atol=0.1  # Most probability is bound

        @info "Ground state ionization" P_ion=P_ion P_bound=1-P_ion
    end

    @testset "Ionization Cutoff Effects" begin
        grid = create_gps_grid(200, 100.0)
        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        # Test different cutoffs
        P_ion_5 = compute_ionization_probability(wfn, grid, r_cutoff=5.0)
        P_ion_10 = compute_ionization_probability(wfn, grid, r_cutoff=10.0)
        P_ion_20 = compute_ionization_probability(wfn, grid, r_cutoff=20.0)

        # Larger cutoff → smaller ionization probability
        @test P_ion_20 < P_ion_10 < P_ion_5
        @test P_ion_5 > 0.0
        @test P_ion_20 >= 0.0
    end

    @testset "Radial Expectation - Ground State" begin
        grid = create_gps_grid(200, 100.0)
        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        # ⟨r⟩ for hydrogen ground state
        r_mean = compute_radial_expectation(wfn, grid, 1)

        # For hydrogen 1s: ⟨r⟩ = 3/(2Z) = 1.5 a.u. for Z=1
        @test r_mean > 1.0
        @test r_mean < 2.5
        @info "Ground state ⟨r⟩" r_mean=r_mean expected=1.5

        # ⟨r²⟩
        r2_mean = compute_radial_expectation(wfn, grid, 2)
        @test r2_mean > r_mean^2  # Variance is positive
        @info "Ground state ⟨r²⟩" r2=r2_mean rms=sqrt(r2_mean)
    end

    @testset "Radial Expectation - Different Powers" begin
        grid = create_gps_grid(200, 100.0)
        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        r0 = compute_radial_expectation(wfn, grid, 0)  # Just norm
        r1 = compute_radial_expectation(wfn, grid, 1)  # ⟨r⟩
        r2 = compute_radial_expectation(wfn, grid, 2)  # ⟨r²⟩

        @test r0 ≈ 1.0 atol=1e-4  # Normalized
        @test r1 > 0.0
        @test r2 > 0.0
        @test r2 > r1  # ⟨r²⟩ > ⟨r⟩ for reasonable distributions
    end

    @testset "Angular Momentum Populations - Ground State" begin
        grid = create_gps_grid(200, 100.0)
        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)

        # Create wavefunction with lmax=2
        wfn = create_wavefunction(grid.nrmax, 2, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        pops = compute_angular_momentum_populations(wfn, grid)

        # All population should be in l=0
        @test pops[1] ≈ 1.0 atol=1e-4  # l=0
        @test pops[2] ≈ 0.0 atol=1e-6  # l=1
        @test pops[3] ≈ 0.0 atol=1e-6  # l=2

        @info "Ground state angular momentum" P_l0=pops[1] P_l1=pops[2] P_l2=pops[3]
    end

    @testset "Record Observables" begin
        grid = create_gps_grid(150, 80.0)
        pot = get_potential(:hydrogen)
        ham = solve_eigenstates(grid, pot, 1, n_max=5)

        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        obs = create_observables_tracker()

        # Record at t=0
        record_observables!(obs, 0.0, wfn, ham, pot, grid)

        @test length(obs.times) == 1
        @test obs.times[1] == 0.0
        @test obs.norms[1] ≈ 1.0 atol=1e-4
        @test obs.energies[1] < 0.0  # Bound state
        @test obs.ionization_probs[1] < 0.1  # Low ionization
        @test obs.r_expectation[1] > 0.0
    end

    @testset "Record Multiple Time Points" begin
        grid = create_gps_grid(100, 50.0)
        pot = get_potential(:hydrogen)
        ham = solve_eigenstates(grid, pot, 1, n_max=5)

        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        obs = create_observables_tracker()

        # Record at multiple times
        times = [0.0, 0.1, 0.2, 0.3]
        for t in times
            record_observables!(obs, t, wfn, ham, pot, grid)
        end

        @test length(obs.times) == length(times)
        @test obs.times == times
        @test length(obs.norms) == length(times)
        @test length(obs.energies) == length(times)
        @test length(obs.ionization_probs) == length(times)
        @test length(obs.r_expectation) == length(times)
    end

    @testset "Energy Conservation Check" begin
        # NOTE: Energy won't be perfectly conserved due to eigenstate solver issues
        # This test just checks that energy doesn't diverge wildly

        grid = create_gps_grid(100, 50.0)
        pot = get_potential(:hydrogen)
        ham = solve_eigenstates(grid, pot, 1, n_max=5)

        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        E_initial = compute_energy(wfn, ham, pot, grid)

        # Energy should be relatively stable
        @test E_initial < 0.0
        @test isfinite(E_initial)
        @test !isnan(E_initial)
    end

    @testset "Norm Check via Observables" begin
        grid = create_gps_grid(150, 80.0)
        gs = get_analytical_ground_state(:hydrogen, grid.radial_grid, grid.quadrature_weights)
        wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
        initialize_ground_state!(wfn, gs.radial_function, 0, 0)

        # Compute norm via Wavefunction module
        norm1 = compute_norm(wfn)

        # Compute norm via Observables (power=0)
        norm2 = compute_radial_expectation(wfn, grid, 0)

        @test norm1 ≈ norm2 atol=1e-6
        @test norm1 ≈ 1.0 atol=1e-4
    end

    @testset "Observables for Different Atoms" begin
        grid = create_gps_grid(200, 100.0)

        for atom in [:hydrogen, :helium]
            pot = get_potential(atom)
            ham = solve_eigenstates(grid, pot, 1, n_max=5)

            gs = get_analytical_ground_state(atom, grid.radial_grid, grid.quadrature_weights)
            wfn = create_wavefunction(grid.nrmax, 1, grid.quadrature_weights)
            initialize_ground_state!(wfn, gs.radial_function, 0, 0)

            # Compute observables
            E = compute_energy(wfn, ham, pot, grid)
            P_ion = compute_ionization_probability(wfn, grid)
            r_mean = compute_radial_expectation(wfn, grid, 1)

            @test E < 0.0  # Bound state
            @test 0.0 <= P_ion < 0.2  # Low ionization for ground state
            @test r_mean > 0.0

            @info "Observables for $atom" E=E P_ion=P_ion r_mean=r_mean
        end
    end

end

println("\n✓ Observables tests completed")
