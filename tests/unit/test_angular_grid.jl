"""
Unit tests for AngularGrid module.

Tests angular grid creation, solid angle integration, and quadrature properties.
"""

using Test

# Include the module
include("../../src/grid/AngularGrid.jl")
using .AngularGrid

@testset "AngularGrid Tests" begin

    @testset "Uniform Grid Creation" begin
        grid = create_angular_grid(180, 60, 10, method=:uniform)

        @test grid.nthmax == 180
        @test grid.nphimax == 60
        @test length(grid.theta_grid) == 180
        @test length(grid.phi_grid) == 60
        @test length(grid.theta_weights) == 180
        @test length(grid.phi_weights) == 60
    end

    @testset "Gauss-Legendre Grid Creation" begin
        grid = create_angular_grid(90, 50, 10, method=:gauss)

        @test grid.nthmax == 90
        @test grid.nphimax == 50
        @test length(grid.theta_grid) == 90
        @test length(grid.phi_grid) == 50
    end

    @testset "Theta Range [0, π]" begin
        # Uniform method
        grid_uniform = create_angular_grid(100, 50, 10, method=:uniform)
        @test grid_uniform.theta_grid[1] ≈ 0.0 atol=1e-10
        @test grid_uniform.theta_grid[end] ≈ π atol=1e-10
        @test all(0.0 .<= grid_uniform.theta_grid .<= π)

        # Gauss method (should cover [0, π] but not include exact endpoints)
        grid_gauss = create_angular_grid(100, 50, 10, method=:gauss)
        @test minimum(grid_gauss.theta_grid) >= 0.0
        @test maximum(grid_gauss.theta_grid) <= π
        @test all(0.0 .<= grid_gauss.theta_grid .<= π)
    end

    @testset "Phi Range [0, 2π)" begin
        grid = create_angular_grid(90, 72, 10, method=:uniform)

        # Should cover [0, 2π) without repeating endpoint
        @test grid.phi_grid[1] ≈ 0.0 atol=1e-10
        @test grid.phi_grid[end] < 2π
        @test all(0.0 .<= grid.phi_grid .< 2π)

        # Spacing should be uniform
        Δφ = 2π / 72
        for i in 2:length(grid.phi_grid)
            @test grid.phi_grid[i] - grid.phi_grid[i-1] ≈ Δφ rtol=1e-6
        end
    end

    @testset "Solid Angle Integration (Should Equal 4π)" begin
        # Uniform method
        grid_uniform = create_angular_grid(180, 120, 10, method=:uniform)
        total_uniform = 0.0
        for i in 1:grid_uniform.nthmax
            for j in 1:grid_uniform.nphimax
                total_uniform += compute_solid_angle_element(grid_uniform, i, j)
            end
        end
        @test total_uniform ≈ 4π rtol=0.01  # Within 1% for uniform grid

        # Gauss method (should be more accurate)
        grid_gauss = create_angular_grid(90, 60, 10, method=:gauss)
        total_gauss = 0.0
        for i in 1:grid_gauss.nthmax
            for j in 1:grid_gauss.nphimax
                total_gauss += compute_solid_angle_element(grid_gauss, i, j)
            end
        end
        @test total_gauss ≈ 4π rtol=1e-6  # Gauss-Legendre should be very accurate
    end

    @testset "Solid Angle Element Formula" begin
        # Test both grid types

        # Uniform grid: dΩ = sin(θ) * dθ * dφ
        grid_uniform = create_angular_grid(90, 60, 10, method=:uniform)
        for i in 1:grid_uniform.nthmax
            for j in 1:grid_uniform.nphimax
                θ = grid_uniform.theta_grid[i]
                w_θ = grid_uniform.theta_weights[i]
                w_φ = grid_uniform.phi_weights[j]

                dΩ_computed = compute_solid_angle_element(grid_uniform, i, j)
                dΩ_expected = sin(θ) * w_θ * w_φ

                @test dΩ_computed ≈ dΩ_expected atol=1e-12
            end
        end

        # Gauss grid: dΩ = w_x * dφ (sin(θ) implicit in transformation)
        grid_gauss = create_angular_grid(90, 60, 10, method=:gauss)
        for i in 1:grid_gauss.nthmax
            for j in 1:grid_gauss.nphimax
                w_θ = grid_gauss.theta_weights[i]
                w_φ = grid_gauss.phi_weights[j]

                dΩ_computed = compute_solid_angle_element(grid_gauss, i, j)
                dΩ_expected = w_θ * w_φ

                @test dΩ_computed ≈ dΩ_expected atol=1e-12
            end
        end
    end

    @testset "Weight Positivity" begin
        grid = create_angular_grid(100, 80, 10, method=:gauss)

        @test all(grid.theta_weights .> 0)
        @test all(grid.phi_weights .> 0)
    end

    @testset "Integration of Spherical Harmonics" begin
        # Test integral of Y₀⁰ = 1/√(4π) over sphere
        # ∫∫ |Y₀⁰|² dΩ = ∫∫ (1/4π) dΩ = 1
        grid = create_angular_grid(90, 60, 10, method=:gauss)

        Y00_squared = 1.0 / (4π)
        integral = 0.0

        for i in 1:grid.nthmax
            for j in 1:grid.nphimax
                dΩ = compute_solid_angle_element(grid, i, j)
                integral += Y00_squared * dΩ
            end
        end

        @test integral ≈ 1.0 rtol=1e-6
    end

    @testset "Theta Monotonicity" begin
        # Uniform grid: θ increasing from 0 to π
        grid_uniform = create_angular_grid(100, 50, 10, method=:uniform)
        for i in 2:length(grid_uniform.theta_grid)
            @test grid_uniform.theta_grid[i] > grid_uniform.theta_grid[i-1]
        end

        # Gauss grid with arccos: θ decreasing (since arccos is decreasing)
        # x goes from -1 to +1, so arccos(x) goes from π to 0
        grid_gauss = create_angular_grid(100, 50, 10, method=:gauss)
        for i in 2:length(grid_gauss.theta_grid)
            @test grid_gauss.theta_grid[i] < grid_gauss.theta_grid[i-1]
        end

        # But all θ values should still be in [0, π]
        @test all(0.0 .<= grid_gauss.theta_grid .<= π)
    end

    @testset "Phi Monotonicity" begin
        grid = create_angular_grid(90, 72, 10, method=:uniform)

        # φ should be monotonically increasing
        for i in 2:length(grid.phi_grid)
            @test grid.phi_grid[i] > grid.phi_grid[i-1]
        end
    end

    @testset "Method Comparison" begin
        nthmax = 90
        nphimax = 60
        lmax = 10

        grid_uniform = create_angular_grid(nthmax, nphimax, lmax, method=:uniform)
        grid_gauss = create_angular_grid(nthmax, nphimax, lmax, method=:gauss)

        # Both should have same dimensions
        @test grid_uniform.nthmax == grid_gauss.nthmax
        @test grid_uniform.nphimax == grid_gauss.nphimax

        # Gauss method should have better integration accuracy
        total_uniform = 0.0
        total_gauss = 0.0

        for i in 1:nthmax
            for j in 1:nphimax
                total_uniform += compute_solid_angle_element(grid_uniform, i, j)
                total_gauss += compute_solid_angle_element(grid_gauss, i, j)
            end
        end

        # Gauss should be closer to 4π
        error_uniform = abs(total_uniform - 4π)
        error_gauss = abs(total_gauss - 4π)
        @test error_gauss < error_uniform
    end

    @testset "Input Validation" begin
        # nthmax too small
        @test_throws Exception create_angular_grid(1, 60, 10)

        # nphimax too small
        @test_throws Exception create_angular_grid(90, 1, 10)

        # Invalid method
        @test_throws Exception create_angular_grid(90, 60, 10, method=:invalid)
    end

    @testset "Small Grid" begin
        # Minimal grid should still work
        grid = create_angular_grid(2, 2, 10, method=:gauss)

        @test grid.nthmax == 2
        @test grid.nphimax == 2
        @test length(grid.theta_grid) == 2
        @test length(grid.phi_grid) == 2
    end

    @testset "Large Grid" begin
        # Larger grid should work efficiently
        grid = create_angular_grid(360, 180, 10, method=:gauss)

        @test grid.nthmax == 360
        @test grid.nphimax == 180

        # Total solid angle should still be 4π
        total = 0.0
        for i in 1:grid.nthmax
            for j in 1:grid.nphimax
                total += compute_solid_angle_element(grid, i, j)
            end
        end
        @test total ≈ 4π rtol=1e-8
    end

    @testset "Reproducibility" begin
        # Same parameters should give same grid
        grid1 = create_angular_grid(90, 60, 10, method=:gauss)
        grid2 = create_angular_grid(90, 60, 10, method=:gauss)

        @test grid1.theta_grid ≈ grid2.theta_grid
        @test grid1.phi_grid ≈ grid2.phi_grid
        @test grid1.theta_weights ≈ grid2.theta_weights
        @test grid1.phi_weights ≈ grid2.phi_weights
    end

    @testset "Integration of cos(θ)" begin
        # ∫∫ cos(θ) sin(θ) dθ dφ = ∫₀^(2π) dφ ∫₀^π cos(θ) sin(θ) dθ
        # = 2π * [sin²(θ)/2]₀^π = 2π * 0 = 0
        grid = create_angular_grid(180, 120, 10, method=:gauss)

        integral = 0.0
        for i in 1:grid.nthmax
            for j in 1:grid.nphimax
                θ = grid.theta_grid[i]
                dΩ = compute_solid_angle_element(grid, i, j)
                integral += cos(θ) * dΩ
            end
        end

        @test abs(integral) < 1e-10  # Should be zero
    end

end

println("\n✓ AngularGrid tests passed")
