"""
Advanced Custom Potential Module Example

This file demonstrates how to define custom atomic potentials using Julia modules
for advanced users who need full control over potential functions, including
analytic derivatives and parameter optimization.

# Usage

1. Define your potential module (see MyCustomPotentials below)
2. Load in Julia REPL or script:
   ```julia
   include("examples/custom_potential_module.jl")
   using .MyCustomPotentials
   pot = MyCustomPotentials.yukawa_potential()
   ```
3. Pass to TDSE solver via Simulation.jl

# Benefits over String Expressions

- Analytic derivatives (more accurate than finite differences)
- Parameter optimization capabilities
- Complex potentials with multiple components
- Better performance (no runtime parsing)
- Type stability
"""

module MyCustomPotentials

export yukawa_potential, soft_core_potential, gaussian_potential

# Import potential wrapper type
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/hamiltonian/Potential.jl")
using .Potential: PotentialFunction

"""
    yukawa_potential(Z::Float64=2.0, λ::Float64=2.0) -> PotentialFunction

Yukawa (screened Coulomb) potential:
    V(r) = -Z * exp(-λr) / r

# Parameters
- `Z::Float64=2.0`: Effective charge
- `λ::Float64=2.0`: Screening parameter (inverse length scale)

# Physical Interpretation
Models screening effects in atomic/nuclear systems. λ controls screening range:
- λ → 0: Pure Coulomb (-Z/r)
- λ → ∞: Short-range interaction

# Example
```julia
pot = yukawa_potential(2.0, 2.0)  # Z=2, λ=2
V_at_5au = pot.V(5.0)
```
"""
function yukawa_potential(Z::Float64=2.0, λ::Float64=2.0)
    # Potential function
    V = function(r)
        return -Z * exp(-λ * r) / r
    end

    # Analytic derivative: dV/dr = Z * exp(-λr) * (1 + λr) / r²
    dVdr = function(r)
        exp_term = exp(-λ * r)
        return Z * exp_term * (1.0 + λ * r) / (r * r)
    end

    description = "Yukawa potential: V(r) = -$Z * exp(-$(λ)r) / r"
    return PotentialFunction(V, dVdr, :yukawa, description)
end

"""
    soft_core_potential(Z::Float64=1.0, a::Float64=0.5) -> PotentialFunction

Soft-core Coulomb potential:
    V(r) = -Z / sqrt(r² + a²)

# Parameters
- `Z::Float64=1.0`: Effective charge
- `a::Float64=0.5`: Softening parameter (prevents singularity at r=0)

# Physical Interpretation
Regularized Coulomb potential avoiding r=0 singularity. Common in molecular dynamics
and strong-field physics. a controls regularization strength:
- a → 0: Approaches Coulomb
- a large: Smoother potential, lower ionization probability

# Example
```julia
pot = soft_core_potential(1.0, 0.5)  # H-like with a=0.5
```
"""
function soft_core_potential(Z::Float64=1.0, a::Float64=0.5)
    # Potential function
    V = function(r)
        return -Z / sqrt(r*r + a*a)
    end

    # Analytic derivative: dV/dr = Z * r / (r² + a²)^(3/2)
    dVdr = function(r)
        denom = (r*r + a*a)^1.5
        return Z * r / denom
    end

    description = "Soft-core potential: V(r) = -$Z / sqrt(r² + $(a)²)"
    return PotentialFunction(V, dVdr, :soft_core, description)
end

"""
    gaussian_potential(V0::Float64=10.0, σ::Float64=2.0) -> PotentialFunction

Gaussian well potential:
    V(r) = -V0 * exp(-r²/(2σ²))

# Parameters
- `V0::Float64=10.0`: Depth of potential well (Ha)
- `σ::Float64=2.0`: Width of Gaussian (a.u.)

# Physical Interpretation
Smooth finite-range attractive potential. Models effective interactions in
ultracold atoms, quantum dots, or simplified nuclear potentials.

# Example
```julia
pot = gaussian_potential(10.0, 2.0)  # 10 Ha deep, 2 a.u. wide
```
"""
function gaussian_potential(V0::Float64=10.0, σ::Float64=2.0)
    # Potential function
    V = function(r)
        return -V0 * exp(-r*r / (2.0 * σ*σ))
    end

    # Analytic derivative: dV/dr = -V0 * exp(-r²/(2σ²)) * (-r/σ²)
    #                              = V0 * r/σ² * exp(-r²/(2σ²))
    dVdr = function(r)
        exp_term = exp(-r*r / (2.0 * σ*σ))
        return V0 * r / (σ*σ) * exp_term
    end

    description = "Gaussian potential: V(r) = -$V0 * exp(-r²/(2*$(σ)²))"
    return PotentialFunction(V, dVdr, :gaussian, description)
end

"""
    print_potential_comparison()

Demonstrate potential comparison at various r values.
"""
function print_potential_comparison()
    println("Custom Potential Comparison")
    println("="^80)

    pots = [
        ("Yukawa (Z=2, λ=2)", yukawa_potential(2.0, 2.0)),
        ("Soft-core (Z=1, a=0.5)", soft_core_potential(1.0, 0.5)),
        ("Gaussian (V0=10, σ=2)", gaussian_potential(10.0, 2.0))
    ]

    r_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

    for (name, pot) in pots
        println("\n$name:")
        println("  r (a.u.)  |  V(r) (Ha)  |  dV/dr (Ha/a.u.)")
        println("  " * "-"^50)
        for r in r_values
            V_r = pot.V(r)
            dV_r = pot.dVdr(r)
            @printf("  %6.2f    |  %10.6f  |  %10.6f\n", r, V_r, dV_r)
        end
    end

    println("\n" * "="^80)
end

end  # module MyCustomPotentials

# Example usage (uncomment to run)
# using .MyCustomPotentials
# pot = yukawa_potential()
# println("Yukawa potential at r=5: ", pot.V(5.0), " Ha")
# print_potential_comparison()
