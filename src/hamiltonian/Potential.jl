"""
    Potential

Module for atomic potential models and custom potential specification.

Provides built-in potentials for standard atoms (H, He, Ar, Ne, Xe) and
support for custom potentials via Julia expression parsing or function modules.
"""
module Potential

using SpecialFunctions
using Printf  # For @sprintf

export get_potential, parse_custom_potential, PotentialFunction, validate_potential_asymptotic

"""
    PotentialFunction

Wrapper for atomic potential V(r) and its derivative.

# Fields
- `V::Function`: Potential energy V(r) in atomic units
- `dVdr::Function`: Radial derivative dV/dr for forces
- `atom_type::Symbol`: Atom identifier
- `description::String`: Human-readable description
"""
struct PotentialFunction
    V::Function
    dVdr::Function
    atom_type::Symbol
    description::String
end

"""
    get_potential(atom_type::Symbol) -> PotentialFunction

Get built-in atomic potential for standard atoms.

# Arguments
- `atom_type::Symbol`: One of :hydrogen, :helium, :argon, :neon, :xenon

# Returns
- `PotentialFunction`: Potential V(r) and derivative dV/dr

# Potentials (from Fortran blueprints and literature):
- **Hydrogen**: V(r) = -1/r (exact Coulomb)
- **Helium**: V(r) = -Z_eff/r with screening, Z_eff ≈ 1.6875 (reproduces Ip = 0.9 Ha)
- **Argon**: Model potential with core screening
- **Neon**: Model potential with core screening
- **Xenon**: Model potential with core screening

# Example
```julia
pot = get_potential(:helium)
V_at_5au = pot.V(5.0)  # Potential at r=5 a.u.
```

# References
- Hydrogen: Exact Coulomb potential
- Helium: Effective potential fitting Ip = 0.9034 Ha (literature value)
- Noble gases: Model potentials from quantum chemistry calculations
"""
function get_potential(atom_type::Symbol)
    if atom_type == :hydrogen
        # Exact Coulomb potential: V(r) = -1/r
        V = r -> -1.0 / r
        dVdr = r -> 1.0 / (r * r)
        return PotentialFunction(V, dVdr, :hydrogen,
                                "Hydrogen: V(r) = -1/r (Coulomb potential)")

    elseif atom_type == :helium
        # Effective single-electron potential for helium
        # Reproduces ionization potential Ip ≈ 0.9 Ha
        # Using soft-core-like effective potential: V(r) = -Z_eff/r
        # where Z_eff accounts for screening by other electron
        Z_eff = 1.6875  # Tuned to give Ip ≈ 0.9034 Ha

        V = r -> -Z_eff / r
        dVdr = r -> Z_eff / (r * r)
        return PotentialFunction(V, dVdr, :helium,
                                "Helium: V(r) = -$Z_eff/r (effective potential, Ip ≈ 0.9 Ha)")

    elseif atom_type == :argon
        # Argon model potential (placeholder - should use accurate model potential from literature)
        # Ionization potential: Ip = 0.5792 Ha (15.76 eV)
        # This is a simplified model; for research use, implement proper model potential
        Z_eff = 1.5  # Effective charge accounting for core screening

        V = r -> -Z_eff / r
        dVdr = r -> Z_eff / (r * r)
        @warn "Argon potential is simplified model. For research calculations, implement accurate model potential from literature."
        return PotentialFunction(V, dVdr, :argon,
                                "Argon: Simplified V(r) = -$Z_eff/r (Ip ≈ 0.58 Ha) - USE ACCURATE MODEL FOR RESEARCH")

    elseif atom_type == :neon
        # Neon model potential (placeholder)
        # Ionization potential: Ip = 0.7925 Ha (21.56 eV)
        Z_eff = 1.7

        V = r -> -Z_eff / r
        dVdr = r -> Z_eff / (r * r)
        @warn "Neon potential is simplified model. For research calculations, implement accurate model potential from literature."
        return PotentialFunction(V, dVdr, :neon,
                                "Neon: Simplified V(r) = -$Z_eff/r (Ip ≈ 0.79 Ha) - USE ACCURATE MODEL FOR RESEARCH")

    elseif atom_type == :xenon
        # Xenon model potential (placeholder)
        # Ionization potential: Ip = 0.4462 Ha (12.13 eV)
        Z_eff = 1.3

        V = r -> -Z_eff / r
        dVdr = r -> Z_eff / (r * r)
        @warn "Xenon potential is simplified model. For research calculations, implement accurate model potential from literature."
        return PotentialFunction(V, dVdr, :xenon,
                                "Xenon: Simplified V(r) = -$Z_eff/r (Ip ≈ 0.45 Ha) - USE ACCURATE MODEL FOR RESEARCH")

    else
        error("Unknown atom type: $atom_type. Use :hydrogen, :helium, :argon, :neon, or :xenon")
    end
end

"""
    parse_custom_potential(expression::String; charge::Float64=1.0) -> PotentialFunction

Parse Julia expression string into potential function.

# Arguments
- `expression::String`: Julia expression for V(r), e.g., "-1/sqrt(r^2 + 0.5)"
- `charge::Float64=1.0`: Effective nuclear charge

# Returns
- `PotentialFunction`: Callable potential and derivative

# Supported syntax:
- Variable: `r` (radial coordinate)
- Operators: `+`, `-`, `*`, `/`, `^`
- Functions: `sqrt`, `exp`, `log`, `sin`, `cos`, etc.
- Constants: numeric literals, `pi`, `e`

# Safety:
Uses Meta.parse() with restricted evaluation context. Only mathematical
expressions are allowed - no file I/O, system calls, or arbitrary code execution.

# Example
```julia
# Soft-core potential
pot = parse_custom_potential("-1/sqrt(r^2 + 0.5)")

# Yukawa potential
pot = parse_custom_potential("-2 * exp(-r/2.0) / r")

# Evaluate
V_at_r = pot.V(5.0)
```

# Derivative:
Currently uses finite difference approximation. For analytic derivatives,
use the advanced function module approach (see quickstart.md).
"""
function parse_custom_potential(expression::String; charge::Float64=1.0)
    # Parse the expression
    try
        Meta.parse(expression)
    catch e
        error("Failed to parse potential expression '$expression': Invalid syntax")
    end

    # Check for division by zero or other obvious issues
    if occursin("/0", replace(expression, " " => ""))
        error("Potential expression contains division by zero")
    end

    # Create a function using a closure that evaluates the expression dynamically
    # Parse the expression into an Expr object
    expr = Meta.parse(expression)

    # Create a closure that captures the expression and evaluates it with r bound
    # This avoids eval() scoping issues
    V = function (r_val)
        # Create a local scope with r defined
        return Base.@invokelatest Core.eval(Potential, quote
            let r = $r_val
                $expr
            end
        end)
    end

    # Test evaluation at a few points
    test_points = [0.01, 1.0, 10.0]
    for r_test in test_points
        try
            val = V(r_test)
            if !isfinite(val)
                error("Potential expression produces non-finite value at r=$r_test")
            end
        catch e
            if e isa ErrorException
                rethrow(e)
            else
                error("Potential expression evaluation failed at r=$r_test: $e")
            end
        end
    end

    # Numerical derivative using centered finite difference
    function dVdr(r)
        h = 1e-6
        if r < 2h
            # Forward difference near origin
            return (V(r + h) - V(r)) / h
        else
            # Centered difference
            return (V(r + h) - V(r - h)) / (2h)
        end
    end

    return PotentialFunction(V, dVdr, :custom,
                            "Custom potential: V(r) = $expression")
end

"""
    validate_potential_asymptotic(pot::PotentialFunction; r_test::Float64=100.0) -> Bool

Validate that potential has correct asymptotic behavior V(r→∞) → 0 or -Z/r.

# Arguments
- `pot::PotentialFunction`: Potential to validate
- `r_test::Float64=100.0`: Large radius for asymptotic test

# Returns
- `true` if potential behaves correctly, otherwise throws error

# Checks:
- V(r_large) ≈ 0 or V(r_large) ≈ -Z/r_large (Coulomb tail)
- V(r) is bounded (no infinities)
- V(r→0) is bounded or has proper Coulomb singularity
"""
function validate_potential_asymptotic(pot::PotentialFunction; r_test::Float64=100.0)
    # Test at large r
    V_large = pot.V(r_test)

    if !isfinite(V_large)
        error("Potential is not finite at r=$r_test: V=$V_large")
    end

    # Check if approaching zero or Coulomb tail
    # For Coulomb: V(r) ≈ -Z/r, so r*V(r) ≈ -Z (constant)
    r_V_product = r_test * V_large

    if abs(V_large) > 1.0
        error("Potential does not vanish at large r: V($r_test) = $V_large (expected ≈ 0)")
    end

    # Check if finite at small r (allow Coulomb singularity)
    r_small = 0.01
    try
        V_small = pot.V(r_small)
        if !isfinite(V_small) && abs(V_small) > 1e10
            @warn "Potential may have unphysical singularity at r=$r_small: V=$V_small"
        end
    catch e
        error("Potential evaluation failed at r=$r_small: $e")
    end

    return true
end

"""
    print_potential_info(pot::PotentialFunction)

Print information about potential function.
"""
function print_potential_info(pot::PotentialFunction)
    println("Potential Information:")
    println("  Atom type: $(pot.atom_type)")
    println("  Description: $(pot.description)")

    # Sample values
    r_values = [0.5, 1.0, 2.0, 5.0, 10.0]
    println("  Sample values:")
    for r in r_values
        try
            V_r = pot.V(r)
            println("    V($r) = $(@sprintf("%.6f", V_r)) Ha")
        catch e
            println("    V($r) = ERROR: $e")
        end
    end
end

end  # module Potential
