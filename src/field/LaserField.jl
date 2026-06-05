"""
    LaserField

Module for laser field calculation following Fortran implementation.

Implements et_sin2 and at_sin2 subroutines from D_inner_out_volkov_3d_with_prob.f90 (lines 1373-1427).

# Field Components
For elliptical/bicircular fields:
- Ex component: oscillates at fundamental frequency ω
- Ey component: oscillates at second harmonic 2ω with phase
- Ez component: typically zero

# Envelope
Trapezoidal envelope with linear ramps over 4 optical cycles (8π/ω).

# Reference
D_inner_out_volkov_3d_with_prob.f90, lines 1373-1427
"""
module LaserField

using Printf

export LaserFieldParameters, create_laser_field, create_custom_laser_field
export compute_electric_field!, compute_vector_potential!
export envelope_function, compute_vector_potential_numerical

"""
    LaserFieldParameters

Container for laser field parameters matching Fortran structure.

# Fields
- `Ex::Float64`: X-component amplitude (a.u.)
- `Ey::Float64`: Y-component amplitude (a.u.)
- `Ez::Float64`: Z-component amplitude (a.u.)
- `omega::Float64`: Fundamental frequency ω (a.u.)
- `tp::Float64`: Pulse duration (a.u.)
- `phase::Float64`: Relative phase between x and y (radians)
- `E_field_x::Vector{Float64}`: Precomputed Ex(t) array
- `E_field_y::Vector{Float64}`: Precomputed Ey(t) array
- `E_field_z::Vector{Float64}`: Precomputed Ez(t) array
- `A_field_x::Vector{Float64}`: Precomputed Ax(t) array
- `A_field_y::Vector{Float64}`: Precomputed Ay(t) array
- `A_field_z::Vector{Float64}`: Precomputed Az(t) array
- `custom_E_func::Union{Nothing, Function}`: Optional custom E(t) function returning [Ex, Ey, Ez]
- `custom_A_func::Union{Nothing, Function}`: Optional custom A(t) function returning [Ax, Ay, Az]
- `is_custom::Bool`: Whether this is a custom field

# Fortran Equivalent
Lines 56, 185: Exx(:), Eyy(:), Ezz(:) arrays
"""
mutable struct LaserFieldParameters
    Ex::Float64
    Ey::Float64
    Ez::Float64
    omega::Float64
    tp::Float64
    phase::Float64
    E_field_x::Vector{Float64}
    E_field_y::Vector{Float64}
    E_field_z::Vector{Float64}
    A_field_x::Vector{Float64}
    A_field_y::Vector{Float64}
    A_field_z::Vector{Float64}
    # Custom field function support
    custom_E_func::Union{Nothing, Function}
    custom_A_func::Union{Nothing, Function}
    is_custom::Bool
end

"""
    envelope_function(t::Float64, tp::Float64, omega::Float64) -> (Float64, Float64)

Compute trapezoidal envelope function and its derivative.

# Arguments
- `t::Float64`: Time (a.u.)
- `tp::Float64`: Pulse duration (a.u.)
- `omega::Float64`: Laser frequency (a.u.)

# Returns
- `ft::Float64`: Envelope function value
- `dft::Float64`: Envelope derivative

# Envelope Shape
- Ramps up linearly over 4 optical cycles (8π/ω)
- Flat top
- Ramps down linearly over 4 optical cycles

# Fortran Reference
Lines 1378-1390, 1405-1417 (identical in both et_sin2 and at_sin2)
"""
function envelope_function(t::Float64, tp::Float64, omega::Float64)
    # Ramp duration: 4 optical cycles = 2 * (2π/ω)
    ramp_duration = 2.0 * 2.0 * π / omega

    t_start = -tp/2.0 + ramp_duration  # End of ramp-up
    t_end = tp/2.0 - ramp_duration     # Start of ramp-down

    if t <= t_start
        # Ramp up
        ft = (t + tp/2.0) / ramp_duration
        dft = 1.0 / ramp_duration
    elseif t > t_start && t <= t_end
        # Flat top
        ft = 1.0
        dft = 0.0
    elseif t > t_end && t <= tp/2.0
        # Ramp down
        ft = -(t - tp/2.0) / ramp_duration
        dft = -1.0 / ramp_duration
    else
        # Outside pulse
        ft = 0.0
        dft = 0.0
    end

    return ft, dft
end

"""
    compute_electric_field(Ex::Float64, Ey::Float64, Ez::Float64,
                          omega::Float64, tp::Float64, t::Float64, phase::Float64) -> Vector{Float64}

Compute electric field components at time t.

# Arguments
- `Ex, Ey, Ez::Float64`: Field amplitudes (a.u.)
- `omega::Float64`: Laser frequency (a.u.)
- `tp::Float64`: Pulse duration (a.u.)
- `t::Float64`: Time (a.u.)
- `phase::Float64`: Relative phase (radians)

# Returns
- `[E_x(t), E_y(t), E_z(t)]`: Electric field vector (a.u.)

# Field Formulas
- E_x(t) = Ex * f(t) * sin(ωt)
- E_y(t) = Ey * f(t) * sin(2ωt + φ)    [Second harmonic!]
- E_z(t) = 0

where f(t) is the trapezoidal envelope function.

# Fortran Reference
Lines 1392-1394 in subroutine et_sin2
"""
function compute_electric_field(Ex::Float64, Ey::Float64, Ez::Float64,
                               omega::Float64, tp::Float64, t::Float64, phase::Float64)
    ft, dft = envelope_function(t, tp, omega)

    # Fortran line 1393-1394
    E = zeros(3)
    E[1] = Ex * ft * sin(omega * t)                    # Fundamental ω
    E[2] = Ey * ft * sin(2.0 * omega * t + phase)      # Second harmonic 2ω
    E[3] = 0.0                                         # Typically zero

    return E
end

"""
    compute_vector_potential(Ex::Float64, Ey::Float64, Ez::Float64,
                            omega::Float64, tp::Float64, t::Float64, phase::Float64) -> Vector{Float64}

Compute vector potential components at time t.

# Arguments
- `Ex, Ey, Ez::Float64`: Field amplitudes (a.u.)
- `omega::Float64`: Laser frequency (a.u.)
- `tp::Float64`: Pulse duration (a.u.)
- `t::Float64`: Time (a.u.)
- `phase::Float64`: Relative phase (radians)

# Returns
- `[A_x(t), A_y(t), A_z(t)]`: Vector potential (a.u.)

# Formulas
A(t) = -∫E(t')dt', integrated with envelope corrections.

# Fortran Reference
Lines 1419-1421 in subroutine at_sin2
"""
function compute_vector_potential(Ex::Float64, Ey::Float64, Ez::Float64,
                                 omega::Float64, tp::Float64, t::Float64, phase::Float64)
    ft, dft = envelope_function(t, tp, omega)

    # Fortran lines 1420-1421
    A = zeros(3)

    # X component (fundamental)
    A[1] = Ex/omega * ft * cos(omega * t) +
           Ex/(omega*omega) * dft * (sin(omega*tp/2.0) - sin(omega*t))

    # Y component (second harmonic)
    omega2 = 2.0 * omega
    A[2] = Ey/omega2 * ft * cos(omega2 * t + phase) +
           Ey/(omega2*omega2) * dft * (sin(omega2*tp/2.0 + phase) - sin(omega2*t + phase))

    # Z component
    A[3] = 0.0

    return A
end

"""
    compute_electric_field!(field::LaserFieldParameters, t::Float64, E_out::Vector{Float64})

In-place version of compute_electric_field.
Supports both built-in and custom field functions.
"""
function compute_electric_field!(field::LaserFieldParameters, t::Float64, E_out::Vector{Float64})
    if field.is_custom && field.custom_E_func !== nothing
        # Use custom field function
        E_custom = field.custom_E_func(t)
        E_out[1] = E_custom[1]
        E_out[2] = E_custom[2]
        E_out[3] = E_custom[3]
    else
        # Use built-in calculation
        E_out .= compute_electric_field(field.Ex, field.Ey, field.Ez,
                                        field.omega, field.tp, t, field.phase)
    end
end

"""
    compute_vector_potential!(field::LaserFieldParameters, t::Float64, A_out::Vector{Float64})

In-place version of compute_vector_potential.
Supports both built-in and custom field functions.
"""
function compute_vector_potential!(field::LaserFieldParameters, t::Float64, A_out::Vector{Float64})
    if field.is_custom && field.custom_A_func !== nothing
        # Use custom vector potential function
        A_custom = field.custom_A_func(t)
        A_out[1] = A_custom[1]
        A_out[2] = A_custom[2]
        A_out[3] = A_custom[3]
    else
        # Use built-in calculation
        A_out .= compute_vector_potential(field.Ex, field.Ey, field.Ez,
                                         field.omega, field.tp, t, field.phase)
    end
end

"""
    create_laser_field(; intensity::Float64=5e14,
                        wavelength::Float64=800.0,
                        duration::Float64=10.0,
                        exrate::Float64=1.0,
                        eyrate::Float64=0.0,
                        ezrate::Float64=0.0,
                        phase::Float64=0.0,
                        dt::Float64=0.1,
                        t_total::Float64=100.0) -> LaserFieldParameters

Create laser field parameters matching Fortran initialization.

# Arguments
- `intensity::Float64=5e14`: Peak intensity (W/cm²)
- `wavelength::Float64=800.0`: Wavelength (nm)
- `duration::Float64=10.0`: Pulse duration FWHM (fs)
- `exrate::Float64=1.0`: X-component amplitude ratio
- `eyrate::Float64=0.0`: Y-component amplitude ratio
- `ezrate::Float64=0.0`: Z-component amplitude ratio
- `phase::Float64=0.0`: Relative phase (radians)
- `dt::Float64=0.1`: Time step (a.u.)
- `t_total::Float64=100.0`: Total propagation time (a.u.)

# Returns
- `LaserFieldParameters`: Field object with precomputed arrays

# Fortran Reference
Lines 133-164, 423-429: Field initialization and array population
"""
function create_laser_field(; intensity::Float64=5e14,
                            wavelength::Float64=800.0,
                            duration::Float64=10.0,
                            exrate::Float64=1.0,
                            eyrate::Float64=0.0,
                            ezrate::Float64=0.0,
                            phase::Float64=0.0,
                            dt::Float64=0.1,
                            t_total::Float64=100.0)

    # Convert to atomic units (matching Fortran lines 133-164)
    intensity_au = intensity / 3.51e16  # W/cm² to a.u.
    E0 = sqrt(intensity_au)

    omega = 0.05695 * 800.0 / wavelength  # a.u.
    tp = duration * 41.341  # fs to a.u.

    # Field component amplitudes (Fortran line 423)
    Ex = E0 * exrate
    Ey = E0 * eyrate
    Ez = E0 * ezrate

    # Precompute field arrays (Fortran lines 425-442)
    n_steps = ceil(Int, t_total / dt) + 1
    E_field_x = zeros(n_steps)
    E_field_y = zeros(n_steps)
    E_field_z = zeros(n_steps)
    A_field_x = zeros(n_steps)
    A_field_y = zeros(n_steps)
    A_field_z = zeros(n_steps)

    # Populate arrays centered at t=0 (Fortran line 427: t = (m - 0.5 - mm*0.5)*dt)
    for m in 1:n_steps
        t = (m - 0.5 - n_steps*0.5) * dt

        E = compute_electric_field(Ex, Ey, Ez, omega, tp, t, phase)
        E_field_x[m] = E[1]
        E_field_y[m] = E[2]
        E_field_z[m] = E[3]

        A = compute_vector_potential(Ex, Ey, Ez, omega, tp, t, phase)
        A_field_x[m] = A[1]
        A_field_y[m] = A[2]
        A_field_z[m] = A[3]
    end

    @info "Created laser field" intensity=intensity wavelength=wavelength duration=duration exrate=exrate eyrate=eyrate phase=phase
    @info "  Field amplitudes" Ex=Ex Ey=Ey Ez=Ez omega=omega tp=tp

    return LaserFieldParameters(Ex, Ey, Ez, omega, tp, phase,
                               E_field_x, E_field_y, E_field_z,
                               A_field_x, A_field_y, A_field_z,
                               nothing, nothing, false)  # Not custom
end

"""
    compute_vector_potential_numerical(E_func::Function, t::Float64, dt::Float64=0.001) -> Vector{Float64}

Compute vector potential A(t) = -∫E(t')dt' using numerical integration.
This is used for custom field functions where analytical A(t) is not available.

# Arguments
- `E_func::Function`: Electric field function E(t) -> [Ex, Ey, Ez]
- `t::Float64`: Time (a.u.)
- `dt::Float64=0.001`: Integration step size (a.u.)

# Returns
- `[Ax(t), Ay(t), Az(t)]`: Vector potential (a.u.)
"""
function compute_vector_potential_numerical(E_func::Function, t::Float64, t_start::Float64=-500.0, dt::Float64=0.01)
    # Integrate from t_start to t using trapezoidal rule
    # A(t) = -∫_{-∞}^{t} E(t') dt'

    A = [0.0, 0.0, 0.0]

    n_steps = ceil(Int, (t - t_start) / dt)
    if n_steps < 1
        return A
    end

    # Trapezoidal integration
    t_current = t_start
    E_prev = E_func(t_current)

    for _ in 1:n_steps
        t_next = min(t_current + dt, t)
        E_next = E_func(t_next)

        # Trapezoidal rule: A += -0.5 * (E_prev + E_next) * dt
        actual_dt = t_next - t_current
        for j in 1:3
            A[j] -= 0.5 * (E_prev[j] + E_next[j]) * actual_dt
        end

        t_current = t_next
        E_prev = E_next
    end

    return A
end

"""
    create_custom_laser_field(E_func::Function;
                              A_func::Union{Nothing, Function}=nothing,
                              dt::Float64=0.1,
                              t_total::Float64=100.0,
                              t_start::Float64=-50.0) -> LaserFieldParameters

Create a laser field from user-provided electric field function.

# Arguments
- `E_func::Function`: Electric field function `E(t) -> [Ex(t), Ey(t), Ez(t)]`
                      where t is in atomic units and E is in atomic units
- `A_func::Union{Nothing, Function}=nothing`: Optional vector potential function `A(t) -> [Ax(t), Ay(t), Az(t)]`
                                              If not provided, A(t) is computed numerically from E(t)
- `dt::Float64=0.1`: Time step for precomputed arrays (a.u.)
- `t_total::Float64=100.0`: Total simulation time (a.u.)
- `t_start::Float64=-50.0`: Start time for numerical integration of A(t)

# Returns
- `LaserFieldParameters`: Field object with custom functions stored

# Example: Linear polarization with Gaussian envelope
```julia
E0 = 0.1  # Field amplitude in a.u.
omega = 0.057  # Frequency in a.u. (800 nm)
sigma = 200.0  # Gaussian width in a.u.

E_func(t) = begin
    envelope = exp(-t^2 / (2*sigma^2))
    [E0 * envelope * sin(omega * t), 0.0, 0.0]
end

field = create_custom_laser_field(E_func, dt=0.1, t_total=500.0)
```

# Example: Circular polarization
```julia
E0 = 0.1
omega = 0.057

E_func(t) = begin
    envelope = t > 0 && t < 400 ? sin(pi*t/400)^2 : 0.0
    [E0 * envelope * sin(omega * t),
     E0 * envelope * cos(omega * t),
     0.0]
end

field = create_custom_laser_field(E_func)
```

# Example: With analytical vector potential
```julia
E0 = 0.1
omega = 0.057

E_func(t) = [E0 * sin(omega * t), 0.0, 0.0]
A_func(t) = [E0/omega * cos(omega * t), 0.0, 0.0]

field = create_custom_laser_field(E_func, A_func=A_func)
```
"""
function create_custom_laser_field(E_func::Function;
                                   A_func::Union{Nothing, Function}=nothing,
                                   dt::Float64=0.1,
                                   t_total::Float64=100.0,
                                   t_start::Float64=-50.0)

    # Create wrapper for numerical A(t) if not provided
    A_func_actual = if A_func !== nothing
        A_func
    else
        t -> compute_vector_potential_numerical(E_func, t, t_start, 0.01)
    end

    # Estimate field parameters from the custom function
    # Sample field at a few points to find approximate amplitude and frequency
    E_sample = E_func(0.0)
    Ex_max = abs(E_sample[1])
    Ey_max = abs(E_sample[2])
    Ez_max = abs(E_sample[3])

    # Scan to find peak amplitude
    for t_test in range(-t_total/4, t_total/4, length=100)
        E_test = E_func(t_test)
        Ex_max = max(Ex_max, abs(E_test[1]))
        Ey_max = max(Ey_max, abs(E_test[2]))
        Ez_max = max(Ez_max, abs(E_test[3]))
    end

    # Set dummy values for omega and tp (not used for custom fields)
    omega_dummy = 0.057  # ~800 nm
    tp_dummy = t_total   # Use total time as pulse duration estimate

    # Precompute field arrays for lookup
    n_steps = ceil(Int, t_total / dt) + 1
    E_field_x = zeros(n_steps)
    E_field_y = zeros(n_steps)
    E_field_z = zeros(n_steps)
    A_field_x = zeros(n_steps)
    A_field_y = zeros(n_steps)
    A_field_z = zeros(n_steps)

    # Populate arrays centered at t=0
    for m in 1:n_steps
        t = (m - 0.5 - n_steps*0.5) * dt

        E = E_func(t)
        E_field_x[m] = E[1]
        E_field_y[m] = E[2]
        E_field_z[m] = E[3]

        A = A_func_actual(t)
        A_field_x[m] = A[1]
        A_field_y[m] = A[2]
        A_field_z[m] = A[3]
    end

    @info "Created custom laser field" Ex_max=Ex_max Ey_max=Ey_max Ez_max=Ez_max t_total=t_total
    @info "  Using $(A_func === nothing ? "numerical" : "analytical") vector potential"

    return LaserFieldParameters(Ex_max, Ey_max, Ez_max, omega_dummy, tp_dummy, 0.0,
                               E_field_x, E_field_y, E_field_z,
                               A_field_x, A_field_y, A_field_z,
                               E_func, A_func_actual, true)  # Custom field
end

"""
    print_field_info(field::LaserFieldParameters)

Print diagnostic information about laser field.
"""
function print_field_info(field::LaserFieldParameters)
    println("Laser Field Parameters:")

    if field.is_custom
        println("  Type: CUSTOM USER-DEFINED FIELD")
        println("  Field amplitudes (estimated max):")
        @printf("    Ex_max = %.6f a.u.\n", field.Ex)
        @printf("    Ey_max = %.6f a.u.\n", field.Ey)
        @printf("    Ez_max = %.6f a.u.\n", field.Ez)
        println("  Custom E(t) function: $(field.custom_E_func !== nothing ? "provided" : "N/A")")
        println("  Custom A(t) function: $(field.custom_A_func !== nothing ? "provided" : "N/A")")
    else
        println("  Type: Built-in ω-2ω field")
        println("  Field amplitudes:")
        @printf("    Ex = %.6f a.u.\n", field.Ex)
        @printf("    Ey = %.6f a.u.\n", field.Ey)
        @printf("    Ez = %.6f a.u.\n", field.Ez)
        @printf("  Frequency: %.6f a.u. (%.2f nm)\n", field.omega, 0.05695*800/field.omega)
        @printf("  Period: %.2f a.u. (%.2f as)\n", 2π/field.omega, 2π/field.omega*24.2)
        @printf("  Pulse duration: %.2f a.u. (%.2f fs)\n", field.tp, field.tp/41.341)
        @printf("  Phase: %.3f rad (%.1f°)\n", field.phase, field.phase*180/π)

        # Polarization type
        if abs(field.Ey) < 1e-10 && abs(field.Ez) < 1e-10
            println("  Polarization: Linear (X)")
        elseif abs(field.Ex) < 1e-10 && abs(field.Ez) < 1e-10
            println("  Polarization: Linear (Y)")
        elseif abs(field.Ez) > 1e-10
            println("  Polarization: 3D")
        else
            ellipticity = abs(field.Ey / field.Ex)
            if ellipticity ≈ 1.0 && abs(field.phase - π/2) < 0.1
                println("  Polarization: Circular")
            else
                @printf("  Polarization: Elliptical (ε = %.3f, φ = %.2f°)\n",
                       ellipticity, field.phase*180/π)
            end
        end
    end

    println("  Precomputed arrays: $(length(field.E_field_x)) time steps")
end

end  # module LaserField
