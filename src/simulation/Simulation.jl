"""
    Simulation

Main orchestrator module for running TDSE simulations.

Provides high-level interface to set up and run time-dependent Schrödinger equation
simulations for atoms in strong laser fields.

# Workflow
1. Define simulation parameters (grid, atom, laser)
2. Initialize system (ground state, propagator)
3. Run time propagation with laser field
4. Record observables during propagation
5. Return results for analysis

# Example
```julia
# Define parameters
params = SimulationParams(
    atom = :hydrogen,
    nrmax = 200,
    rmax = 100.0,
    lmax = 2,
    dt = 0.1,
    t_total = 100.0,
    laser_wavelength = 800.0,  # nm
    laser_intensity = 1e14,    # W/cm²
    laser_duration = 50.0      # a.u.
)

# Run simulation
results = run_simulation(params)

# Access results
plot(results.observables.times, results.observables.ionization_probs)
```
"""
module Simulation

using Printf
using LinearAlgebra

# Import all required modules
using ..PhysicalUnits: wavelength_nm_to_frequency_au, intensity_SI_to_au
using ..GPSGrid: GPSGridData, create_gps_grid
using ..AngularGrid: AngularGridData, create_angular_grid
using ..Potential: PotentialFunction, get_potential, parse_custom_potential
using ..Wavefunction: WavefunctionData, create_wavefunction, initialize_ground_state!, compute_norm
using ..Hamiltonian: HamiltonianData, solve_eigenstates, get_ground_state
using ..Hamiltonian: load_eigenstates_from_hdf5, save_eigenstates_to_hdf5, print_hamiltonian_info
using ..Propagator: PropagatorData, create_propagator, propagate_step!
using ..Observables: ObservablesData, create_observables_tracker, record_observables!
using ..Observables: print_observables_summary
using ..HHG: compute_dipole_eigenstate!, compute_dipole_acceleration!, compute_dipole_length!
using ..LaserField: compute_electric_field, LaserFieldParameters, compute_electric_field!
using ..RegionSplit: RegionSplitterData, create_region_splitter, split_wavefunction!
using ..VolkovProjection: MomentumGrid, create_momentum_grid, VolkovProjectorData, create_volkov_projector, project_volkov!

export SimulationParams, SimulationResults
export run_simulation, create_default_params, compute_laser_field
export print_simulation_params
export generate_eigenstates  # For standalone eigenstate generation

"""
    SimulationParams

Container for all simulation parameters.

# Eigenstate Mode Control
**Three operational modes via `eigenstate_mode` switch:**

1. **Generate only** (`eigenstate_mode = :generate_only`):
   - Compute eigenstates and save to HDF5 file
   - Exit without running simulation
   - Use for pre-computing eigenstates on a cluster

2. **Runtime** (`eigenstate_mode = :runtime`, default):
   - Compute eigenstates at runtime
   - Build S-matrix and run simulation
   - No HDF5 file needed

3. **Load from file** (`eigenstate_mode = :load_from_file`):
   - Load pre-computed eigenstates from HDF5 file
   - Build S-matrix from loaded data
   - Fastest startup for repeated simulations

# Analysis Control Switches (Independent)
**Four operational modes via two Boolean flags:**

1. **Ionization only** (fastest):
   - `enable_momentum_analysis = false, enable_hhg = false`

2. **Ionization + Photoelectron momentum** (ATI spectra):
   - `enable_momentum_analysis = true, enable_hhg = false`

3. **Ionization + HHG spectrum** (harmonic generation):
   - `enable_momentum_analysis = false, enable_hhg = true`

4. **Full analysis** (all physics):
   - `enable_momentum_analysis = true, enable_hhg = true`

See `docs/ANALYSIS_MODES.md` for detailed usage examples.

# Grid Parameters
- `nrmax::Int`: Number of radial grid points
- `rmax::Float64`: Maximum radial distance (a.u.)
- `L::Float64`: GPS grid mapping parameter (α = 2L/rmax computed automatically)
- `lmax::Int`: Maximum angular momentum

# Atomic Parameters
- `atom::Symbol`: Atom type (:hydrogen, :helium, etc.)

# Time Propagation
- `dt::Float64`: Time step (a.u.)
- `t_total::Float64`: Total simulation time (a.u.)
- `obs_interval::Int`: Observables recording interval (steps)

# Laser Field Parameters
- `laser_wavelength::Float64`: Wavelength (nm)
- `laser_intensity::Float64`: Peak intensity (W/cm²)
- `laser_duration::Float64`: Pulse duration FWHM (a.u.)
- `laser_cep::Float64`: Carrier-envelope phase (radians)
- `laser_polarization::Vector{Float64}`: Polarization vector [Ex, Ey, Ez] (simple mode)

# Bicircular Field Parameters (optional, requires use_bicircular=true)
- `use_bicircular::Bool`: Use bicircular field (ω + 2ω) from LaserField module
- `exrate::Float64`: X-component amplitude ratio (ω component)
- `eyrate::Float64`: Y-component amplitude ratio (2ω component)
- `ezrate::Float64`: Z-component amplitude ratio
- `phase_offset::Float64`: Relative phase between ω and 2ω (radians)

# Custom Laser Field (optional, overrides other laser params when set)
- `custom_laser_field::Union{LaserFieldParameters, Nothing}`: User-defined laser field
  Created with `create_custom_laser_field(E_func)` where E_func(t) returns [Ex, Ey, Ez]

# Numerical Parameters
- `energy_cutoff::Float64`: Energy cutoff for S-matrix (Ha)
- `ionization_cutoff::Float64`: Radial cutoff for ionization (a.u.)

# Eigenstate Mode Parameters
- `eigenstate_mode::Symbol`: Mode for eigenstate handling (default: :runtime)
  - `:generate_only`: Compute and save to HDF5, no simulation
  - `:runtime`: Compute at runtime (default)
  - `:load_from_file`: Load from HDF5 file
- `eigenstate_file::String`: HDF5 file path for eigenstates (default: "eigenstates.h5")

# Region Splitting and Momentum Analysis Parameters
- `enable_momentum_analysis::Bool`: Enable region splitting and momentum distribution calculation
- `msplit::Int`: Split every msplit time steps (Fortran default: 50)
- `R_c::Float64`: Region splitting radius (a.u., Fortran default: 100.0)
- `delta_split::Float64`: Splitting smoothness parameter (a.u., Fortran default: 5.0)
- `p_max::Float64`: Maximum momentum for momentum grid (a.u., default: 2.0)
- `n_p::Int`: Number of momentum grid points (default: 30)
- `n_theta::Int`: Number of θ grid points for momentum space (default: 20)
- `n_phi::Int`: Number of φ grid points for momentum space (default: 20)

# HHG Spectrum Calculation Parameters
- `enable_hhg::Bool`: Enable HHG dipole moment tracking (default: false)
- `hhg_method::Symbol`: Dipole calculation method (default: :eigenstate)
  - `:eigenstate`: Eigenstate expansion (Fortran rescatteing+hhg-he.f90, lines 891-899, active)
  - `:acceleration`: Acceleration form (Fortran lines 905-909, commented)
  - `:length`: Length form
- `hhg_record_interval::Int`: Record dipole every N steps (default: 1)
"""
Base.@kwdef struct SimulationParams
    # Grid
    nrmax::Int = 200
    rmax::Float64 = 100.0
    L::Float64 = 30.0
    # Note: α is computed automatically as α = 2L/rmax in create_gps_grid
    lmax::Int = 2

    # Atom
    atom::Symbol = :hydrogen
    custom_potential::Union{String, Function, Nothing} = nothing  # Custom potential for atom=:custom

    # Time
    dt::Float64 = 0.1
    t_total::Float64 = 100.0
    obs_interval::Int = 10

    # Laser (defaults: 800 nm, 1e14 W/cm², 50 a.u. duration)
    laser_wavelength::Float64 = 800.0  # nm
    laser_intensity::Float64 = 1.0e14  # W/cm²
    laser_duration::Float64 = 50.0     # a.u. FWHM
    laser_cep::Float64 = 0.0           # radians
    laser_polarization::Vector{Float64} = [0.0, 0.0, 1.0]  # z-polarized (simple mode)

    # Bicircular field (optional - requires use_bicircular=true)
    use_bicircular::Bool = false       # Enable bicircular field (ω + 2ω)
    exrate::Float64 = 1.0              # X amplitude ratio (ω component)
    eyrate::Float64 = 0.0              # Y amplitude ratio (2ω component)
    ezrate::Float64 = 0.0              # Z amplitude ratio
    phase_offset::Float64 = 0.0        # Relative phase (radians)

    # Custom laser field (optional - overrides other laser parameters when set)
    # Use create_custom_laser_field() to create, then pass here
    custom_laser_field::Union{LaserFieldParameters, Nothing} = nothing

    # Numerical
    energy_cutoff::Float64 = 50.0
    ionization_cutoff::Float64 = 10.0

    # ============================================================================
    # EIGENSTATE MODE CONTROL
    # ============================================================================
    # Controls how eigenstates are obtained for S-matrix construction
    eigenstate_mode::Symbol = :runtime           # :generate_only, :runtime (default), :load_from_file
    eigenstate_file::String = "eigenstates.h5"   # HDF5 file path for eigenstates

    # ============================================================================
    # ANALYSIS CONTROL SWITCHES (Independent - enable any combination)
    # ============================================================================

    # Region splitting and momentum analysis (Fortran lines 162-164, 786-794, 909-923)
    # Set to TRUE for photoelectron momentum distributions (ATI spectra)
    enable_momentum_analysis::Bool = false     # Enable region splitting and Volkov projection
    msplit::Int = 50                           # Split every msplit time steps
    R_c::Float64 = 100.0                       # Splitting radius (a.u.)
    delta_split::Float64 = 5.0                 # Smoothness parameter (a.u.)
    p_max::Float64 = 2.0                       # Maximum momentum (a.u.)
    n_p::Int = 30                              # Momentum grid points
    n_theta::Int = 20                          # θ grid points for momentum space
    n_phi::Int = 20                            # φ grid points for momentum space

    # HHG spectrum calculation (Fortran: rescatteing+hhg-he.f90)
    # Set to TRUE for high harmonic generation spectrum
    enable_hhg::Bool = false                   # Enable HHG dipole moment tracking
    hhg_method::Symbol = :eigenstate           # :eigenstate (Fortran 891-899), :acceleration (905-909), :length
    hhg_record_interval::Int = 1               # Record dipole every N steps (1 = every step for accurate spectrum)
end

"""
    SimulationResults

Container for simulation results.

# Fields
- `params::SimulationParams`: Input parameters
- `grid::GPSGridData`: Radial grid
- `hamiltonian::HamiltonianData`: Hamiltonian with eigenstates
- `wavefunction::WavefunctionData`: Final wavefunction
- `observables::ObservablesData`: Time series of observables
- `wall_time::Float64`: Computation time (seconds)
- `momentum_grid::Union{MomentumGrid,Nothing}`: Momentum space grid (if momentum analysis enabled)
- `psai_p::Union{Array{ComplexF64,3},Nothing}`: Momentum space wavefunction ψ(φ,θ,p) (if momentum analysis enabled)
"""
struct SimulationResults
    params::SimulationParams
    grid::GPSGridData
    hamiltonian::HamiltonianData
    wavefunction::WavefunctionData
    observables::ObservablesData
    wall_time::Float64
    momentum_grid::Union{MomentumGrid,Nothing}
    psai_p::Union{Array{ComplexF64,3},Nothing}
end

"""
    create_default_params(; kwargs...) -> SimulationParams

Create simulation parameters with sensible defaults.

Allows overriding specific parameters via keyword arguments.

# Example
```julia
params = create_default_params(
    atom = :helium,
    laser_intensity = 5e13,
    t_total = 200.0
)
```
"""
function create_default_params(; kwargs...)
    return SimulationParams(; kwargs...)
end

"""
    compute_laser_field(t::Float64, params::SimulationParams) -> Vector{Float64}

Compute electric field vector at time t.

# Mode Selection
- If `params.use_bicircular = false`: Uses simple Gaussian envelope (default)
- If `params.use_bicircular = true`: Uses bicircular field (ω + 2ω) from LaserField module

# Simple Mode (Gaussian Envelope)
Uses Gaussian envelope:
    E(t) = E₀ exp(-2ln(2) (t-t₀)²/τ²) cos(ω(t-t₀) + φ_CEP)

where:
- E₀: Peak field amplitude
- τ: FWHM duration
- ω: Carrier frequency
- φ_CEP: Carrier-envelope phase
- t₀: Pulse center time = t_total/2

Direction set by `laser_polarization` vector.

# Custom Mode
When `custom_laser_field` is provided, uses the custom E(t) function directly.

# Bicircular Mode (ω + 2ω)
Uses trapezoidal envelope with bicircular field from LaserField module:
- Ex(t) = Ex * f(t) * sin(ωt)              [Fundamental ω]
- Ey(t) = Ey * f(t) * sin(2ωt + φ)         [Second harmonic 2ω]
- Ez(t) = 0

where:
- Ex, Ey, Ez: Field amplitudes from exrate, eyrate, ezrate
- f(t): Trapezoidal envelope (ramps over 4 optical cycles)
- φ: phase_offset (relative phase between ω and 2ω)

Pulse centered at t = 0.

# Returns
- `Vector{Float64}`: [Ex, Ey, Ez] in atomic units

# References
- Simple mode: Standard Gaussian pulse
- Bicircular mode: D_inner_out_volkov_3d_with_prob.f90 lines 1373-1427
"""
function compute_laser_field(t::Float64, params::SimulationParams)
    # ===== Custom Mode (user-provided E(t) function) =====
    if params.custom_laser_field !== nothing
        E_out = zeros(3)
        compute_electric_field!(params.custom_laser_field, t, E_out)
        return E_out
    elseif params.use_bicircular
        # ===== Bicircular Mode (LaserField module) =====
        # Convert parameters to atomic units
        ω = wavelength_nm_to_frequency_au(params.laser_wavelength)
        I_au = intensity_SI_to_au(params.laser_intensity)
        E0 = sqrt(I_au)  # Field amplitude from intensity

        # Field component amplitudes (using exrate, eyrate, ezrate)
        Ex = E0 * params.exrate
        Ey = E0 * params.eyrate
        Ez = E0 * params.ezrate

        # Pulse duration and phase
        tp = params.laser_duration
        phase = params.phase_offset

        # Compute field using LaserField module
        # Pulse is centered at t=0 (matching Fortran convention)
        E_field = compute_electric_field(Ex, Ey, Ez, ω, tp, t, phase)

        return E_field
    else
        # ===== Simple Mode (Gaussian envelope) =====
        # Convert laser parameters to atomic units
        ω = wavelength_nm_to_frequency_au(params.laser_wavelength)
        I_au = intensity_SI_to_au(params.laser_intensity)
        E0 = sqrt(2.0 * I_au)  # Electric field amplitude from intensity
        τ = params.laser_duration

        # Pulse centered at t_total/2
        t0 = params.t_total / 2.0

        # Gaussian envelope
        envelope = exp(-2.0 * log(2.0) * (t - t0)^2 / τ^2)

        # Carrier with CEP
        carrier = cos(ω * (t - t0) + params.laser_cep)

        # Total field amplitude
        E_amplitude = E0 * envelope * carrier

        # Apply polarization
        E_field = E_amplitude .* params.laser_polarization

        return E_field
    end
end

"""
    run_simulation(params::SimulationParams; verbose::Bool=true) -> SimulationResults

Run complete TDSE simulation.

# Algorithm
1. Initialize grid, potential, Hamiltonian
2. Initialize wavefunction with analytical ground state
3. Create propagator with pre-computed S-matrices
4. Time propagation loop:
   - Compute laser field E(t)
   - Propagate wavefunction: ψ(t+Δt) = Û(Δt) ψ(t)
   - Record observables at specified intervals
5. Return results

# Arguments
- `params::SimulationParams`: Simulation parameters
- `verbose::Bool=true`: Print progress information

# Returns
- `SimulationResults`: Complete simulation results

# Example
```julia
params = create_default_params(atom=:hydrogen, t_total=100.0)
results = run_simulation(params)
println("Final ionization: \$(results.observables.ionization_probs[end])")
```
"""
function run_simulation(params::SimulationParams; verbose::Bool=true)
    start_time = time()

    if verbose
        @info "Starting TDSE simulation" atom=params.atom nrmax=params.nrmax dt=params.dt t_total=params.t_total
    end

    # ===== Step 1: Initialize grid =====
    if verbose
        @info "Creating radial grid"
    end
    grid = create_gps_grid(params.nrmax, params.rmax, L=params.L)

    # ===== Step 2: Initialize potential and Hamiltonian =====
    if verbose
        @info "Setting up Hamiltonian" eigenstate_mode=params.eigenstate_mode
    end

    # Get potential function
    if params.atom == :custom
        if isnothing(params.custom_potential)
            error("Custom atom type requires custom_potential specification")
        end

        if typeof(params.custom_potential) == String
            pot = parse_custom_potential(params.custom_potential)
            if verbose
                @info "Parsed custom potential from expression" expression=params.custom_potential
            end
        elseif typeof(params.custom_potential) <: Function
            # Assume it's already a PotentialFunction
            pot = params.custom_potential
            if verbose
                @info "Using provided custom potential function"
            end
        else
            error("custom_potential must be a String expression or PotentialFunction")
        end
    else
        pot = get_potential(params.atom)
    end

    # Handle eigenstates based on eigenstate_mode
    ham = nothing

    if params.eigenstate_mode == :generate_only
        # Mode 1: Compute eigenstates and save to HDF5, then return early
        if verbose
            @info "Mode: generate_only - Computing eigenstates and saving to HDF5"
        end

        # Use n_max=nrmax and E_cutoff=Inf to get complete basis for unitary S-matrix
        ham = solve_eigenstates(grid, pot, params.lmax,
                               n_max=grid.nrmax, E_cutoff=Inf)

        # Save to HDF5
        save_eigenstates_to_hdf5(params.eigenstate_file, ham)

        if verbose
            print_hamiltonian_info(ham)
            @info "Eigenstates saved to $(params.eigenstate_file). Exiting."
        end

        # Return early with minimal results (no simulation run)
        wall_time = time() - start_time
        wfn = create_wavefunction(grid.nrmax, params.lmax, grid.quadrature_weights)
        obs = create_observables_tracker(enable_hhg=params.enable_hhg)

        return SimulationResults(params, grid, ham, wfn, obs, wall_time, nothing, nothing)

    elseif params.eigenstate_mode == :load_from_file
        # Mode 3: Load pre-computed eigenstates from HDF5
        if verbose
            @info "Mode: load_from_file - Loading eigenstates from $(params.eigenstate_file)"
        end

        ham = load_eigenstates_from_hdf5(params.eigenstate_file, grid, pot)

        if verbose
            @info "Eigenstates loaded successfully"
        end

    else
        # Mode 2: Runtime computation (default)
        if verbose
            @info "Mode: runtime - Computing eigenstates at runtime"
        end

        # Use n_max=nrmax and E_cutoff=Inf to get complete basis for unitary S-matrix
        ham = solve_eigenstates(grid, pot, params.lmax,
                               n_max=grid.nrmax, E_cutoff=Inf)
    end

    # ===== Step 3: Initialize wavefunction with ground state =====
    if verbose
        @info "Initializing wavefunction with ground state"
    end

    # Get numerical ground state from Hamiltonian (all atoms use this now)
    E_ground, φ_ground, n_gs, l_gs = get_ground_state(ham)
    if verbose
        @info "Using numerical ground state" E=E_ground n=n_gs l=l_gs
    end

    wfn = create_wavefunction(grid.nrmax, params.lmax, grid.quadrature_weights)
    initialize_ground_state!(wfn, φ_ground, 0, 0)  # l=0, m=0

    # Verify normalization
    norm_initial = compute_norm(wfn)
    if verbose
        @info "Initial state prepared" norm=norm_initial energy=E_ground
    end

    # ===== Step 4: Create propagator =====
    if verbose
        @info "Creating propagator with S-matrices"
    end
    ang_grid = create_angular_grid(90, 60, params.lmax)  # θ × φ grid with lmax for transformations
    prop = create_propagator(ham, ang_grid, params.dt)

    # ===== Step 5: Initialize observables tracker =====
    obs = create_observables_tracker(enable_hhg=params.enable_hhg)

    # Record initial observables
    record_observables!(obs, 0.0, wfn, ham, pot, grid,
                       r_cutoff=params.ionization_cutoff)

    # ===== Step 5b: Initialize momentum analysis (if enabled) =====
    splitter = nothing
    wfn_outer = nothing
    volkov_projector = nothing
    A_accumulated = zeros(Float64, 3)  # Accumulated vector potential

    if params.enable_momentum_analysis
        if verbose
            @info "Initializing region splitting and momentum analysis" R_c=params.R_c delta=params.delta_split msplit=params.msplit
        end

        # Create region splitter
        splitter = create_region_splitter(grid, params.R_c, params.delta_split)

        # Create outer wavefunction container
        wfn_outer = create_wavefunction(grid.nrmax, params.lmax, grid.quadrature_weights)

        # Create momentum grid (using keyword arguments)
        momentum_grid = create_momentum_grid(p_max=params.p_max, n_p=params.n_p,
                                             n_theta=params.n_theta, n_phi=params.n_phi)

        # Create Volkov projector (now requires angular_grid)
        volkov_projector = create_volkov_projector(grid, params.lmax, momentum_grid, ang_grid)

        if verbose
            @info "Momentum analysis ready" p_max=params.p_max n_p=params.n_p n_theta=params.n_theta n_phi=params.n_phi
        end
    end

    # ===== Step 6: Time propagation loop =====
    n_steps = Int(ceil(params.t_total / params.dt))

    if verbose
        @info "Starting time propagation" n_steps=n_steps obs_interval=params.obs_interval
    end

    t = 0.0
    ω = wavelength_nm_to_frequency_au(params.laser_wavelength)  # For Volkov phase

    for step in 1:n_steps
        # Compute laser field at current time
        E_field = compute_laser_field(t, params)

        # Propagate one time step
        propagate_step!(wfn, prop, E_field, t)

        # Advance time
        t += params.dt

        # ===== Region splitting and Volkov projection (if enabled) =====
        if params.enable_momentum_analysis && (step % params.msplit == 0)
            # Split wavefunction into outer region
            split_wavefunction!(wfn_outer, wfn, splitter)

            # Accumulate vector potential: A(t) = -∫ E(t') dt'
            A_accumulated .+= -E_field .* params.dt

            # Project to momentum space with Volkov state (now requires dt)
            project_volkov!(volkov_projector, wfn_outer, A_accumulated, t, params.dt)
        end

        # ===== HHG dipole moment calculation (if enabled) =====
        if params.enable_hhg && (step % params.hhg_record_interval == 0)
            if params.hhg_method == :eigenstate
                # Method 1: Eigenstate expansion (Fortran lines 891-899, active)
                compute_dipole_eigenstate!(obs, wfn, ham, grid, t)
            elseif params.hhg_method == :acceleration
                # Method 2a: Acceleration form (Fortran lines 905-909, commented)
                compute_dipole_acceleration!(obs, wfn, pot, grid, E_field)
            elseif params.hhg_method == :length
                # Method 2b: Length form
                compute_dipole_length!(obs, wfn, grid)
            else
                @warn "Unknown HHG method: $(params.hhg_method), skipping dipole calculation"
            end
        end

        # Record observables at specified intervals
        if step % params.obs_interval == 0
            record_observables!(obs, t, wfn, ham, pot, grid,
                              r_cutoff=params.ionization_cutoff)

            if verbose && step % (10 * params.obs_interval) == 0
                progress = 100.0 * step / n_steps
                @info @sprintf("Progress: %.1f%%  t=%.2f a.u.  norm=%.4f  P_ion=%.4e",
                              progress, t, obs.norms[end], obs.ionization_probs[end])
            end
        end
    end

    # Final observables
    record_observables!(obs, t, wfn, ham, pot, grid,
                       r_cutoff=params.ionization_cutoff)

    # ===== Step 7: Finalize =====
    wall_time = time() - start_time

    if verbose
        @info @sprintf("Simulation completed in %.2f seconds", wall_time)
        @info "Final state" norm=obs.norms[end] P_ion=obs.ionization_probs[end]

        if params.enable_momentum_analysis
            @info "Momentum projection accumulated" n_projections=volkov_projector.n_projections
        end

        println()
        print_observables_summary(obs)
    end

    # Extract momentum data if available
    momentum_grid = params.enable_momentum_analysis ? volkov_projector.momentum_grid : nothing
    psai_p = params.enable_momentum_analysis ? volkov_projector.psai_pt : nothing  # Already accumulated with phases

    return SimulationResults(params, grid, ham, wfn, obs, wall_time, momentum_grid, psai_p)
end

"""
    print_simulation_params(params::SimulationParams)

Print formatted summary of simulation parameters.
"""
function print_simulation_params(params::SimulationParams)
    println("Simulation Parameters:")
    println("=" ^ 60)

    println("\nGrid:")
    println("  Radial points: $(params.nrmax)")
    println("  Max radius: $(params.rmax) a.u.")
    println("  Angular momentum: lmax = $(params.lmax)")

    println("\nAtom:")
    println("  Type: $(params.atom)")

    println("\nTime propagation:")
    println("  Time step: $(params.dt) a.u. ($(@sprintf("%.2f", params.dt * 24.2)) as)")
    println("  Total time: $(params.t_total) a.u. ($(@sprintf("%.2f", params.t_total * 24.2)) as)")
    println("  Steps: $(Int(ceil(params.t_total / params.dt)))")

    println("\nLaser field:")
    println("  Mode: $(params.use_bicircular ? "Bicircular (ω + 2ω)" : "Simple (Gaussian)")")
    println("  Wavelength: $(params.laser_wavelength) nm")
    println("  Intensity: $(@sprintf("%.2e", params.laser_intensity)) W/cm²")
    println("  Duration: $(params.laser_duration) a.u. ($(@sprintf("%.2f", params.laser_duration * 24.2)) as)")

    if params.use_bicircular
        # Bicircular mode parameters
        println("  Amplitude ratios: Ex=$(params.exrate), Ey=$(params.eyrate), Ez=$(params.ezrate)")
        println("  Phase offset: $(params.phase_offset) rad ($(@sprintf("%.1f", params.phase_offset * 180/π))°)")
        println("  Envelope: Trapezoidal (4 optical cycle ramps)")

        # Determine polarization type
        if abs(params.eyrate) < 1e-10 && abs(params.ezrate) < 1e-10
            println("  Type: Linear (X-axis)")
        elseif abs(params.exrate) < 1e-10 && abs(params.ezrate) < 1e-10
            println("  Type: Linear (Y-axis, 2ω)")
        elseif abs(params.ezrate) > 1e-10
            println("  Type: 3D field")
        else
            ellipticity = abs(params.eyrate / params.exrate)
            if ellipticity ≈ 1.0 && abs(params.phase_offset - π/2) < 0.1
                println("  Type: Circular (bicircular)")
            else
                @printf("  Type: Elliptical/Bicircular (ε = %.3f, φ = %.1f°)\n",
                       ellipticity, params.phase_offset*180/π)
            end
        end
    else
        # Simple mode parameters
        println("  CEP: $(params.laser_cep) rad")
        println("  Polarization: $(params.laser_polarization)")
        println("  Envelope: Gaussian")
    end

    # Compute derived quantities
    ω = wavelength_nm_to_frequency_au(params.laser_wavelength)
    I_au = intensity_SI_to_au(params.laser_intensity)
    E0 = params.use_bicircular ? sqrt(I_au) : sqrt(2.0 * I_au)
    Up = E0^2 / (4 * ω^2)  # Ponderomotive energy

    println("\nDerived quantities:")
    println("  Photon energy: $(@sprintf("%.3f", ω * 27.2114)) eV")
    println("  Peak field: $(@sprintf("%.4f", E0)) a.u.")
    println("  Ponderomotive energy: $(@sprintf("%.3f", Up * 27.2114)) eV")

    println("\nNumerical:")
    println("  Energy cutoff: $(params.energy_cutoff) Ha")
    println("  Ionization cutoff: $(params.ionization_cutoff) a.u.")
    println("=" ^ 60)
end

"""
    generate_eigenstates(params::SimulationParams; verbose::Bool=true) -> HamiltonianData

Generate eigenstates and save to HDF5 file without running simulation.

This is a convenience wrapper for `run_simulation` with `eigenstate_mode = :generate_only`.
Use this function for pre-computing eigenstates on a cluster or workstation.

# Arguments
- `params::SimulationParams`: Simulation parameters (eigenstate_file determines output location)
- `verbose::Bool=true`: Print progress information

# Returns
- `HamiltonianData`: The computed Hamiltonian with eigenstates

# Example
```julia
# Pre-compute eigenstates for later use
params = create_default_params(
    atom = :hydrogen,
    nrmax = 400,
    lmax = 50,
    energy_cutoff = 10.0,
    eigenstate_file = "hydrogen_eigenstates.h5"
)

ham = generate_eigenstates(params)

# Later, run simulations using the pre-computed eigenstates:
params_sim = create_default_params(
    atom = :hydrogen,
    nrmax = 400,  # Must match!
    lmax = 50,    # Must match!
    eigenstate_mode = :load_from_file,
    eigenstate_file = "hydrogen_eigenstates.h5",
    t_total = 100.0
)
results = run_simulation(params_sim)
```
"""
function generate_eigenstates(params::SimulationParams; verbose::Bool=true)
    # Create a modified params with eigenstate_mode = :generate_only
    new_params = SimulationParams(;
        # Copy all fields from params
        nrmax = params.nrmax,
        rmax = params.rmax,
        L = params.L,
        lmax = params.lmax,
        atom = params.atom,
        custom_potential = params.custom_potential,
        dt = params.dt,
        t_total = params.t_total,
        obs_interval = params.obs_interval,
        laser_wavelength = params.laser_wavelength,
        laser_intensity = params.laser_intensity,
        laser_duration = params.laser_duration,
        laser_cep = params.laser_cep,
        laser_polarization = params.laser_polarization,
        use_bicircular = params.use_bicircular,
        exrate = params.exrate,
        eyrate = params.eyrate,
        ezrate = params.ezrate,
        phase_offset = params.phase_offset,
        energy_cutoff = params.energy_cutoff,
        ionization_cutoff = params.ionization_cutoff,
        # Override eigenstate_mode
        eigenstate_mode = :generate_only,
        eigenstate_file = params.eigenstate_file,
        # Analysis switches
        enable_momentum_analysis = params.enable_momentum_analysis,
        msplit = params.msplit,
        R_c = params.R_c,
        delta_split = params.delta_split,
        p_max = params.p_max,
        n_p = params.n_p,
        n_theta = params.n_theta,
        n_phi = params.n_phi,
        enable_hhg = params.enable_hhg,
        hhg_method = params.hhg_method,
        hhg_record_interval = params.hhg_record_interval
    )

    results = run_simulation(new_params, verbose=verbose)
    return results.hamiltonian
end

end  # module Simulation
