"""
Master include file for all TDSE solver modules.

Include this file to load all modules in the correct dependency order.
"""

# Foundational utilities
include("utils/PhysicalUnits.jl")
using .PhysicalUnits

include("utils/Validation.jl")
using .Validation

# Grid modules
include("grid/GPSGrid.jl")
using .GPSGrid

include("grid/AngularGrid.jl")
using .AngularGrid

# Laser field module (before Propagator)
include("field/LaserField.jl")
using .LaserField

# Hamiltonian components
include("hamiltonian/Potential.jl")
using .Potential

# Wavefunction representation
include("wavefunction/Wavefunction.jl")
using .Wavefunction

# Hamiltonian (depends on GPSGrid, Potential)
# Import at parent scope so Hamiltonian can reference them
import .GPSGrid: GPSGridData
import .Potential: PotentialFunction

include("hamiltonian/Hamiltonian.jl")
using .Hamiltonian

# Coordinate transformation for field interaction (uses Wavefunction, AngularGrid, GPSGrid)
include("propagator/CoordinateTransform.jl")
using .CoordinateTransform

# Propagator (depends on Hamiltonian, AngularGrid, Wavefunction, LaserField, CoordinateTransform)
include("propagator/Propagator.jl")
using .Propagator

# Absorbing boundary for HHG calculations (depends on GPSGrid, Wavefunction)
include("propagator/AbsorbingBoundary.jl")
using .AbsorbingBoundary

# Observables (depends on Wavefunction, Hamiltonian, Potential, GPSGrid)
include("observables/Observables.jl")
using .Observables

# HHG module (depends on Observables, Wavefunction, Hamiltonian, GPSGrid, Potential)
include("observables/HHG.jl")
using .HHG

# Region splitting for continuum analysis (depends on Wavefunction, GPSGrid)
include("region_split/RegionSplit.jl")
using .RegionSplit

# Continuum analysis modules
include("continuum/VolkovProjection.jl")
using .VolkovProjection

include("continuum/MomentumDistribution.jl")
using .MomentumDistribution

# Simulation orchestrator (depends on all modules)
include("simulation/Simulation.jl")
using .Simulation

# I/O modules
include("io/ConfigParser.jl")
using .ConfigParser

include("io/ResultsIO.jl")
using .ResultsIO

# Parameter scanning module (depends on Simulation and ConfigParser)
include("scan/ParameterScan.jl")
using .ParameterScan

println("✓ All TDSE solver modules loaded successfully")
