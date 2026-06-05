# TDSESolver.jl

**Time-Dependent Schrödinger Equation Solver in Julia**

A modern Julia implementation of TDSE solver for computing atomic ionization, photoelectron momentum distributions, and high-harmonic generation (HHG) spectra in strong laser fields.

## Status

> **⚠️ WARNING: This project is unfinished and at a very early stage.**
>
> The codebase implements the core physics infrastructure but has **not been validated** against the Fortran reference benchmarks it was ported from. Key validation gates (correlation > 0.999 with Fortran outputs, ground state energy accuracy, norm conservation under strong fields) remain unconfirmed. **This code is not yet ready for research use.**

✅ **Infrastructure Complete** (Phases 1-6)
📝 **Documentation Complete** (Phases 7-8)
⚠️ **Physical Validation Pending** (See Known Limitations)

## Features

- 🎯 **Complete Pipeline**: Configuration → Propagation → Analysis → Results
- ⚡ **Multiple Analysis Modes**:
  - Ionization probability
  - Photoelectron momentum distributions (ATI spectra)
  - High-harmonic generation (HHG) spectra
  - Combined analysis
- 🔬 **GPS Discretization**: Gauss Pseudospectral method with algebraic mapping
- 📊 **Split-Operator Propagation**: S-matrix (field-free) + coordinate-space (field interaction)
- 🔧 **Flexible Configuration**: TOML files + programmatic API
- 🎨 **Custom Potentials**: User-defined V(r) via string expressions or Julia modules
- 🔄 **Parameter Scanning**: Systematic exploration with parallel execution support

## Quick Start

### Installation
```bash
# Clone repository
git clone <repository-url>
cd TDSE

# No external packages required (stdlib only: LinearAlgebra, FFTW, Distributed)
```

### Basic Usage
```julia
include("src/TDSESolver.jl")
using .Simulation

# Define parameters
params = create_default_params(
    atom = :hydrogen,
    laser_wavelength = 800.0,   # nm
    laser_intensity = 1.0e14,   # W/cm²
    laser_duration = 50.0,      # a.u.
    t_total = 200.0
)

# Run simulation
results = run_simulation(params, verbose=true)

# Access results
final_ionization = results.observables.ionization_probs[end]
println("Ionization probability: $final_ionization")
```

### Using TOML Configuration
```julia
include("src/TDSESolver.jl")
using .ConfigParser, .ParameterScan

# Single simulation
config = load_config("examples/helium_800nm.toml")
results = run_parameter_scan(config)

# Parameter scan
scan_config = load_config("examples/intensity_scan.toml")
scan_results = run_parameter_scan(scan_config, parallel=false)
print_scan_summary(scan_results)
```

## Documentation

### Quick References
- **[Quick Reference](docs/QUICK_REFERENCE.md)** ⭐ - Fast lookup for common tasks
- **[Project Summary](docs/PROJECT_COMPLETION_SUMMARY.md)** - Complete feature overview
- **[Known Limitations](docs/UNSOLVED_PROBLEMS.md)** ⚠️ - GPS eigenstate issues

### Technical Details
- **[Algorithm References](docs/algorithm_references.md)** - Mathematical foundations
- **[GPS Investigation](docs/GPS_INVESTIGATION.md)** - Eigenstate accuracy analysis
- **[Project Status Report](docs/PROJECT_STATUS_REPORT.md)** - Detailed implementation status

## Examples

See [examples/](examples/) directory for:
- `helium_800nm.toml` - Basic helium ionization
- `helium_hhg.toml` - High-harmonic generation
- `intensity_scan.toml` - Parameter scans
- `soft_core.toml` - Custom potential
- `custom_potential.jl` - Advanced Julia module approach
- `elliptical_polarization.jl` - Elliptical and bicircular field examples
- `bicircular_test.jl` - Bicircular field integration tests

## Citation

If you use TDSESolver.jl in your research, please cite:

```bibtex
@software{tdsesolver_jl,
  title = {TDSESolver.jl: Research-Grade TDSE Solver with User-Friendly Julia API},
  author = {TDSE Research Team},
  year = {2025},
  url = {https://github.com/your-org/TDSESolver.jl}
}
```

## License

MIT License

## Acknowledgments

This implementation follows algorithms from:
- Tong & Chu (1997) - Split-operator method
- PRA 74.031405(R)(2006) - Region splitting method
- Original Fortran blueprints in `explore/` directory
