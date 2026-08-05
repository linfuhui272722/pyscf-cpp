# PySCF-C++ Tutorial: Writing C++ Programs for Molecular Calculations

> **⚠️ DISCLAIMER**: This tutorial was written by AI and does not guarantee safety, stability, or scientific accuracy. This code is for learning and research purposes only. Do not use it for production environments or scientific publications. For serious quantum chemistry calculations, please use validated software such as Gaussian, ORCA, PySCF, etc.

## Table of Contents

1. [Introduction](#1-introduction)
2. [Development Environment Setup](#2-development-environment-setup)
3. [Project Structure](#3-project-structure)
4. [Core Concepts](#4-core-concepts)
5. [Basic Example: Water Molecule Energy Calculation](#5-basic-example-water-molecule-energy-calculation)
6. [Advanced Example: Geometry Optimization](#6-advanced-example-geometry-optimization)
7. [Advanced Example: Mayer Bond Order Analysis](#7-advanced-example-mayer-bond-order-analysis)
8. [Complete Example Programs](#8-complete-example-programs)
9. [Basis Set Handling](#9-basis-set-handling)
10. [Common Issues and Debugging](#10-common-issues-and-debugging)

---

## 1. Introduction

This project is a C++ implementation of density functional theory (DFT) calculations, inspired by PySCF's design. Main features include:

- **B3LYP Functional**: The most commonly used hybrid functional
- **Multiple Basis Sets**: STO-3G, 6-31G, 6-31G*, 6-31G**, 3-21G, etc.
- **Geometry Optimization**: BFGS-based molecular structure optimization
- **Bond Order Analysis**: Mayer bond orders and Wiberg indices

### 1.1 System Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    User Program (main.cpp)                │
├─────────────────────────────────────────────────────────┤
│  Molecule │ RKS (SCF) │ Optimizer │ BondAnalyzer        │
├─────────────────────────────────────────────────────────┤
│                    Core Calculation Layer                 │
├─────────────────────────────────────────────────────────┤
│  GTO Integrals │ Lebedev Grid │ B3LYP XC │ BLAS/LAPACK │
├─────────────────────────────────────────────────────────┤
│              libcint (Integral Library) / OpenBLAS      │
└─────────────────────────────────────────────────────────┘
```

---

## 2. Development Environment Setup

### 2.1 Dependencies

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    cmake \
    g++ \
    libopenblas-dev \
    liblapack-dev \
    libcint-dev \
    make
```

### 2.2 Build from Source (libcint)

```bash
git clone https://github.com/sunqm/libcint.git
cd libcint
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make -j4
sudo make install
```

### 2.3 Build Project

```bash
cd pyscf-cpp
mkdir build && cd build
cmake ..
make -j4
```

---

## 3. Project Structure

```
pyscf-cpp/
├── include/                    # Header files
│   ├── molecule.h              # Molecular structure
│   ├── scf.h                  # SCF interface
│   ├── optimizer.h             # Geometry optimizer
│   ├── mayer_bond.h            # Bond order analysis
│   ├── b3lyp_xc.h             # B3LYP functional
│   └── lebedev_grid.h          # Numerical integration grid
├── src/                       # Source files
│   ├── molecule.cpp
│   ├── scf.cpp
│   ├── gto/                    # GTO integrals
│   ├── dft/                    # DFT calculations
│   ├── geom/                   # Geometry optimization
│   └── bond/                   # Bond analysis
├── basis/                     # Basis set data files
├── examples/                  # Example programs
└── CMakeLists.txt
```

---

## 4. Core Concepts

### 4.1 Atomic Units

Quantum chemistry typically uses atomic units (a.u.):
- Length: 1 Bohr = 0.529177 Å
- Energy: 1 Eh (Hartree) = 27.2114 eV

### 4.2 Gaussian-Type Orbitals (GTO)

Basis functions are expressed as linear combinations of Gaussian functions:

```
χ(r) = Σ c_i × g_i(r)
g_i(r) = (x-ax)^l (y-ay)^m (z-az)^n × exp(-α_i × |r-R|^2)
```

### 4.3 Density Matrix

For closed-shell systems:

```
P_μν = 2 × Σ_i(occ) C_μi × C_νi
```

### 4.4 B3LYP Functional

B3LYP = 20% HF + 72% LDA/GGA + 8% GGA exchange

```
E_XC^B3LYP = (1-a) × E_X^LDA + a × E_X^HF 
            + c × ΔE_X^B88 + (1-c) × E_C^VWN 
            + d × ΔE_C^LYP
where a=0.20, c=0.72, d=0.81
```

---

## 5. Basic Example: Water Molecule Energy Calculation

### 5.1 Complete Code

```cpp
/**
 * Basic Example: Water Molecule Energy Calculation
 * 
 * This program demonstrates how to perform basic DFT energy calculations
 * using PySCF-C++
 */

#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

#include "molecule.h"
#include "scf.h"

using namespace pyscf;

int main() {
    std::cout << "========================================\n";
    std::cout << "  Water Molecule B3LYP/6-31G Energy\n";
    std::cout << "========================================\n\n";
    
    // ========== Step 1: Create molecule object ==========
    // 1. Define molecular structure
    // Coordinates are in Bohr (1 Bohr = 0.529 Å)
    
    auto mol = std::make_shared<Molecule>();
    
    // Add atoms: atomic_number, x, y, z (in Bohr)
    mol->add_atom(8,  0.0,  0.0,  0.0);       // O at origin
    mol->add_atom(1,  1.4,  1.2,  0.0);      // H1
    mol->add_atom(1, -1.4,  1.2,  0.0);      // H2
    
    // 2. Set basis set
    mol->set_basis("6-31g");
    
    // Print molecule info
    std::cout << "Molecule Info:\n";
    std::cout << "  Atoms: " << mol->num_atoms() << "\n";
    std::cout << "  Electrons: " << mol->num_electrons() << "\n";
    std::cout << "  Basis functions: " << mol->num_basis_functions() << "\n";
    std::cout << "  Nuclear repulsion: " << mol->nuclear_repulsion_energy() << " Eh\n\n";
    
    // ========== Step 2: Create SCF object ==========
    auto scf = std::make_shared<dft::RKS>(mol);
    
    // Set exchange-correlation functional
    scf->set_xc_functional("b3lyp");
    
    // Set SCF options (optional)
    dft::SCFOptions opts;
    opts.conv_tol = 1e-8;
    opts.max_cycle = 100;
    opts.verbose = true;
    scf->set_options(opts);
    
    // ========== Step 3: Run calculation ==========
    std::cout << "Starting SCF iterations...\n\n";
    
    scf->compute();
    
    // ========== Step 4: Get results ==========
    std::cout << "\n========================================\n";
    std::cout << "            Results\n";
    std::cout << "========================================\n\n";
    
    std::cout << "Total energy: " << std::scientific << std::setprecision(8)
              << scf->get_total_energy() << " Eh\n";
    std::cout << "Electronic energy: " << scf->get_electronic_energy() << " Eh\n";
    std::cout << "Nuclear repulsion: " << scf->get_nuclear_repulsion() << " Eh\n\n";
    
    return 0;
}
```

### 5.2 Code Explanation

#### 5.2.1 Include Headers

```cpp
#include "molecule.h"   // Molecule definition
#include "scf.h"        // SCF calculation
```

#### 5.2.2 Create Molecule

```cpp
auto mol = std::make_shared<Molecule>();
mol->add_atom(8, 0.0, 0.0, 0.0);  // Atomic number 8 = Oxygen
```

**Atomic Number Reference:**
| Element | Z | Element | Z |
|---------|---|---------|---|
| H | 1 | S | 16 |
| C | 6 | Cl | 17 |
| N | 7 | Br | 35 |
| O | 8 | I | 53 |

#### 5.2.3 Set Basis

```cpp
mol->set_basis("sto-3g");    // Minimal basis
mol->set_basis("6-31g");     // Pople basis
mol->set_basis("6-31g*");   // With polarization
mol->set_basis("6-31g**");  // Double polarization
```

#### 5.2.4 Create SCF and Compute

```cpp
auto scf = std::make_shared<dft::RKS>(mol);
scf->set_xc_functional("b3lyp");
scf->compute();
```

---

## 6. Advanced Example: Geometry Optimization

### 6.1 Complete Code

```cpp
/**
 * Geometry Optimization Example
 * 
 * Using BFGS algorithm to optimize molecular geometry
 */

#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

#include "molecule.h"
#include "scf.h"
#include "optimizer.h"

using namespace pyscf;

int main() {
    std::cout << "========================================\n";
    std::cout << "  H2 Geometry Optimization (BFGS)\n";
    std::cout << "========================================\n\n";
    
    // ========== Create molecule ==========
    auto mol = std::make_shared<Molecule>();
    
    // Initial guess: H-H distance 1.4 Bohr (experimental ~1.4 Bohr)
    mol->add_atom(1, 0.0, 0.0, 0.0);      // H1
    mol->add_atom(1, 1.4, 0.0, 0.0);     // H2
    mol->set_basis("sto-3g");
    
    // ========== Create SCF ==========
    auto scf = std::make_shared<dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    
    // ========== Create optimizer ==========
    geom::BFGSOptimizer optimizer;
    
    // Set convergence criteria
    geom::OptConvergence conv;
    conv.gradient_max = 4.5e-4;   // Max gradient (Eh/Bohr)
    conv.gradient_rms = 3.0e-4;   // RMS gradient
    conv.step_max = 1.8e-3;       // Max step (Bohr)
    conv.step_rms = 1.2e-3;       // RMS step
    conv.energy_tol = 1.0e-6;      // Energy tolerance
    conv.max_iterations = 50;
    conv.verbose = true;
    optimizer.set_convergence(conv);
    
    // ========== Run optimization ==========
    std::cout << "Starting geometry optimization...\n\n";
    
    auto result = optimizer.optimize(scf, mol);
    
    // ========== Output results ==========
    std::cout << "\n========================================\n";
    std::cout << "           Optimization Results\n";
    std::cout << "========================================\n\n";
    
    std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final energy: " << std::scientific << std::setprecision(8)
              << result.final_energy << " Eh\n";
    
    return 0;
}
```

### 6.2 BFGS Algorithm Principles

BFGS is a quasi-Newton method that accelerates convergence by approximating the Hessian:

```
Search direction: p_k = -H_k × g_k
Step length: α_k determined by line search
Update: H_{k+1} = H_k + (y_k × y_k^T) / (y_k^T × s_k) 
                - (H_k × s_k × s_k^T × H_k) / (s_k^T × H_k × s_k)
where: s_k = x_{k+1} - x_k
       y_k = g_{k+1} - g_k
```

### 6.3 Convergence Criteria

| Parameter | Default | Description |
|-----------|---------|-------------|
| gradient_max | 4.5e-4 | Max gradient component (Eh/Bohr) |
| gradient_rms | 3.0e-4 | RMS gradient |
| step_max | 1.8e-3 | Max displacement (Bohr) |
| step_rms | 1.2e-3 | RMS displacement |
| energy_tol | 1.0e-6 | Energy change |

---

## 7. Advanced Example: Mayer Bond Order Analysis

### 7.1 Complete Code

```cpp
/**
 * Mayer Bond Order Analysis Example
 * 
 * Calculate bond orders and types in molecules
 */

#include <iostream>
#include <iomanip>
#include <memory>

#include "molecule.h"
#include "scf.h"
#include "mayer_bond.h"

using namespace pyscf;

int main() {
    std::cout << "========================================\n";
    std::cout << "  Mayer Bond Order Analysis Example\n";
    std::cout << "========================================\n\n";
    
    // ========== Water molecule ==========
    std::cout << "--- Water Molecule (H2O) ---\n\n";
    
    auto mol = std::make_shared<Molecule>();
    mol->add_atom(8,  0.0,  0.0,  0.0);   // O
    mol->add_atom(1,  1.4,  1.2,  0.0);   // H1
    mol->add_atom(1, -1.4,  1.2,  0.0);   // H2
    mol->set_basis("sto-3g");
    
    // SCF calculation
    auto scf = std::make_shared<dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    scf->compute();
    
    std::cout << "Total energy: " << scf->get_total_energy() << " Eh\n\n";
    
    // Bond order analysis
    auto result = bond::analyze_bonds(scf, mol);
    bond::print_bond_analysis(result, *mol);
    
    return 0;
}
```

### 7.2 Mayer Bond Order Formula

Mayer bond order is defined as:

```
W_AB = Σ_μ∈A Σ_ν∈B P_μν × S_μν

where:
- P_μν: Density matrix element
- S_μν: Overlap matrix element
```

### 7.3 Bond Order Reference Values

| Molecule | Bond | Mayer BO | Type |
|----------|------|----------|------|
| H2 | H-H | ~1.0 | Single |
| H2O | O-H | ~0.8 | Single |
| N2 | N≡N | ~3.0 | Triple |
| CO | C≡O | ~2.5 | Triple |
| C2H4 | C=C | ~2.0 | Double |

---

## 8. Complete Example Programs

### 8.1 Energy Calculation Template

```cpp
#include <iostream>
#include <memory>
#include "molecule.h"
#include "scf.h"

int main() {
    // 1. Create molecule
    auto mol = std::make_shared<pyscf::Molecule>();
    // mol->add_atom(atomic_number, x, y, z);
    mol->set_basis("basis_name");
    
    // 2. Create SCF calculator
    auto scf = std::make_shared<pyscf::dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    
    // 3. Run calculation
    scf->compute();
    
    // 4. Get results
    std::cout << "Total energy: " << scf->get_total_energy() << " Eh\n";
    
    return 0;
}
```

### 8.2 Common Molecular Geometries

```cpp
// Water molecule (experimental)
O  0.0000  0.0000  0.0000
H  0.7570  0.5860  0.0000
H -0.7570  0.5860  0.0000
(Units: Angstrom)

// Methane (tetrahedral)
C  0.0000  0.0000  0.0000
H  0.6290  0.6290  0.6290
H -0.6290 -0.6290  0.6290
H -0.6290  0.6290 -0.6290
H  0.6290 -0.6290 -0.6290

// Ammonia (trigonal pyramidal)
N  0.0000  0.0000  0.0000
H  0.0000  0.9390  0.8110
H  0.8130 -0.4700  0.8110
H -0.8130 -0.4700  0.8110
```

### 8.3 Compilation

```bash
# Compile
g++ -std=c++17 -O3 -I./include \
    main.cpp \
    src/molecule.cpp \
    src/dft/scf.cpp \
    src/dft/b3lyp_xc.cpp \
    src/dft/lebedev_grid.cpp \
    src/gto/gto_integrals.cpp \
    -o my_program \
    -L/usr/local/lib -lcint \
    -lopenblas -llapack -lm -lpthread

# Run
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH ./my_program
```

---

## 9. Basis Set Handling

### 9.1 Supported Basis Sets

| Basis | Type | Description |
|-------|------|-------------|
| STO-3G | Minimal | Earliest and simplest basis |
| 3-21G | Split valence | Two Gaussians for inner, one for valence |
| 6-31G | Split valence | More precise inner and valence |
| 6-31G* | Polarization | Adds d polarization |
| 6-31G** | Double pol. | Adds p polarization to H |
| 6-311G | Triple split | Triple split valence |

### 9.2 Basis Set File Format

```text
# STO-3G Basis Example
H S
  3
  18.7311370  0.03349460
   2.8253937  0.23472695
   0.6401217  0.81375733
  1
   0.1612778  1.00000000
```

---

## 10. Common Issues and Debugging

### 10.1 SCF Not Converging

**Problem**: Calculation produces NaN or SCF doesn't converge

**Possible causes**:
1. Initial geometry is unreasonable (atoms too close)
2. Basis set unsuitable for the molecule
3. SCF parameters need adjustment

**Solutions**:
```cpp
// Lower convergence requirements
opts.conv_tol = 1e-6;
opts.max_cycle = 200;
scf->set_options(opts);
```

### 10.2 Out of Memory

**Problem**: Too large basis set causes memory overflow

**Solutions**:
- Use smaller basis (sto-3g instead of cc-pVTZ)
- Reduce integral precision
- Use density fitting approximations

### 10.3 Compilation Errors

**Problem**: Cannot find libcint

**Solutions**:
```bash
# Check installation
ldconfig -p | grep cint
# Or manually specify
cmake -DCINT_LIBRARY=/path/to/libcint.so ..
```

### 10.4 Abnormal Energy Values

**Problem**: Calculated energy differs significantly from expectation

**Check items**:
- [ ] Coordinate units correct (Bohr)
- [ ] Basis set loaded correctly
- [ ] Nuclear charges correct
- [ ] Electron count correct

---

## Appendix: Complete API Reference

### A.1 Molecule Class

```cpp
class Molecule {
public:
    void add_atom(int atomic_number, double x, double y, double z);
    void set_basis(const std::string& basis_name);
    int num_atoms() const;
    int num_electrons() const;
    int num_basis_functions() const;
    double nuclear_repulsion_energy() const;
    std::vector<double> get_coordinates() const;
    void set_coordinates(const std::vector<double>& coords);
};
```

### A.2 RKS Class

```cpp
class RKS {
public:
    RKS(std::shared_ptr<Molecule> mol);
    void set_xc_functional(const std::string& xc);
    SCFResult compute();
    double get_total_energy() const;
    double get_electronic_energy() const;
    double get_nuclear_repulsion() const;
    const std::vector<double>& get_density_matrix() const;
    const std::vector<double>& get_overlap_matrix() const;
    const std::vector<double>& get_mo_coefficients() const;
};
```

### A.3 BFGSOptimizer Class

```cpp
class BFGSOptimizer {
public:
    OptimizationResult optimize(
        std::shared_ptr<dft::RKS> scf,
        std::shared_ptr<Molecule> mol
    );
    void set_convergence(const OptConvergence& conv);
    void set_verbose(bool v);
};
```

### A.4 Bond Order Functions

```cpp
BondAnalysisResult analyze_bonds(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol
);
void print_bond_analysis(
    const BondAnalysisResult& result,
    const Molecule& mol
);
```

---

## References

1. Becke, A. D. (1993). "A new mixing of Hartree-Fock and local density-functional theories". The Journal of Chemical Physics, 98(2), 1372.
2. Lee, C., Yang, W., & Parr, R. G. (1988). "Development of the Colle-Salvetti correlation-energy formula into a functional of the electron density". Physical Review B, 37(2), 785.
3. Mayer, I. (1986). "Bond order and valence indices: a personal account". Journal of Computational Chemistry, 28(1), 204-221.
4. Nocedal, J., & Wright, S. J. (2006). Numerical Optimization. Springer.

---

**⚠️ DISCLAIMER**: This tutorial was written by AI (OpenHands) and is for reference only. The code does not guarantee:
- Scientific calculation accuracy
- Numerical stability
- Production environment safety

For serious scientific research, please use validated commercial or open-source software.

*Last updated: August 2026*
