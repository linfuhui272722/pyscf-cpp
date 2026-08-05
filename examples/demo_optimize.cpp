/**
 * Demo: Geometry Optimization
 * 
 * Demonstrates the BFGS geometry optimizer on H2 molecule.
 * Compares with PySCF results.
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
    std::cout << "  Geometry Optimization Demo (BFGS)\n";
    std::cout << "========================================\n\n";
    
    // Create H2 molecule with initial bond length
    auto mol = std::make_shared<Molecule>();
    mol->add_atom(1, 0.0, 0.0, 0.0);      // H at origin
    mol->add_atom(1, 1.4, 0.0, 0.0);      // H at 1.4 Bohr (experiment ~1.4 Bohr = 0.74 Ang)
    mol->set_basis("sto-3g");
    
    std::cout << "Molecule: H2\n";
    std::cout << "Basis set: STO-3G\n";
    std::cout << "Initial H-H distance: " << 1.4 << " Bohr\n\n";
    
    // Create optimizer FIRST (this creates a copy of the molecule)
    geom::BFGSOptimizer optimizer;
    geom::OptConvergence conv;
    conv.gradient_max = 4.5e-4;
    conv.gradient_rms = 3.0e-4;
    conv.step_max = 1.8e-3;
    conv.step_rms = 1.2e-3;
    conv.energy_tol = 1.0e-6;
    conv.max_iterations = 50;
    conv.verbose = true;
    optimizer.set_convergence(conv);
    
    // Create SCF object with the molecule that optimizer will use
    // Note: optimizer.optimize() creates its own copy of mol, but we pass it here
    // to initialize the SCF. The optimizer will update the same mol object.
    auto scf = std::make_shared<dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    
    // Initial energy
    scf->compute();
    std::cout << "Initial energy: " << std::scientific << std::setprecision(8)
              << scf->get_total_energy() << " Eh\n\n";
    
    // Run optimization
    std::cout << "Starting geometry optimization...\n\n";
    auto result = optimizer.optimize(scf, mol);
    
    // Results
    std::cout << "\n========================================\n";
    std::cout << "           Optimization Results\n";
    std::cout << "========================================\n\n";
    std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final energy: " << std::scientific << std::setprecision(8)
              << result.final_energy << " Eh\n";
    
    // Compute final distances
    auto final_mol = std::make_shared<Molecule>(*mol);
    final_mol->set_coordinates(result.coordinates);
    auto atom1 = final_mol->get_atom(0);
    auto atom2 = final_mol->get_atom(1);
    double dx = atom2.x - atom1.x;
    double dy = atom2.y - atom1.y;
    double dz = atom2.z - atom1.z;
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "Final H-H distance: " << std::fixed << std::setprecision(4)
              << distance << " Bohr (" << distance * 0.529177 << " Angstrom)\n";
    
    std::cout << "\nEnergy history:\n";
    for (size_t i = 0; i < result.energy_history.size(); ++i) {
        std::cout << "  Cycle " << i << ": E = " << std::scientific 
                  << std::setprecision(8) << result.energy_history[i];
        std::cout << "  |g| = " << result.gradient_norms[i] << "\n";
    }
    
    // Compare with PySCF reference
    std::cout << "\n========================================\n";
    std::cout << "           PySCF Reference\n";
    std::cout << "========================================\n";
    std::cout << "Expected H-H distance: ~1.85 Bohr (0.98 Angstrom)\n";
    std::cout << "Expected energy: ~-1.028 Eh (STO-3G)\n";
    
    return 0;
}
