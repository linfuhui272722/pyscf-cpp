/**
 * PySCF C++ Demo - Using libcint for accurate GTO integrals
 */

#include <iostream>
#include <memory>
#include <string>
#include "pyscf_cpp.h"

using namespace pyscf;

int main(int argc, char** argv) {
    std::cout << "========================================\n";
    std::cout << "   PySCF C++ B3LYP Demo (libcint)\n";
    std::cout << "========================================\n\n";

    // Get basis set from command line or use default
    std::string basis = (argc > 1) ? argv[1] : "sto-3g";
    
    // Test with hydrogen molecule - small system for testing
    std::cout << "Creating hydrogen molecule (H2)...\n";
    auto mol = std::make_shared<Molecule>();

    // H2 at equilibrium bond length (0.735 Angstrom = 1.389 Bohr)
    // Convert Angstrom to Bohr (1 Angstrom = 1.889726 Bohr)
    double bohr = 1.889726;
    mol->add_atom(1, 0.0, 0.0, 0.0);       // H1 at origin
    mol->add_atom(1, 0.735 * bohr, 0.0, 0.0);  // H2 at 0.735 Angstrom in Bohr

    std::cout << "  Atoms added: " << mol->num_atoms() << "\n";
    std::cout << "  Electrons: " << mol->num_electrons() << "\n";
    std::cout << "  Nuclear repulsion: " << mol->nuclear_repulsion_energy() << " Eh\n\n";

    // Set basis set
    std::cout << "Loading basis set (" << basis << ")...\n";
    mol->set_basis(basis);
    std::cout << "  Basis functions: " << mol->num_basis_functions() << "\n\n";

    if (mol->num_basis_functions() == 0) {
        std::cerr << "ERROR: No basis functions loaded!\n";
        return 1;
    }

    // Create RKS object and set options
    auto scf = std::make_shared<dft::RKS>(mol);

    dft::SCFOptions options;
    options.conv_tol = 1e-6;
    options.max_cycle = 50;
    options.verbose = true;
    scf->set_options(options);

    // Run SCF calculation
    std::cout << "Running B3LYP SCF calculation...\n";
    auto result = scf->compute();

    std::cout << "\nSCF Results:\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Total energy: " << result.energy << " Eh\n";
    std::cout << "  Electronic energy: " << scf->get_electronic_energy() << " Eh\n";
    std::cout << "  Nuclear repulsion: " << scf->get_nuclear_repulsion() << " Eh\n\n";

    // Compare with reference values
    if (basis == "sto-3g") {
        std::cout << "Reference H2/STO-3G energy: -1.066 Eh\n";
    }

    std::cout << "\n========================================\n";
    std::cout << "Demo complete!\n";
    std::cout << "========================================\n";

    return 0;
}
