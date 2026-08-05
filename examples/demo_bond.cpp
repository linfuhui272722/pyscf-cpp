/**
 * Demo: Mayer Bond Order Analysis
 * 
 * Demonstrates bond order analysis on various molecules:
 * - H2 (single bond)
 * - H2O (two O-H bonds)
 * - CO (triple bond)
 * 
 * Reference:
 * - Mayer bond orders: J. Chem. Phys. 85, 5583 (1986)
 * - Wiberg indices: Tetrahedron 24, 1083 (1968)
 */

#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>

#include "molecule.h"
#include "scf.h"
#include "mayer_bond.h"

using namespace pyscf;

// Run bond analysis on a molecule
void analyze_molecule(
    const std::string& name,
    const std::string& formula,
    const std::string& basis,
    const std::vector<std::tuple<int, double, double, double>>& atoms
) {
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "                    " << name << " (" << formula << ")\n";
    std::cout << "================================================================================\n";
    
    // Create molecule
    auto mol = std::make_shared<Molecule>();
    for (const auto& atom : atoms) {
        mol->add_atom(std::get<0>(atom),
                      std::get<1>(atom),
                      std::get<2>(atom),
                      std::get<3>(atom));
    }
    mol->set_basis(basis);
    
    std::cout << "Basis: " << basis << "\n";
    std::cout << "Atoms: " << mol->num_atoms() << "\n";
    std::cout << "Basis functions: " << mol->num_basis_functions() << "\n\n";
    
    // Create SCF object
    auto scf = std::make_shared<dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    
    // Run SCF
    scf->compute();
    double energy = scf->get_total_energy();
    std::cout << "SCF energy: " << std::scientific << std::setprecision(8)
              << energy << " Eh\n\n";
    
    // Perform bond analysis
    auto result = bond::analyze_bonds(scf, mol);
    
    // Print detailed results
    bond::print_bond_analysis(result, *mol);
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Mayer Bond Order Analysis Demo\n";
    std::cout << "========================================\n";
    std::cout << "\nThis program computes Mayer bond orders and Wiberg indices\n";
    std::cout << "for various molecules to analyze chemical bonding.\n";
    
    // ========================================================================
    // H2 - Single bond
    // ========================================================================
    analyze_molecule(
        "Hydrogen Molecule",
        "H2",
        "sto-3g",
        {
            {1, 0.0, 0.0, 0.0},           // H at origin
            {1, 1.4, 0.0, 0.0}            // H at 1.4 Bohr
        }
    );
    
    // ========================================================================
    // H2O - Water molecule (bent)
    // ========================================================================
    // Geometry: O at origin, H's at ~104.5 degree angle
    double angle = 104.5 * M_PI / 180.0;
    double OH = 1.81;  // Bohr
    double hx = OH * std::sin(angle / 2);
    double hy = OH * std::cos(angle / 2);
    
    analyze_molecule(
        "Water Molecule",
        "H2O",
        "sto-3g",
        {
            {8,  0.0,  0.0,  0.0},        // O at origin
            {1, -hx,   hy,  0.0},         // H1
            {1,  hx,   hy,  0.0}          // H2
        }
    );
    
    // ========================================================================
    // CH4 - Methane (tetrahedral)
    // ========================================================================
    // Tetrahedral geometry, C-H bond length ~2.05 Bohr
    double CH = 2.05;
    double t = std::sqrt(8.0/9.0);
    double u = 1.0/3.0;
    double v = std::sqrt(1.0 - t*t - u*u);
    
    analyze_molecule(
        "Methane",
        "CH4",
        "sto-3g",
        {
            {6,  0.0,  0.0,  0.0},         // C at origin
            {1,  CH*t,  CH*u,  CH*v},      // H1
            {1, -CH*t,  CH*u,  CH*v},      // H2
            {1,  0.0, -CH*t,  CH*v},       // H3
            {1,  0.0,  0.0, -CH*t - CH*u}  // H4
        }
    );
    
    // ========================================================================
    // N2 - Nitrogen molecule (triple bond)
    // ========================================================================
    analyze_molecule(
        "Nitrogen Molecule",
        "N2",
        "sto-3g",
        {
            {7, 0.0, 0.0, 0.0},            // N at origin
            {7, 2.0, 0.0, 0.0}             // N at 2.0 Bohr
        }
    );
    
    // ========================================================================
    // CO - Carbon monoxide (triple bond)
    // ========================================================================
    analyze_molecule(
        "Carbon Monoxide",
        "CO",
        "sto-3g",
        {
            {6, 0.0, 0.0, 0.0},            // C at origin
            {8, 2.0, 0.0, 0.0}            // O at 2.0 Bohr
        }
    );
    
    // ========================================================================
    // Summary
    // ========================================================================
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "           Summary of Bond Orders\n";
    std::cout << "========================================\n\n";
    
    std::cout << "Molecule | Bond | Mayer | Wiberg | Type\n";
    std::cout << "---------|------|-------|--------|------\n";
    std::cout << "H2       | H-H  | ~1.0  | ~1.0   | Single\n";
    std::cout << "H2O      | O-H  | ~0.8  | ~0.8   | Single\n";
    std::cout << "CH4      | C-H  | ~0.9  | ~0.9   | Single\n";
    std::cout << "N2       | N-N  | ~3.0  | ~3.0   | Triple\n";
    std::cout << "CO       | C-O  | ~2.5  | ~2.5   | Triple\n";
    
    std::cout << "\n";
    std::cout << "Note: Mayer bond orders > 1.5 typically indicate double bonds,\n";
    std::cout << "      and > 2.5 indicate triple bonds.\n";
    
    return 0;
}
