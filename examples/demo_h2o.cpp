#include <iostream>
#include <memory>
#include <string>
#include "pyscf_cpp.h"

using namespace pyscf;

int main(int argc, char** argv) {
    std::string basis = (argc > 1) ? argv[1] : "sto-3g";
    
    std::cout << "========================================\n";
    std::cout << "   PySCF C++ B3LYP Demo (H2O)\n";
    std::cout << "========================================\n\n";
    
    auto mol = std::make_shared<Molecule>();
    
    double bohr = 1.889726;
    double r_oh = 0.958 * bohr;
    double theta = 104.5 * M_PI / 180.0;
    double h_x = r_oh * sin(theta/2);
    double h_z = r_oh * cos(theta/2);
    
    mol->add_atom(8, 0.0, 0.0, 0.0);
    mol->add_atom(1, h_x, 0.0, h_z);
    mol->add_atom(1, -h_x, 0.0, h_z);
    
    std::cout << "Creating water molecule...\n";
    std::cout << "  Atoms: " << mol->num_atoms() << ", Electrons: " << mol->num_electrons() << "\n";
    
    mol->set_basis(basis);
    std::cout << "  Basis functions: " << mol->num_basis_functions() << "\n\n";
    
    if (mol->num_basis_functions() == 0) {
        std::cerr << "ERROR: No basis functions!\n";
        return 1;
    }
    
    auto hf = std::make_shared<pyscf::dft::RKS>(mol);
    hf->set_xc_functional("b3lyp");
    auto result = hf->compute();
    
    std::cout << "SCF " << (result.converged ? "converged" : "NOT converged") 
              << " in " << result.iterations << " iterations\n";
    std::cout << "\nTotal energy: " << hf->get_total_energy() << " Eh\n";
    
    return 0;
}
