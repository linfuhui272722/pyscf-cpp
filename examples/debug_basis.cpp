#include <iostream>
#include <memory>
#include "pyscf_cpp.h"

using namespace pyscf;

int main() {
    auto mol = std::make_shared<Molecule>();
    mol->add_atom(1, 0.0, 0.0, 0.0);
    mol->add_atom(1, 1.389, 0.0, 0.0);
    
    mol->set_basis("6-31g*");
    
    std::cout << "Total basis functions: " << mol->num_basis_functions() << std::endl;
    
    return 0;
}
