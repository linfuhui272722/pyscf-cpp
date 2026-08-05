#ifndef PYSCF_CPP_BASIS_PARSER_H
#define PYSCF_CPP_BASIS_PARSER_H

#include "molecule.h"
#include <string>
#include <map>

namespace pyscf {
namespace utils {

// Basis set data structure
struct BasisSet {
    std::string name;
    std::map<int, std::vector<Shell>> shells_by_atom;
};

// Parse Gaussian basis set file
BasisSet parse_basis_file(const std::string& filename);

// Load built-in basis set by name
BasisSet load_basis_set(const std::string& name);

// Available built-in basis sets
std::vector<std::string> available_basis_sets();

// Check if basis set is available
bool has_basis_set(const std::string& name);

// Get basis set for a specific element
std::vector<Shell> get_basis_for_element(const BasisSet& basis, int atomic_number);

// Compute AO overlap matrix S
std::vector<double> compute_overlap_matrix(const Molecule& mol);

// Compute kinetic energy matrix T
std::vector<double> compute_kinetic_matrix(const Molecule& mol);

// Compute nuclear attraction matrix V
std::vector<double> compute_nuclear_matrix(const Molecule& mol);

// Compute one-electron integrals for H_core = T + V
std::vector<double> compute_h_core(const Molecule& mol);

// Evaluate AO basis functions at a point
std::vector<double> eval_ao_at_point(const Molecule& mol, 
                                      double x, double y, double z);

// Evaluate AO basis functions at grid points
void eval_ao_at_grid(const Molecule& mol,
                    const std::vector<double>& coords,
                    std::vector<double>& ao_values);

} // namespace utils
} // namespace pyscf

#endif // PYSCF_CPP_BASIS_PARSER_H
