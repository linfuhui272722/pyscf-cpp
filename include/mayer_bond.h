#ifndef PYSCF_CPP_MAYER_BOND_H
#define PYSCF_CPP_MAYER_BOND_H

#include "molecule.h"
#include "scf.h"
#include <vector>
#include <map>
#include <memory>
#include <string>

namespace pyscf {
namespace bond {

// Mayer Bond Order Matrix
// Reference: Mayer, J. Chem. Phys. 85, 5583 (1986)
struct BondOrderMatrix {
    int n_atoms;
    std::vector<double> data;  // [i * n_atoms + j]
    
    BondOrderMatrix() : n_atoms(0) {}
    BondOrderMatrix(int n) : n_atoms(n), data(n * n, 0.0) {}
    
    double& at(int i, int j) { return data[i * n_atoms + j]; }
    const double& at(int i, int j) const { return data[i * n_atoms + j]; }
    
    void resize(int n) {
        n_atoms = n;
        data.assign(n * n, 0.0);
    }
};

// Wiberg Bond Index
// Reference: Wiberg, Tetrahedron 24, 1083 (1968)
struct WibergIndex {
    int n_atoms;
    std::vector<double> data;
    
    WibergIndex() : n_atoms(0) {}
    WibergIndex(int n) : n_atoms(n), data(n * n, 0.0) {}
    
    double& at(int i, int j) { return data[i * n_atoms + j]; }
    const double& at(int i, int j) const { return data[i * n_atoms + j]; }
    
    void resize(int n) {
        n_atoms = n;
        data.assign(n * n, 0.0);
    }
};

// Bond analysis result
struct BondAnalysisResult {
    BondOrderMatrix bo_matrix;           // Mayer bond orders
    WibergIndex wi_matrix;               // Wiberg indices
    std::vector<double> valences;         // Total valence per atom
    std::vector<double> free_valences;   // Free valence per atom
    std::vector<std::tuple<int, int, double>> bonds;  // (i, j, order) sorted
    std::vector<std::tuple<int, int, double>> pi_bonds;  // pi bond contributions
};

// Compute Mayer bond orders
// W_AB = sum_{mu in A} sum_{nu in B} (P_{mu,nu} * S_{mu,nu})
// where P is the density matrix and S is the overlap matrix
BondOrderMatrix compute_mayer_bond_orders(
    const std::vector<double>& density_matrix,
    const std::vector<double>& overlap_matrix,
    int n_atoms,
    int nao,
    const std::vector<int>& ao_to_atom  // Maps AO index to atom index
);

// Compute Wiberg bond indices
// WI_AB = sum_{mu in A} sum_{nu in B} P_{mu,nu}^2
WibergIndex compute_wiberg_index(
    const std::vector<double>& density_matrix,
    int n_atoms,
    int nao,
    const std::vector<int>& ao_to_atom
);

// Compute total valence for each atom
// V_A = sum_B W_AB
std::vector<double> compute_valences(
    const BondOrderMatrix& bo
);

// Compute free valence (difference from maximum)
// F_A = sqrt( (N_valence[A] - V_A)^2 )
// For main group: max valence ~ 2 * n_unpaired_electrons
std::vector<double> compute_free_valences(
    const BondOrderMatrix& bo,
    const Molecule& mol
);

// Complete bond analysis
BondAnalysisResult analyze_bonds(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol
);

// Print bond analysis in a formatted way
void print_bond_analysis(
    const BondAnalysisResult& result,
    const Molecule& mol
);

// Identify bonds based on threshold
std::vector<std::tuple<int, int, double>> identify_bonds(
    const BondOrderMatrix& bo,
    double threshold = 0.1
);

// Compute AO-to-atom mapping
// Returns vector where ao_to_atom[mu] = atom index containing AO mu
std::vector<int> compute_ao_to_atom_map(
    std::shared_ptr<Molecule> mol,
    int nao
);

// Popov bond orders (for delocalized systems)
// Reference: Popov, J. Chem. Phys. 83, 2947 (1985)
struct PopovBondOrders {
    int n_atoms;
    std::vector<double> data;
    std::vector<double> valences;
    
    double& at(int i, int j) { return data[i * n_atoms + j]; }
    const double& at(int i, int j) const { return data[i * n_atoms + j]; }
};

PopovBondOrders compute_popov_bond_orders(
    const std::vector<double>& density_matrix,
    const std::vector<double>& overlap_matrix,
    const std::vector<double>& orthog_matrix,
    int n_atoms,
    int nao,
    const std::vector<int>& ao_to_atom
);

} // namespace bond
} // namespace pyscf

#endif // PYSCF_CPP_MAYER_BOND_H
