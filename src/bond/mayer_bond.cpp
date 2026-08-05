/**
 * Mayer Bond Order Analysis Implementation
 * 
 * Computes various bond order indices for molecular bond analysis:
 * - Mayer bond orders
 * - Wiberg bond indices
 * - Valence analysis
 * 
 * Reference:
 * - Mayer, J. Chem. Phys. 85, 5583 (1986) - Mayer bond orders
 * - Wiberg, Tetrahedron 24, 1083 (1968) - Wiberg indices
 * - PySCF lo/nao.py - AO localization
 * - PySCF lo/iao.py - Intrinsic Atomic Orbitals
 */

#include "mayer_bond.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

extern "C" {
// BLAS matrix multiplication
void dgemm_(char* transa, char* transb, int* m, int* n, int* k,
           double* alpha, double* a, int* lda, double* b, int* ldb,
           double* beta, double* c, int* ldc);
}

namespace pyscf {
namespace bond {

// ============================================================================
// Mayer Bond Order Calculation
// ============================================================================

BondOrderMatrix compute_mayer_bond_orders(
    const std::vector<double>& density_matrix,
    const std::vector<double>& overlap_matrix,
    int n_atoms,
    int nao,
    const std::vector<int>& ao_to_atom
) {
    BondOrderMatrix bo(n_atoms);
    
    if (nao <= 0 || n_atoms <= 0) {
        return bo;
    }
    
    // W_AB = sum_{mu in A} sum_{nu in B} P_{mu,nu} * S_{mu,nu}
    for (int mu = 0; mu < nao; ++mu) {
        int atom_i = ao_to_atom[mu];
        if (atom_i < 0 || atom_i >= n_atoms) continue;
        
        for (int nu = 0; nu < nao; ++nu) {
            int atom_j = ao_to_atom[nu];
            if (atom_j < 0 || atom_j >= n_atoms) continue;
            
            // Get density and overlap matrix elements
            double P_munu = 0.0;
            double S_munu = 0.0;
            
            // Density matrix is stored as [i * nao + j]
            if (mu < nao && nu < nao) {
                P_munu = density_matrix[mu * nao + nu];
                S_munu = overlap_matrix[mu * nao + nu];
            }
            
            // Add contribution to bond order
            double contribution = P_munu * S_munu;
            bo.at(atom_i, atom_j) += contribution;
        }
    }
    
    return bo;
}

// ============================================================================
// Wiberg Bond Index Calculation
// ============================================================================

WibergIndex compute_wiberg_index(
    const std::vector<double>& density_matrix,
    int n_atoms,
    int nao,
    const std::vector<int>& ao_to_atom
) {
    WibergIndex wi(n_atoms);
    
    if (nao <= 0 || n_atoms <= 0) {
        return wi;
    }
    
    // WI_AB = sum_{mu in A} sum_{nu in B} P_{mu,nu}^2
    for (int mu = 0; mu < nao; ++mu) {
        int atom_i = ao_to_atom[mu];
        if (atom_i < 0 || atom_i >= n_atoms) continue;
        
        for (int nu = 0; nu < nao; ++nu) {
            int atom_j = ao_to_atom[nu];
            if (atom_j < 0 || atom_j >= n_atoms) continue;
            
            // Get density matrix element
            double P_munu = 0.0;
            if (mu < nao && nu < nao) {
                P_munu = density_matrix[mu * nao + nu];
            }
            
            // Add contribution (P^2 for Wiberg)
            double contribution = P_munu * P_munu;
            wi.at(atom_i, atom_j) += contribution;
        }
    }
    
    return wi;
}

// ============================================================================
// Valence Calculation
// ============================================================================

std::vector<double> compute_valences(const BondOrderMatrix& bo) {
    std::vector<double> valences(bo.n_atoms, 0.0);
    
    for (int i = 0; i < bo.n_atoms; ++i) {
        for (int j = 0; j < bo.n_atoms; ++j) {
            valences[i] += bo.at(i, j);
        }
    }
    
    return valences;
}

// Get expected valence for an element
// Based on common oxidation states and bonding patterns
static double get_expected_valence(int atomic_number) {
    switch (atomic_number) {
        case 1:  return 1.0;   // H
        case 6:  return 4.0;   // C
        case 7:  return 3.0;   // N
        case 8:  return 2.0;   // O
        case 9:  return 1.0;   // F
        case 15: return 3.0;   // P
        case 16: return 2.0;   // S
        case 17: return 1.0;   // Cl
        case 35: return 1.0;   // Br
        case 53: return 1.0;   // I
        default: return 2.0;   // Default
    }
}

std::vector<double> compute_free_valences(
    const BondOrderMatrix& bo,
    const Molecule& mol
) {
    std::vector<double> free_valences(bo.n_atoms, 0.0);
    auto valences = compute_valences(bo);
    
    for (int i = 0; i < bo.n_atoms; ++i) {
        int atomic_num = mol.get_atom(i).atomic_number;
        double expected = get_expected_valence(atomic_num);
        double diff = expected - valences[i];
        free_valences[i] = std::max(0.0, diff * diff);  // Squared difference
        free_valences[i] = std::sqrt(free_valences[i]);  // Root
    }
    
    return free_valences;
}

// ============================================================================
// AO-to-Atom Mapping
// ============================================================================

std::vector<int> compute_ao_to_atom_map(
    std::shared_ptr<Molecule> mol,
    int nao
) {
    std::vector<int> ao_to_atom(nao, -1);
    
    if (!mol || nao <= 0) {
        return ao_to_atom;
    }
    
    // For each atom, we need to know which AOs belong to it
    // This is typically stored in the molecule's shell structure
    // For a simplified version, we'll distribute AOs evenly or
    // use a naive approach based on the fact that AOs are ordered by atom
    
    // Get shells from molecule
    const auto& shells = mol->get_shells();
    
    if (shells.empty()) {
        // Fallback: assume even distribution
        int n_atoms = mol->num_atoms();
        if (n_atoms > 0) {
            int ao_per_atom = nao / n_atoms;
            for (int i = 0; i < nao; ++i) {
                ao_to_atom[i] = std::min(i / ao_per_atom, n_atoms - 1);
            }
        }
        return ao_to_atom;
    }
    
    // Use shell information
    // Each shell has an atom index and angular momentum
    // Count AOs up to each shell's atom
    int current_ao = 0;
    int current_atom = 0;
    
    // Group shells by atom
    std::map<int, int> atom_ao_count;
    for (const auto& shell : shells) {
        int atom_idx = shell.atom_index;
        int l = static_cast<int>(shell.l);
        
        // Number of basis functions in this shell
        // s: 1, p: 3, d: 6, f: 10 (spherical)
        // Contracted Gaussians have multiple primitives but same number of functions
        int n_funcs = (l == 0) ? 1 : (2 * l + 1);
        
        atom_ao_count[atom_idx] += n_funcs;
    }
    
    // Now assign AOs to atoms
    std::map<int, int> atom_ao_start;
    int ao_counter = 0;
    for (int i = 0; i < mol->num_atoms(); ++i) {
        atom_ao_start[i] = ao_counter;
        ao_counter += atom_ao_count[i];
    }
    
    // Fill the mapping
    ao_counter = 0;
    for (const auto& shell : shells) {
        int atom_idx = shell.atom_index;
        int l = static_cast<int>(shell.l);
        int n_funcs = (l == 0) ? 1 : (2 * l + 1);
        
        for (int i = 0; i < n_funcs && ao_counter < nao; ++i) {
            ao_to_atom[ao_counter++] = atom_idx;
        }
    }
    
    return ao_to_atom;
}

// ============================================================================
// Bond Identification
// ============================================================================

std::vector<std::tuple<int, int, double>> identify_bonds(
    const BondOrderMatrix& bo,
    double threshold
) {
    std::vector<std::tuple<int, int, double>> bonds;
    
    for (int i = 0; i < bo.n_atoms; ++i) {
        for (int j = i + 1; j < bo.n_atoms; ++j) {
            double order = bo.at(i, j);
            if (order > threshold) {
                bonds.emplace_back(i, j, order);
            }
        }
    }
    
    // Sort by bond order (descending)
    std::sort(bonds.begin(), bonds.end(),
        [](const auto& a, const auto& b) {
            return std::get<2>(a) > std::get<2>(b);
        });
    
    return bonds;
}

// ============================================================================
// Complete Bond Analysis
// ============================================================================

BondAnalysisResult analyze_bonds(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol
) {
    BondAnalysisResult result;
    
    int n_atoms = mol->num_atoms();
    result.bo_matrix.resize(n_atoms);
    result.wi_matrix.resize(n_atoms);
    result.valences.resize(n_atoms);
    result.free_valences.resize(n_atoms);
    
    // Get density and overlap matrices
    auto density_matrix = scf->get_density_matrix();
    auto overlap_matrix = scf->get_overlap_matrix();
    
    int nao = scf->get_num_orbitals();
    
    // Compute AO-to-atom mapping
    auto ao_to_atom = compute_ao_to_atom_map(mol, nao);
    
    // Compute bond orders
    result.bo_matrix = compute_mayer_bond_orders(
        density_matrix, overlap_matrix, n_atoms, nao, ao_to_atom);
    
    // Compute Wiberg indices
    result.wi_matrix = compute_wiberg_index(
        density_matrix, n_atoms, nao, ao_to_atom);
    
    // Compute valences
    result.valences = compute_valences(result.bo_matrix);
    result.free_valences = compute_free_valences(result.bo_matrix, *mol);
    
    // Identify bonds
    result.bonds = identify_bonds(result.bo_matrix, 0.1);
    
    return result;
}

// ============================================================================
// Print Functions
// ============================================================================

void print_bond_analysis(
    const BondAnalysisResult& result,
    const Molecule& mol
) {
    int n_atoms = mol.num_atoms();
    
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "                        BOND ORDER ANALYSIS (Mayer)\n";
    std::cout << "================================================================================\n\n";
    
    // Header
    std::cout << "Atom coordinates (Bohr):\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    for (int i = 0; i < n_atoms; ++i) {
        auto atom = mol.get_atom(i);
        std::cout << "  " << std::setw(2) << element_symbol(atom.atomic_number) << std::setw(2) << (i + 1);
        std::cout << "  " << std::fixed << std::setprecision(6);
        std::cout << "  " << std::setw(12) << atom.x;
        std::cout << "  " << std::setw(12) << atom.y;
        std::cout << "  " << std::setw(12) << atom.z;
        std::cout << "\n";
    }
    
    std::cout << "\n";
    std::cout << "Mayer Bond Order Matrix:\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    // Column headers
    std::cout << "        ";
    for (int j = 0; j < n_atoms; ++j) {
        std::cout << "  " << std::setw(8) << element_symbol(mol.get_atom(j).atomic_number) << std::setw(2) << (j+1);
    }
    std::cout << "\n";
    
    // Matrix rows
    for (int i = 0; i < n_atoms; ++i) {
        std::cout << "  " << std::setw(2) << element_symbol(mol.get_atom(i).atomic_number) << std::setw(2) << (i+1);
        for (int j = 0; j < n_atoms; ++j) {
            std::cout << "  " << std::fixed << std::setprecision(4);
            std::cout << "  " << std::setw(8) << result.bo_matrix.at(i, j);
        }
        std::cout << "\n";
    }
    
    std::cout << "\n";
    std::cout << "Valence Analysis:\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << "  " << std::setw(12) << "Atom";
    std::cout << "  " << std::setw(12) << "Valence";
    std::cout << "  " << std::setw(12) << "Free Valence";
    std::cout << "\n";
    for (int i = 0; i < n_atoms; ++i) {
        std::cout << "  " << std::setw(2) << element_symbol(mol.get_atom(i).atomic_number) << std::setw(2) << (i + 1);
        std::cout << "  " << std::fixed << std::setprecision(4);
        std::cout << "  " << std::setw(12) << result.valences[i];
        std::cout << "  " << std::setw(12) << result.free_valences[i];
        std::cout << "\n";
    }
    
    std::cout << "\n";
    std::cout << "Bond List (order > 0.1):\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    if (result.bonds.empty()) {
        std::cout << "  No bonds found\n";
    } else {
        std::cout << "  " << std::setw(8) << "Bond";
        std::cout << "  " << std::setw(12) << "Type";
        std::cout << "  " << std::setw(12) << "Order";
        std::cout << "\n";
        for (const auto& bond : result.bonds) {
            int i = std::get<0>(bond);
            int j = std::get<1>(bond);
            double order = std::get<2>(bond);
            
            std::string bond_type;
            if (order > 2.5) bond_type = "Triple";
            else if (order > 1.5) bond_type = "Double";
            else if (order > 0.5) bond_type = "Single";
            else bond_type = "Weak";
            
            std::cout << "  " << std::setw(2) << element_symbol(mol.get_atom(i).atomic_number) << std::setw(2) << (i + 1);
            std::cout << "-";
            std::cout << std::setw(2) << element_symbol(mol.get_atom(j).atomic_number) << std::setw(2) << (j + 1);
            std::cout << "  " << std::setw(12) << bond_type;
            std::cout << "  " << std::fixed << std::setprecision(4);
            std::cout << "  " << std::setw(12) << order;
            std::cout << "\n";
        }
    }
    
    std::cout << "\n";
    std::cout << "================================================================================\n\n";
}

// ============================================================================
// Popov Bond Orders (Stub)
// ============================================================================

PopovBondOrders compute_popov_bond_orders(
    const std::vector<double>& density_matrix,
    const std::vector<double>& overlap_matrix,
    const std::vector<double>& orthog_matrix,
    int n_atoms,
    int nao,
    const std::vector<int>& ao_to_atom
) {
    // Popov bond orders use orthogonalized basis
    // For now, just return the Mayer bond orders
    PopovBondOrders result;
    result.n_atoms = n_atoms;
    result.data.resize(n_atoms * n_atoms, 0.0);
    result.valences.resize(n_atoms, 0.0);
    
    auto bo = compute_mayer_bond_orders(
        density_matrix, overlap_matrix, n_atoms, nao, ao_to_atom);
    
    for (int i = 0; i < n_atoms; ++i) {
        for (int j = 0; j < n_atoms; ++j) {
            result.at(i, j) = bo.at(i, j);
        }
    }
    
    for (int i = 0; i < n_atoms; ++i) {
        result.valences[i] = bo.at(i, i);
    }
    
    return result;
}

} // namespace bond
} // namespace pyscf
