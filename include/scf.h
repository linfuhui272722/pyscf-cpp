#ifndef PYSCF_CPP_SCF_H
#define PYSCF_CPP_SCF_H

#include "molecule.h"
#include "lebedev_grid.h"
#include "b3lyp_xc.h"
#include <vector>
#include <memory>

namespace pyscf {

// Forward declaration for IntegralEngine
namespace gto { class IntegralEngine; }

namespace dft {

// SCF convergence options
struct SCFOptions {
    double conv_tol = 1e-8;
    double max_cycle = 100;
    double diis_start_cycle = 3;
    int diis_nspace = 8;
    bool verbose = true;
};

// SCF result
struct SCFResult {
    double energy;
    std::vector<double> mo_energies;
    std::vector<double> mo_coefficients;
    std::vector<double> density_matrix;
    int iterations;
    bool converged;
};

// Restricted Kohn-Sham SCF
class RKS {
public:
    RKS(std::shared_ptr<Molecule> mol);
    ~RKS();
    
    void set_xc_functional(const std::string& xc);
    void set_options(const SCFOptions& opts);
    SCFResult compute();
    
    // Get energy components
    double get_electronic_energy() const { return energy_; }
    double get_nuclear_repulsion() const { return nuc_energy_; }
    double get_total_energy() const { return energy_ + nuc_energy_; }
    
    // Get MO coefficients
    const std::vector<double>& get_mo_coefficients() const { return mo_coefficients_; }

    // Get density matrix
    const std::vector<double>& get_density_matrix() const { return density_matrix_; }

    // Get overlap matrix
    const std::vector<double>& get_overlap_matrix() const { return overlap_matrix_; }

    // Get number of orbitals (NAOs)
    int get_num_orbitals() const { return mol_ ? mol_->num_basis_functions() : 0; }
    
    // AO evaluation at grid points
    std::vector<double> eval_ao(const std::vector<dft::GridPoint>& grid);
    
private:
    std::shared_ptr<Molecule> mol_;
    B3LYPFunctional xc_func_;
    SCFOptions options_;
    
    // Current state
    double energy_ = 0.0;
    double nuc_energy_ = 0.0;
    std::vector<double> mo_coefficients_;
    std::vector<double> mo_energies_;
    std::vector<double> density_matrix_;
    
    // One-electron integrals
    std::vector<double> overlap_matrix_;
    std::vector<double> kinetic_matrix_;
    std::vector<double> nuclear_matrix_;
    
    // Integral engine using libcint
    std::unique_ptr<gto::IntegralEngine> integral_engine_;
    
    // DIIS storage
    std::vector<std::vector<double>> diis_focks_;
    std::vector<std::vector<double>> diis_errors_;
    
    void initialize_matrices();
    void compute_one_electron_integrals();
    double compute_electronic_energy(const std::vector<double>& dm);
    double compute_electronic_energy(const std::vector<double>& dm, 
                                     const std::vector<double>& mo_energies);
    std::vector<double> compute_h_core();
    std::vector<double> compute_j_matrix(const std::vector<double>& dm);
    std::vector<double> compute_k_matrix(const std::vector<double>& dm);
    std::vector<double> compute_vxc(const std::vector<double>& dm);
    
    // DIIS extrapolation
    std::vector<double> diis_extrapolate(const std::vector<double>& fock,
                                          const std::vector<double>& error);
    
    // Linear algebra helpers
    void diagonalize(const std::vector<double>& matrix, 
                     std::vector<double>& eigenvalues,
                     std::vector<double>& eigenvectors);
    void generalized_eig(const std::vector<double>& matrix, 
                        const std::vector<double>& overlap,
                        std::vector<double>& eigenvalues,
                        std::vector<double>& eigenvectors);
    void dgemm(char transa, char transb, int m, int n, int k,
              double alpha, const std::vector<double>& a, int lda,
              const std::vector<double>& b, int ldb,
              double beta, std::vector<double>& c, int ldc);
};

} // namespace dft
} // namespace pyscf

#endif // PYSCF_CPP_SCF_H
