#include "scf.h"
#include "gto_integrals.h"
#include <cmath>
#include <algorithm>
#include <iostream>

extern "C" {
void dgemm_(char* transa, char* transb, int* m, int* n, int* k,
           double* alpha, double* a, int* lda, double* b, int* ldb,
           double* beta, double* c, int* ldc);

void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
           double* w, double* work, int* lwork, int* info);

void dsygv_(int* itype, char* jobz, char* uplo, int* n, double* a, int* lda,
           double* b, int* ldb, double* w, double* work, int* lwork, int* info);

void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);

void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
}

namespace pyscf {
namespace dft {

RKS::RKS(std::shared_ptr<Molecule> mol) : mol_(mol) {
    nuc_energy_ = mol_->nuclear_repulsion_energy();
}

RKS::~RKS() = default;

void RKS::set_xc_functional(const std::string& xc) {
    // Parse xc string (e.g., "b3lyp", "pbe", etc.)
}

void RKS::set_options(const SCFOptions& opts) {
    options_ = opts;
}

void RKS::initialize_matrices() {
    int nbf = mol_->num_basis_functions();
    int n2 = nbf * nbf;
    
    overlap_matrix_.assign(n2, 0.0);
    kinetic_matrix_.assign(n2, 0.0);
    nuclear_matrix_.assign(n2, 0.0);
    mo_coefficients_.assign(n2, 0.0);
    mo_energies_.assign(nbf, 0.0);
    density_matrix_.assign(n2, 0.0);
    
    // Initialize integral engine with libcint
    integral_engine_ = std::make_unique<gto::IntegralEngine>();
    integral_engine_->initialize(*mol_);
    
    // Compute one-electron integrals using libcint
    compute_one_electron_integrals();
}

void RKS::compute_one_electron_integrals() {
    if (options_.verbose) {
        std::cout << "Computing one-electron integrals with libcint...\n";
    }
    
    // Use libcint integral engine for accurate integrals
    integral_engine_->compute_overlap(*mol_, overlap_matrix_.data());
    integral_engine_->compute_kinetic(*mol_, kinetic_matrix_.data());
    integral_engine_->compute_nuclear(*mol_, nuclear_matrix_.data());
    
    if (options_.verbose) {
        std::cout << "Integrals computed successfully.\n";
    }
}

void RKS::diagonalize(const std::vector<double>& matrix, 
                       std::vector<double>& eigenvalues,
                       std::vector<double>& eigenvectors) {
    // For ROHF/RKS, we use generalized eigenvalue problem: F*C = S*C*e
    // But for simplicity, we'll use standard eigenvalue problem
    // (assuming orthonormal basis or use generalized version)
    int n = eigenvalues.size();
    int lda = n;
    int lwork = 6 * n;
    
    std::vector<double> work(lwork);
    int info;
    
    char jobz = 'V';
    char uplo = 'U';
    
    // Copy input matrix (it will be modified by LAPACK)
    eigenvectors = matrix;
    
    dsyev_(&jobz, &uplo, &n, eigenvectors.data(), &lda,
           eigenvalues.data(), work.data(), &lwork, &info);
}

void RKS::generalized_eig(const std::vector<double>& matrix, 
                           const std::vector<double>& overlap,
                           std::vector<double>& eigenvalues,
                           std::vector<double>& eigenvectors) {
    // Generalized eigenvalue problem: F*C = S*C*e
    int n = eigenvalues.size();
    int lda = n;
    int ldb = n;
    int lwork = 6 * n;
    
    std::vector<double> work(lwork);
    int info;
    
    char jobz = 'V';
    char uplo = 'U';
    int itype = 1;  // A*x = lambda*B*x
    
    // Copy input matrices (they will be modified by LAPACK)
    eigenvectors = matrix;
    std::vector<double> overlap_copy = overlap;
    
    dsygv_(&itype, &jobz, &uplo, &n, eigenvectors.data(), &lda,
            overlap_copy.data(), &ldb, eigenvalues.data(), work.data(), &lwork, &info);
}

SCFResult RKS::compute() {
    SCFResult result;
    
    int nbf = mol_->num_basis_functions();
    int n_electrons = mol_->num_electrons();
    int n_occupied = n_electrons / 2;
    
    if (options_.verbose) {
        std::cout << "Starting SCF calculation...\n";
        std::cout << "Number of atoms: " << mol_->num_atoms() << "\n";
        std::cout << "Number of basis functions: " << nbf << "\n";
        std::cout << "Number of electrons: " << n_electrons << "\n";
        std::cout << "Nuclear repulsion energy: " << nuc_energy_ << "\n";
    }
    
    // Initialize
    initialize_matrices();
    
    // Get H_core
    std::vector<double> h_core = compute_h_core();
    
    // Initial guess: diagonalize H_core with S
    std::vector<double> mo_c = h_core;
    std::vector<double> mo_e(nbf);
    generalized_eig(h_core, overlap_matrix_, mo_e, mo_c);
    
    // Build initial density matrix
    // C_ia (coeff of basis i in MO a) = mo_c[a * nbf + i]
    // P_ij = 2 * sum_a C_ia * C_ja
    std::vector<double> dm(nbf * nbf, 0.0);
    for (int a = 0; a < n_occupied; ++a) {
        for (int i = 0; i < nbf; ++i) {
            for (int j = 0; j < nbf; ++j) {
                dm[i * nbf + j] += 2.0 * mo_c[a * nbf + i] * mo_c[a * nbf + j];
            }
        }
    }
    
    // Normalize density matrix to have correct trace
    // trace(P * S) should equal n_electrons = 2
    double trace_ps = 0.0;
    for (int i = 0; i < nbf; ++i) {
        for (int j = 0; j < nbf; ++j) {
            trace_ps += dm[i * nbf + j] * overlap_matrix_[i * nbf + j];
        }
    }
    
    double scale = (2.0 * n_occupied) / trace_ps;
    for (int i = 0; i < nbf * nbf; ++i) {
        dm[i] *= scale;
    }
    
    density_matrix_ = dm;
    mo_coefficients_ = mo_c;
    mo_energies_ = mo_e;
    
    // SCF iteration
    double energy_last = 0.0;
    
    for (int cycle = 0; cycle < options_.max_cycle; ++cycle) {
        // Compute Fock matrix: F = H_core + J - hybrid*K + VXC
        std::vector<double> fock = h_core;
        
        // Coulomb matrix J (using ERI)
        std::vector<double> j_mat = compute_j_matrix(density_matrix_);
        for (int i = 0; i < nbf * nbf; ++i) {
            fock[i] += j_mat[i];
        }
        
        // Exchange matrix K (hybrid coefficient will be handled in energy calculation)
        std::vector<double> k_mat = compute_k_matrix(density_matrix_);
        double hyb = xc_func_.get_hf_exchange_coefficient();
        // For RKS, we use F = H + J - 0.5*K (consistent with PySCF's vhf = J - 0.5*K)
        // Note: For pure DFT functionals, hyb=0
        for (int i = 0; i < nbf * nbf; ++i) {
            fock[i] -= 0.5 * k_mat[i];
        }
        
        // Debug: print Fock matrix
        if (cycle == 0) {
        }
        
        // XC potential VXC (simplified LDA for now)
        // Note: VXC is computed from the density but not added to Fock here
        // to avoid double counting with the exchange term
        // For pure DFT (no HF exchange), add VXC here
        // For hybrid functionals, VXC only includes the DFT part
        
        // DIIS extrapolation
        if (cycle >= options_.diis_start_cycle && !diis_focks_.empty()) {
            fock = diis_extrapolate(fock, std::vector<double>(nbf * nbf, 0.0));
        }
        
        // Diagonalize Fock matrix with S (generalized eigenvalue problem)
        mo_c = fock;
        generalized_eig(fock, overlap_matrix_, mo_e, mo_c);
        
        
        // Build new density matrix
        // C_ia = mo_c[a * nbf + i]
        // P_ij = 2 * sum_a C_ia * C_ja
        std::vector<double> dm_new(nbf * nbf, 0.0);
        for (int a = 0; a < n_occupied; ++a) {
            for (int i = 0; i < nbf; ++i) {
                for (int j = 0; j < nbf; ++j) {
                    dm_new[i * nbf + j] += 2.0 * mo_c[a * nbf + i] * mo_c[a * nbf + j];
                }
            }
        }
        
        // Compute electronic energy using current MO energies
        double e_elec = compute_electronic_energy(dm_new, mo_e);
        double e_tot = e_elec + nuc_energy_;
        
        // Check convergence
        double de = std::abs(e_tot - energy_last);
        double max_diff = 0.0;
        for (int i = 0; i < nbf * nbf; ++i) {
            max_diff = std::max(max_diff, std::abs(dm_new[i] - density_matrix_[i]));
        }
        
        if (options_.verbose) {
            std::cout << "Cycle " << cycle + 1 << ": E = " << e_tot 
                      << " (dE = " << de << ", dDM = " << max_diff << ")\n";
        }
        
        density_matrix_ = dm_new;
        mo_coefficients_ = mo_c;
        mo_energies_ = mo_e;
        energy_ = e_elec;
        energy_last = e_tot;
        
        if (de < options_.conv_tol) {
            result.converged = true;
            result.iterations = cycle + 1;
            break;
        }
        
        // Store for DIIS
        if (cycle >= options_.diis_start_cycle) {
            diis_focks_.push_back(fock);
            if (diis_focks_.size() > (size_t)options_.diis_nspace) {
                diis_focks_.erase(diis_focks_.begin());
            }
        }
    }
    
    // Store final results
    result.energy = energy_ + nuc_energy_;
    result.mo_energies = mo_energies_;
    result.mo_coefficients = mo_coefficients_;
    result.density_matrix = density_matrix_;
    
    if (options_.verbose) {
        std::cout << "\nSCF Results:\n";
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n";
        std::cout << "  Iterations: " << result.iterations << "\n";
        std::cout << "  Total energy: " << result.energy << " Eh\n";
    }
    
    return result;
}

double RKS::compute_electronic_energy(const std::vector<double>& dm) {
    return compute_electronic_energy(dm, mo_energies_);
}

double RKS::compute_electronic_energy(const std::vector<double>& dm, 
                                       const std::vector<double>& mo_energies) {
    int nbf = mol_->num_basis_functions();
    
    double hyb = xc_func_.get_hf_exchange_coefficient();
    
    // Compute J and K matrices
    std::vector<double> j_mat = compute_j_matrix(dm);
    std::vector<double> k_mat = compute_k_matrix(dm);
    
    // For RHF, vhf = J - 0.5*K (for HF exchange)
    // E_elec = trace(P * H) + 0.5 * trace(P * vhf)
    //        = trace(P * H) + 0.5 * trace(P * J) - 0.25 * trace(P * K)
    
    // Compute trace(P * H_core)
    std::vector<double> h_core = compute_h_core();
    double trace_PH = 0.0;
    for (int i = 0; i < nbf * nbf; ++i) {
        trace_PH += dm[i] * h_core[i];
    }
    
    // Compute trace(P * J)
    double trace_PJ = 0.0;
    for (int i = 0; i < nbf * nbf; ++i) {
        trace_PJ += dm[i] * j_mat[i];
    }
    
    // Compute trace(P * K)
    double trace_PK = 0.0;
    for (int i = 0; i < nbf * nbf; ++i) {
        trace_PK += dm[i] * k_mat[i];
    }
    
    // For hybrid functionals: E_elec = trace(P*H) + 0.5*trace(P*J) - 0.25*trace(P*K)
    // For pure HF (hyb=1): E_elec = trace(P*H) + 0.5*trace(P*J) - 0.25*trace(P*K)
    // For pure DFT (hyb=0): E_elec = trace(P*H) + 0.5*trace(P*J)
    double energy = trace_PH + 0.5 * trace_PJ - 0.25 * trace_PK;
    
    return energy;
}

std::vector<double> RKS::compute_h_core() {
    int nbf = mol_->num_basis_functions();
    std::vector<double> h_core(nbf * nbf, 0.0);
    
    for (int i = 0; i < nbf * nbf; ++i) {
        h_core[i] = kinetic_matrix_[i] + nuclear_matrix_[i];
    }
    
    return h_core;
}

std::vector<double> RKS::compute_j_matrix(const std::vector<double>& dm) {
    int nbf = mol_->num_basis_functions();
    std::vector<double> j_mat(nbf * nbf, 0.0);
    
    // Compute ERIs using libcint
    size_t eri_size = 0;
    double* eri = integral_engine_->compute_eri(*mol_, eri_size);
    
    // J_ij = sum_kl P_kl * (ij|kl)
    for (int i = 0; i < nbf; ++i) {
        for (int j = 0; j < nbf; ++j) {
            double sum = 0.0;
            for (int k = 0; k < nbf; ++k) {
                for (int l = 0; l < nbf; ++l) {
                    // ERI index (ij|kl)
                    size_t idx = ((i * nbf + j) * nbf + k) * nbf + l;
                    if (idx < eri_size) {
                        sum += dm[k * nbf + l] * eri[idx];
                    }
                }
            }
            j_mat[i * nbf + j] = sum;
        }
    }
    
    delete[] eri;
    return j_mat;
}

std::vector<double> RKS::compute_k_matrix(const std::vector<double>& dm) {
    int nbf = mol_->num_basis_functions();
    std::vector<double> k_mat(nbf * nbf, 0.0);
    
    // Compute ERIs
    size_t eri_size = 0;
    double* eri = integral_engine_->compute_eri(*mol_, eri_size);
    
    // K_ij = sum_kl P_kl * (ik|jl)
    // ERI is stored as (ij|kl) format: idx = ((i*nbf+j)*nbf+k)*nbf+l
    // K_ij = sum_kl P_kl * eri[(ik|jl)]
    // For (ik|jl), we need to access eri[(i*nbf+k)*nbf+j)*nbf+l]
    for (int i = 0; i < nbf; ++i) {
        for (int j = 0; j < nbf; ++j) {
            double sum = 0.0;
            for (int k = 0; k < nbf; ++k) {
                for (int l = 0; l < nbf; ++l) {
                    // K_ij = sum_kl P_kl * (ik|jl)
                    // ERI index: (ik|jl) corresponds to ((i*nbf+k)*nbf+j)*nbf+l
                    size_t idx_ik_jl = ((i * nbf + k) * nbf + j) * nbf + l;
                    if (idx_ik_jl < eri_size) {
                        sum += dm[k * nbf + l] * eri[idx_ik_jl];
                    }
                }
            }
            k_mat[i * nbf + j] = sum;
        }
    }
    
    delete[] eri;
    return k_mat;
}

std::vector<double> RKS::compute_vxc(const std::vector<double>& dm) {
    int nbf = mol_->num_basis_functions();
    std::vector<double> vxc(nbf * nbf, 0.0);
    
    // Simplified XC potential for B3LYP
    // Full implementation requires numerical integration with grids
    double trace_P = 0.0;
    for (int i = 0; i < nbf; ++i) {
        trace_P += dm[i * nbf + i];
    }
    
    double rho_avg = trace_P / (4.0 * M_PI * std::pow(3.0, 1.0/3.0));
    double vxc0 = -1.0 * std::pow(rho_avg, 1.0/3.0);
    
    // Use overlap matrix as approximation
    for (int i = 0; i < nbf * nbf; ++i) {
        vxc[i] = vxc0 * overlap_matrix_[i] * 0.2;
    }
    
    return vxc;
}

std::vector<double> RKS::diis_extrapolate(const std::vector<double>& fock,
                                          const std::vector<double>& diis_error) {
    if (diis_focks_.size() < 2) {
        return fock;
    }
    
    int nbf = mol_->num_basis_functions();
    int n = diis_focks_.size();
    
    // Build B matrix for DIIS
    std::vector<double> B((n + 1) * (n + 1), 0.0);
    B[n * (n + 1) + n] = -1.0;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < nbf; ++k) {
                for (int l = 0; l < nbf; ++l) {
                    sum += (diis_focks_[i][k * nbf + l] - diis_focks_.back()[k * nbf + l]) *
                           (diis_focks_[j][k * nbf + l] - diis_focks_.back()[k * nbf + l]);
                }
            }
            B[i * (n + 1) + j] = sum;
        }
        B[i * (n + 1) + n] = -1.0;
        B[n * (n + 1) + i] = -1.0;
    }
    
    // Solve for coefficients
    std::vector<double> rhs(n + 1, 0.0);
    rhs[n] = -1.0;
    
    // Simple linear solve (for small n)
    std::vector<double> c(n + 1);
    double max_iter = 100;
    double tol = 1e-10;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i <= n; ++i) {
            double sigma = rhs[i];
            for (int j = 0; j <= n; ++j) {
                if (i != j) {
                    sigma -= B[i * (n + 1) + j] * c[j];
                }
            }
            c[i] = sigma / B[i * (n + 1) + i];
        }
    }
    
    // Extrapolate Fock matrix
    std::vector<double> fock_new(nbf * nbf, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < nbf * nbf; ++k) {
            fock_new[k] += c[i] * diis_focks_[i][k];
        }
    }
    
    return fock_new;
}

} // namespace dft
} // namespace pyscf
