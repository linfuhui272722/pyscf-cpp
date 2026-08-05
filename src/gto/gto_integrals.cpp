#include "gto_integrals.h"
#include "molecule.h"
#include <cmath>
#include <cstring>
#include <iostream>

// libcint headers
extern "C" {
#include <cint.h>
}

// Helper functions for GTO normalization
namespace {

// Compute (2n-1)!! = 1 * 3 * 5 * ... * (2n-1)
double double_factorial(int n) {
    if (n <= 0) return 1.0;
    double result = 1.0;
    for (int i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

// Normalize contraction coefficient for libcint
// libcint expects normalized primitive GTO coefficients
double normalize_gto(double alpha, int l, int m, int n, double c) {
    int lmn = l + m + n;
    const double PI = 3.14159265358979323846;
    
    // Normalization constant for primitive GTO:
    // N(l,m,n,alpha) = sqrt(2^(2l+2m+2n+3) * alpha^(l+m+n+3/2) / 
    //                  (pi^(3/2) * (2l+2m+2n+1)!!))
    double two_alpha = 2.0 * alpha;
    double df = double_factorial(2 * lmn + 1);
    double norm = std::sqrt(
        std::pow(two_alpha, lmn + 1.5) / 
        (std::pow(PI, 1.5) * df)
    );
    
    return c * norm;
}

} // anonymous namespace

// Use libcint function declarations
extern "C" {
// One-electron integrals
FINT cint1e_ovlp_cart(double *buf, FINT *shls,
                       FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
FINT cint1e_kin_cart(double *buf, FINT *shls,
                      FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
FINT cint1e_nuc_cart(double *buf, FINT *shls,
                     FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

// Two-electron integrals
FINT cint2e_cart(double *buf, FINT *shls,
                 FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                 CINTOpt *opt);

// Utility functions
FINT CINTcgtos_cart(FINT bas_id, const FINT *bas);
FINT CINTtot_cgto_cart(const FINT *bas, FINT nbas);
void CINTshells_cart_offset(FINT ao_loc[], const FINT *bas, FINT nbas);
void CINTinit_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env);
void CINTdel_optimizer(CINTOpt **opt);
}

namespace pyscf {
namespace gto {

// Internal implementation using libcint
struct IntegralEngine::Impl {
    // libcint data structures
    std::vector<FINT> atm;      // atom data
    std::vector<FINT> bas;     // basis data
    std::vector<double> env;   // environment data
    CINTOpt* opt = nullptr;
    int natm = 0;
    int nbas = 0;
    int nbf = 0;               // number of basis functions
    std::vector<FINT> ao_loc;  // AO location index
    
    ~Impl() {
        if (opt) {
            CINTdel_optimizer(&opt);
        }
    }
};

IntegralEngine::IntegralEngine() : pimpl_(new Impl) {}

IntegralEngine::~IntegralEngine() = default;

void IntegralEngine::initialize(const Molecule& mol) {
    auto& impl = *pimpl_;
    
    // Clear previous data
    impl.atm.clear();
    impl.bas.clear();
    impl.env.clear();
    impl.natm = mol.num_atoms();
    impl.nbas = mol.num_shells();
    
    // Environment: PTR_EXPCUTOFF needs to be set for integral prescreening
    // Value is approximately ln(threshold) for the cutoff
    impl.env.resize(20);  // PTR_ENV_START = 20
    impl.env[PTR_EXPCUTOFF] = -100.0;  // Cutoff threshold (very permissive)
    
    size_t env_offset = impl.env.size();
    
    // Add atom data
    impl.atm.resize(ATM_SLOTS * impl.natm);
    for (int i = 0; i < impl.natm; ++i) {
        const auto& atom = mol.get_atom(i);
        impl.atm[ATM_SLOTS * i + CHARGE_OF] = atom.atomic_number;
        impl.atm[ATM_SLOTS * i + PTR_COORD] = env_offset;
        impl.atm[ATM_SLOTS * i + NUC_MOD_OF] = POINT_NUC;
        impl.atm[ATM_SLOTS * i + 3] = 0;  // External charge pointer
        impl.atm[ATM_SLOTS * i + 4] = 0;  // Reserved
        impl.atm[ATM_SLOTS * i + 5] = 0;  // Reserved
        
        // Add coordinates to environment
        impl.env.push_back(atom.x);
        impl.env.push_back(atom.y);
        impl.env.push_back(atom.z);
        env_offset = impl.env.size();
    }
    
    // Add basis set data
    impl.bas.resize(BAS_SLOTS * impl.nbas);
    for (int i = 0; i < impl.nbas; ++i) {
        const auto& shell = mol.get_shell(i);
        impl.bas[BAS_SLOTS * i + ATOM_OF] = shell.atom_index;
        impl.bas[BAS_SLOTS * i + ANG_OF] = shell.l + shell.m + shell.n;  // Total angular momentum
        impl.bas[BAS_SLOTS * i + NPRIM_OF] = shell.exponents.size();
        impl.bas[BAS_SLOTS * i + NCTR_OF] = 1;  // One contraction
        impl.bas[BAS_SLOTS * i + KAPPA_OF] = 0;  // 0 = Cartesian, 1 = spherical
        
        size_t exp_offset = impl.env.size();
        impl.bas[BAS_SLOTS * i + PTR_EXP] = exp_offset;
        
        // Add exponents
        for (double exp : shell.exponents) {
            impl.env.push_back(exp);
        }
        
        size_t coeff_offset = impl.env.size();
        impl.bas[BAS_SLOTS * i + PTR_COEFF] = coeff_offset;
        
        // Normalize contraction coefficients for libcint
        // libcint expects normalized primitive GTO coefficients
        // Normalization: c_norm = c_raw * sqrt(gaussian_int(2l+2, 2*alpha))
        // where gaussian_int(n, alpha) = gamma((n+1)/2) / (2 * alpha^((n+1)/2))
        int l = shell.l + shell.m + shell.n;
        for (size_t j = 0; j < shell.exponents.size(); ++j) {
            double alpha = shell.exponents[j];
            double c_raw = shell.contractions[j];
            // gto_norm = 1/sqrt(gaussian_int(2l+2, 2*alpha))
            // gaussian_int = gamma(l+1.5) / (2 * (2*alpha)^(l+1.5))
            double l_half = l + 1.5;
            double gaussian_int = std::exp(std::lgamma(l_half)) / (2.0 * std::pow(2.0 * alpha, l_half));
            double norm = 1.0 / std::sqrt(gaussian_int);
            impl.env.push_back(c_raw * norm);
        }
    }
    
    // Initialize optimizer
    if (impl.opt) {
        CINTdel_optimizer(&impl.opt);
    }
    CINTinit_optimizer(&impl.opt, impl.atm.data(), impl.natm,
                       impl.bas.data(), impl.nbas, impl.env.data());
    
    // Calculate total number of Cartesian GTOs
    impl.nbf = CINTtot_cgto_cart(impl.bas.data(), impl.nbas);
    
    // Calculate AO offsets - manually compute for Cartesian GTOs
    impl.ao_loc.resize(impl.nbas + 1);
    impl.ao_loc[0] = 0;
    for (int i = 0; i < impl.nbas; ++i) {
        int n_cart = CINTcgtos_cart(i, impl.bas.data());
        impl.ao_loc[i + 1] = impl.ao_loc[i] + n_cart;
    }
}

void IntegralEngine::compute_overlap(const Molecule& mol, double* S) {
    auto& impl = *pimpl_;
    std::fill(S, S + impl.nbf * impl.nbf, 0.0);
    
    // Compute overlap for each shell pair
    for (int i = 0; i < impl.nbas; ++i) {
        for (int j = 0; j < impl.nbas; ++j) {
            FINT shls[2] = {i, j};
            
            int ni = CINTcgtos_cart(i, impl.bas.data());
            int nj = CINTcgtos_cart(j, impl.bas.data());
            int ao_i = impl.ao_loc[i];
            int ao_j = impl.ao_loc[j];
            
            std::vector<double> buf(ni * nj);
            cint1e_ovlp_cart(buf.data(), shls, impl.atm.data(), impl.natm,
                            impl.bas.data(), impl.nbas, impl.env.data());
            
            // Copy to output matrix (both ij and ji for symmetry)
            for (int ii = 0; ii < ni; ++ii) {
                for (int jj = 0; jj < nj; ++jj) {
                    S[(ao_i + ii) * impl.nbf + (ao_j + jj)] = buf[ii * nj + jj];
                    S[(ao_j + jj) * impl.nbf + (ao_i + ii)] = buf[ii * nj + jj];
                }
            }
        }
    }
}

void IntegralEngine::compute_kinetic(const Molecule& mol, double* T) {
    auto& impl = *pimpl_;
    std::fill(T, T + impl.nbf * impl.nbf, 0.0);
    
    for (int i = 0; i < impl.nbas; ++i) {
        for (int j = 0; j < impl.nbas; ++j) {
            FINT shls[2] = {i, j};
            
            int ni = CINTcgtos_cart(i, impl.bas.data());
            int nj = CINTcgtos_cart(j, impl.bas.data());
            int ao_i = impl.ao_loc[i];
            int ao_j = impl.ao_loc[j];
            
            std::vector<double> buf(ni * nj);
            cint1e_kin_cart(buf.data(), shls, impl.atm.data(), impl.natm,
                           impl.bas.data(), impl.nbas, impl.env.data());
            
            for (int ii = 0; ii < ni; ++ii) {
                for (int jj = 0; jj < nj; ++jj) {
                    T[(ao_i + ii) * impl.nbf + (ao_j + jj)] = buf[ii * nj + jj];
                    T[(ao_j + jj) * impl.nbf + (ao_i + ii)] = buf[ii * nj + jj];
                }
            }
        }
    }
}

void IntegralEngine::compute_nuclear(const Molecule& mol, double* V) {
    auto& impl = *pimpl_;
    std::fill(V, V + impl.nbf * impl.nbf, 0.0);
    
    for (int i = 0; i < impl.nbas; ++i) {
        for (int j = 0; j < impl.nbas; ++j) {
            FINT shls[2] = {i, j};
            
            int ni = CINTcgtos_cart(i, impl.bas.data());
            int nj = CINTcgtos_cart(j, impl.bas.data());
            int ao_i = impl.ao_loc[i];
            int ao_j = impl.ao_loc[j];
            
            std::vector<double> buf(ni * nj);
            cint1e_nuc_cart(buf.data(), shls, impl.atm.data(), impl.natm,
                           impl.bas.data(), impl.nbas, impl.env.data());
            
            for (int ii = 0; ii < ni; ++ii) {
                for (int jj = 0; jj < nj; ++jj) {
                    V[(ao_i + ii) * impl.nbf + (ao_j + jj)] = buf[ii * nj + jj];
                    V[(ao_j + jj) * impl.nbf + (ao_i + ii)] = buf[ii * nj + jj];
                }
            }
        }
    }
}

double* IntegralEngine::compute_eri(const Molecule& mol, size_t& size) {
    auto& impl = *pimpl_;
    int nbf = impl.nbf;
    
    // Allocate full (ij|kl) ERI buffer
    size_t neri = nbf * nbf * nbf * nbf;
    size = neri;
    double* eri = new double[neri];
    std::fill(eri, eri + neri, 0.0);
    
    // Precompute AO offsets
    std::vector<int> ao_offsets(impl.nbas + 1);
    for (int i = 0; i <= impl.nbas; ++i) {
        ao_offsets[i] = impl.ao_loc[i];
    }
    
    // Compute ERIs shell by shell
    for (int i = 0; i < impl.nbas; ++i) {
        int ni = CINTcgtos_cart(i, impl.bas.data());
        int off_i = ao_offsets[i];
        
        for (int j = 0; j < impl.nbas; ++j) {
            int nj = CINTcgtos_cart(j, impl.bas.data());
            int off_j = ao_offsets[j];
            
            for (int k = 0; k < impl.nbas; ++k) {
                int nk = CINTcgtos_cart(k, impl.bas.data());
                int off_k = ao_offsets[k];
                
                for (int l = 0; l < impl.nbas; ++l) {
                    int nl = CINTcgtos_cart(l, impl.bas.data());
                    int off_l = ao_offsets[l];
                    
                    FINT shls[4] = {i, j, k, l};
                    std::vector<double> buf(ni * nj * nk * nl);
                    
                    FINT not0 = cint2e_cart(buf.data(), shls, impl.atm.data(), impl.natm,
                                          impl.bas.data(), impl.nbas, impl.env.data(), nullptr);
                    
                    if (not0) {
                        // Copy to full ERI buffer with correct indexing
                        // libcint returns (ij|kl) with buf[i,j,k,l] = buf[i*nj*nk*nl + j*nk*nl + k*nl + l]
                        for (int ii = 0; ii < ni; ++ii) {
                            for (int jj = 0; jj < nj; ++jj) {
                                for (int kk = 0; kk < nk; ++kk) {
                                    for (int ll = 0; ll < nl; ++ll) {
                                        int i_ao = off_i + ii;
                                        int j_ao = off_j + jj;
                                        int k_ao = off_k + kk;
                                        int l_ao = off_l + ll;
                                        
                                        // libcint buffer indexing: [i][j][k][l]
                                        size_t buf_idx = ((ii * nj + jj) * nk + kk) * nl + ll;
                                        double val = buf[buf_idx];
                                        
                                        // Store as (ij|kl)
                                        size_t idx = ((i_ao * nbf + j_ao) * nbf + k_ao) * nbf + l_ao;
                                        eri[idx] = val;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Progress indicator
        if ((i + 1) % 5 == 0) {
            std::cerr << ".";
        }
    }
    std::cerr << std::endl;
    
    return eri;
}

void IntegralEngine::eval_ao(const Molecule& mol, const double* grid, int n_points, double* ao_values) {
    // Simplified AO evaluation - for full implementation,
    // use libcint's int1e_grid_ao or similar function
    auto& impl = *pimpl_;
    
    int nbf = impl.nbf;
    std::fill(ao_values, ao_values + nbf * n_points, 0.0);
    
    // For each grid point and each shell, evaluate contracted GTOs
    for (int ipoint = 0; ipoint < n_points; ++ipoint) {
        double x = grid[3 * ipoint];
        double y = grid[3 * ipoint + 1];
        double z = grid[3 * ipoint + 2];
        
        int ao_offset = 0;
        for (int ishell = 0; ishell < impl.nbas; ++ishell) {
            const auto& shell = mol.get_shell(ishell);
            const auto& atom = mol.get_atom(shell.atom_index);
            
            double dx = x - atom.x;
            double dy = y - atom.y;
            double dz = z - atom.z;
            double r2 = dx*dx + dy*dy + dz*dz;
            
            int n_ao = CINTcgtos_cart(ishell, impl.bas.data());
            
            // Evaluate contracted GTO
            double gto_val = 0.0;
            for (size_t i = 0; i < shell.exponents.size(); ++i) {
                gto_val += shell.contractions[i] * std::exp(-shell.exponents[i] * r2);
            }
            
            // Multiply by angular part
            for (int i = 0; i < n_ao; ++i) {
                ao_values[ipoint * nbf + ao_offset + i] = gto_val;
            }
            
            ao_offset += n_ao;
        }
    }
}

} // namespace gto
} // namespace pyscf
