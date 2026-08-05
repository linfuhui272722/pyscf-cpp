#include "gto_integrals.h"
#include "molecule.h"
#include <cmath>
#include <cstring>
#include <iostream>

// For iOS builds without libcint, provide stub implementations
#ifdef __APPLE__
#ifndef HAVE_LIBCINT
#define STUB_IMPLEMENTATION 1
#endif
#endif

#ifdef HAVE_LIBCINT
extern "C" {
#include <cint.h>
}
#endif

namespace pyscf {
namespace gto {

// Stub implementation for iOS builds
#ifdef STUB_IMPLEMENTATION

struct IntegralEngine::Impl {
    int nbf = 0;
};

IntegralEngine::IntegralEngine() : pimpl_(new Impl) {}
IntegralEngine::~IntegralEngine() = default;

void IntegralEngine::initialize(const Molecule& mol) {
    pimpl_->nbf = mol.num_shells();  // Rough approximation
    std::cerr << "Warning: Using stub integral engine (libcint not available)" << std::endl;
}

void IntegralEngine::compute_overlap(const Molecule& mol, double* S) {
    int n = pimpl_->nbf;
    std::fill(S, S + n * n, 0.0);
    // Identity matrix as stub
    for (int i = 0; i < n; ++i) S[i * n + i] = 1.0;
}

void IntegralEngine::compute_kinetic(const Molecule& mol, double* T) {
    int n = pimpl_->nbf;
    std::fill(T, T + n * n, 0.0);
}

void IntegralEngine::compute_nuclear(const Molecule& mol, double* V) {
    int n = pimpl_->nbf;
    std::fill(V, V + n * n, 0.0);
}

double* IntegralEngine::compute_eri(const Molecule& mol, size_t& size) {
    int n = pimpl_->nbf;
    size = n * n * n * n;
    double* eri = new double[size];
    std::fill(eri, eri + size, 0.0);
    return eri;
}

void IntegralEngine::eval_ao(const Molecule& mol, const double* grid, int n_points, double* ao_values) {
    int nbf = pimpl_->nbf;
    std::fill(ao_values, ao_values + nbf * n_points, 0.0);
}

#else  // HAVE_LIBCINT

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
double normalize_gto(double alpha, int l, int m, int n, double c) {
    int lmn = l + m + n;
    const double PI = 3.14159265358979323846;
    
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

// Internal implementation using libcint
struct IntegralEngine::Impl {
    std::vector<FINT> atm;
    std::vector<FINT> bas;
    std::vector<double> env;
    CINTOpt* opt = nullptr;
    int natm = 0;
    int nbas = 0;
    int nbf = 0;
    std::vector<FINT> ao_loc;
    
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
    
    impl.atm.clear();
    impl.bas.clear();
    impl.env.clear();
    impl.natm = mol.num_atoms();
    impl.nbas = mol.num_shells();
    
    impl.env.resize(20);
    impl.env[PTR_EXPCUTOFF] = -100.0;
    
    size_t env_offset = impl.env.size();
    
    impl.atm.resize(ATM_SLOTS * impl.natm);
    for (int i = 0; i < impl.natm; ++i) {
        const auto& atom = mol.get_atom(i);
        impl.atm[ATM_SLOTS * i + CHARGE_OF] = atom.atomic_number;
        impl.atm[ATM_SLOTS * i + PTR_COORD] = env_offset;
        impl.atm[ATM_SLOTS * i + NUC_MOD_OF] = POINT_NUC;
        impl.atm[ATM_SLOTS * i + 3] = 0;
        impl.atm[ATM_SLOTS * i + 4] = 0;
        impl.atm[ATM_SLOTS * i + 5] = 0;
        
        impl.env.push_back(atom.x);
        impl.env.push_back(atom.y);
        impl.env.push_back(atom.z);
        env_offset = impl.env.size();
    }
    
    impl.bas.resize(BAS_SLOTS * impl.nbas);
    for (int i = 0; i < impl.nbas; ++i) {
        const auto& shell = mol.get_shell(i);
        impl.bas[BAS_SLOTS * i + ATOM_OF] = shell.atom_index;
        impl.bas[BAS_SLOTS * i + ANG_OF] = shell.l + shell.m + shell.n;
        impl.bas[BAS_SLOTS * i + NPRIM_OF] = shell.exponents.size();
        impl.bas[BAS_SLOTS * i + NCTR_OF] = 1;
        impl.bas[BAS_SLOTS * i + KAPPA_OF] = 0;
        
        size_t exp_offset = impl.env.size();
        impl.bas[BAS_SLOTS * i + PTR_EXP] = exp_offset;
        
        for (double exp : shell.exponents) {
            impl.env.push_back(exp);
        }
        
        size_t coeff_offset = impl.env.size();
        impl.bas[BAS_SLOTS * i + PTR_COEFF] = coeff_offset;
        
        int l = shell.l + shell.m + shell.n;
        for (size_t j = 0; j < shell.exponents.size(); ++j) {
            double alpha = shell.exponents[j];
            double c_raw = shell.contractions[j];
            double l_half = l + 1.5;
            double gaussian_int = std::exp(std::lgamma(l_half)) / (2.0 * std::pow(2.0 * alpha, l_half));
            double norm = 1.0 / std::sqrt(gaussian_int);
            impl.env.push_back(c_raw * norm);
        }
    }
    
    if (impl.opt) {
        CINTdel_optimizer(&impl.opt);
    }
    CINTinit_optimizer(&impl.opt, impl.atm.data(), impl.natm,
                       impl.bas.data(), impl.nbas, impl.env.data());
    
    impl.nbf = CINTtot_cgto_cart(impl.bas.data(), impl.nbas);
    
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
    
    size_t neri = nbf * nbf * nbf * nbf;
    size = neri;
    double* eri = new double[neri];
    std::fill(eri, eri + neri, 0.0);
    
    std::vector<int> ao_offsets(impl.nbas + 1);
    for (int i = 0; i <= impl.nbas; ++i) {
        ao_offsets[i] = impl.ao_loc[i];
    }
    
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
                        for (int ii = 0; ii < ni; ++ii) {
                            for (int jj = 0; jj < nj; ++jj) {
                                for (int kk = 0; kk < nk; ++kk) {
                                    for (int ll = 0; ll < nl; ++ll) {
                                        int i_ao = off_i + ii;
                                        int j_ao = off_j + jj;
                                        int k_ao = off_k + kk;
                                        int l_ao = off_l + ll;
                                        
                                        size_t buf_idx = ((ii * nj + jj) * nk + kk) * nl + ll;
                                        double val = buf[buf_idx];
                                        
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
        
        if ((i + 1) % 5 == 0) {
            std::cerr << ".";
        }
    }
    std::cerr << std::endl;
    
    return eri;
}

void IntegralEngine::eval_ao(const Molecule& mol, const double* grid, int n_points, double* ao_values) {
    auto& impl = *pimpl_;
    
    int nbf = impl.nbf;
    std::fill(ao_values, ao_values + nbf * n_points, 0.0);
    
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
            
            double gto_val = 0.0;
            for (size_t i = 0; i < shell.exponents.size(); ++i) {
                gto_val += shell.contractions[i] * std::exp(-shell.exponents[i] * r2);
            }
            
            for (int i = 0; i < n_ao; ++i) {
                ao_values[ipoint * nbf + ao_offset + i] = gto_val;
            }
            
            ao_offset += n_ao;
        }
    }
}

#endif  // HAVE_LIBCINT

} // namespace gto
} // namespace pyscf
