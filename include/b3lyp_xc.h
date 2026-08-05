#ifndef PYSCF_CPP_B3LYP_XC_H
#define PYSCF_CPP_B3LYP_XC_H

#include <vector>
#include <cmath>

namespace pyscf {
namespace dft {

// B3LYP functional coefficients
struct B3LYPCoefficients {
    // Exchange: 0.8*B88 + 0.2*HF
    double alpha_hf = 0.2;
    double alpha_slater = 0.8;
    double beta_becke88 = 0.42;
    
    // Correlation: 0.81*LYP + 0.19*VWN
    double alpha_lyp = 0.81;
    double alpha_vwn = 0.19;
    
    B3LYPCoefficients() = default;
};

// LDA correlation (VWN formula V)
double lda_c_vwn(double rho, double beta = 0.003);

// LDA correlation potential (derivative of VWN)
double lda_c_vwn_polarization(double zeta);

// Slater-Dirac exchange energy density
double slater_exchange(double rho);

// Becke88 gradient correction to exchange
double becke88_exchange(double rho, double gamma, double beta);

// LYP correlation energy density
double lyp_correlation(double rhoa, double rhob, 
                       double gamma_aa, double gamma_ab, double gamma_bb);

// Compute XC energy and potential for B3LYP
// Input: electron density rho, gradient gamma = |nabla rho|^2
// Output: exc (energy per unit volume), vxc (potential d(rho*exc)/d(rho))
struct XCResult {
    double exc;       // Exchange-correlation energy density
    double vxc;       // XC potential (zero for GGA)
    double vsigma;    // Derivative w.r.t. gamma (for GGA potential)
};

XCResult evaluate_b3lyp_xc(double rho, double gamma);

// GGA evaluation with both alpha and beta densities
struct GGAResult {
    double exc;        // Exchange-correlation energy density
    double vrho_a;    // Potential w.r.t. alpha density
    double vrho_b;    // Potential w.r.t. beta density
    double vsigma_aa;  // Potential w.r.t. gamma_aa
    double vsigma_ab;  // Potential w.r.t. gamma_ab
    double vsigma_bb;  // Potential w.r.t. gamma_bb
};

GGAResult evaluate_b3lyp_gga(double rhoa, double rhob,
                               double gamma_aa, double gamma_ab, double gamma_bb);

// B3LYP XC functional class
class B3LYPFunctional {
public:
    B3LYPFunctional();
    
    // Evaluate XC contribution on a grid
    void compute_xc_energy(const std::vector<double>& rho,
                          const std::vector<double>& gamma,
                          double& exc, std::vector<double>& vxc);
    
    // Compute GGA potential
    void compute_gga_potential(const std::vector<double>& rhoa,
                               const std::vector<double>& rhob,
                               const std::vector<double>& gamma_aa,
                               const std::vector<double>& gamma_ab,
                               const std::vector<double>& gamma_bb,
                               std::vector<double>& vrhoa,
                               std::vector<double>& vrhob,
                               std::vector<double>& vsigma_aa,
                               std::vector<double>& vsigma_ab,
                               std::vector<double>& vsigma_bb);
    
    // HF exchange contribution (20% in B3LYP)
    double get_hf_exchange_coefficient() const { return coef_.alpha_hf; }
    
private:
    B3LYPCoefficients coef_;
};

// Compute electron density and its gradients from AO basis
// rho(r) = sum_i |phi_i(r)|^2 (for restricted, multiply by 2 for occupied)
void compute_density(const std::vector<double>& ao_values,
                     const std::vector<double>& mo_coefficients,
                     int n_occupied,
                     std::vector<double>& rho);

// Compute gradient of density
void compute_density_gradient(const std::vector<double>& ao_values,
                              const std::vector<double>& ao_grad_x,
                              const std::vector<double>& ao_grad_y,
                              const std::vector<double>& ao_grad_z,
                              const std::vector<double>& mo_coefficients,
                              int n_occupied,
                              std::vector<double>& gamma);

} // namespace dft
} // namespace pyscf

#endif // PYSCF_CPP_B3LYP_XC_H
