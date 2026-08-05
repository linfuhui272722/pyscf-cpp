#include "b3lyp_xc.h"
#include <cmath>

namespace {
constexpr double PI = 3.14159265358979323846;
}

namespace pyscf {
namespace dft {

// B3LYP parameters
constexpr double B3LYP_A = 0.05;
constexpr double B3LYP_B = 0.05;
constexpr double B3LYP_C = 0.8;
constexpr double B3LYP_GAMMA = 0.004;

// VWN5 LDA correlation parameters
constexpr double VWN_A = 0.0310907;
constexpr double VWN_B = 3.72744;
constexpr double VWN_C = 7.06079;
constexpr double VWN_X0 = -1.06840;

// Becke88 beta parameter
constexpr double BECKE88_BETA = 0.0042;

// LYP parameters
constexpr double LYP_C1 = 3.0;
constexpr double LYP_C2 = 1.0/3.0;
constexpr double LYP_C3 = 14.0;
constexpr double LYP_C4 = 1.0/3.0;
constexpr double LYP_C5 = 1.0/3.0;
constexpr double LYP_C6 = 0.00330;
constexpr double LYP_C7 = 1.0/3.0;
constexpr double LYP_C8 = 2.0/3.0;
constexpr double LYP_C9 = 0.001667;

// LDA exchange (Slater-Dirac)
double slater_exchange(double rho) {
    if (rho <= 0.0) return 0.0;
    return -std::pow(3.0 * rho / PI, 1.0/3.0);
}

// LDA correlation (VWN formula V)
double vwn_polarization_function(double zeta) {
    double r = 1.0 + zeta;
    double u = 1.0 - zeta;
    return 0.5 * (std::pow(r, 4.0/3.0) + std::pow(u, 4.0/3.0)) - 
           2.0 * (std::pow(r, 2.0/3.0) + std::pow(u, 2.0/3.0)) + 
           (r + u) - 2.0;
}

double lda_c_vwn(double rho, double zeta, double beta) {
    if (rho <= 0.0) return 0.0;
    
    double r = std::pow(rho, 1.0/3.0);
    double X = r * r + beta * r + VWN_B;
    double Q = std::sqrt(4.0 * VWN_C - VWN_B * VWN_B);
    double X0 = VWN_X0;
    
    double vwn_param = (std::log(r * r / X) - VWN_B * X0 / (X0 * X0 + VWN_B * X0 + VWN_C) + 
                        2.0 * VWN_B / Q * std::atan(Q / (2.0 * r + VWN_B))) / (2.0 * PI);
    
    // Unpolarized
    if (std::abs(zeta) < 1e-10) {
        return VWN_A * vwn_param;
    }
    
    // Polarized
    double alpha = -VWN_A * std::log(1.0 + 1.0/VWN_A);
    double z4 = zeta * zeta * zeta * zeta;
    
    return vwn_param * (1.0 - z4) + 
           alpha * vwn_polarization_function(zeta) * z4;
}

// LDA correlation energy density (epsilon_c)
double lda_c_vwn_energy(double rho, double zeta) {
    if (rho <= 0.0) return 0.0;
    return lda_c_vwn(rho, zeta, VWN_C);
}

// Becke88 gradient correction to exchange
double becke88_exchange(double rho, double gamma, double beta) {
    if (rho <= 0.0 || gamma <= 0.0) return 0.0;
    
    double x = std::sqrt(gamma) / std::pow(rho, 4.0/3.0);
    double g = beta * x * x / (1.0 + 6.0 * beta * x * std::atanh(std::sqrt(gamma) * x / rho));
    
    return -beta * gamma * g / rho;
}

// LYP correlation energy density
double lyp_correlation(double rhoa, double rhob, 
                       double gamma_aa, double gamma_ab, double gamma_bb) {
    if (rhoa <= 0.0 || rhob <= 0.0) return 0.0;
    
    double rho = rhoa + rhob;
    if (rho <= 0.0) return 0.0;
    
    double omega_ab = std::sqrt(rhoa * rhob);
    double C1 = LYP_C1 * rhoa * rhob / rho;
    double C2 = LYP_C2 * (gamma_aa * rhoa * rhoa + 2.0 * gamma_ab * rhoa * rhob + 
                         gamma_bb * rhob * rhob);
    double C3 = -LYP_C3 * omega_ab * (std::pow(rhoa, LYP_C4) + std::pow(rhob, LYP_C4));
    double C4 = LYP_C5 * rho * rho * rho * rho;
    double C5 = LYP_C6 * rho * rho * rho;
    double C6 = LYP_C7 * rho * rho;
    double C7 = LYP_C8 * rho;
    double C8 = LYP_C9 * omega_ab * (rhoa * rhoa + rhob * rhob);
    
    double d = C4 + C5 * (rhoa * rhoa + rhob * rhob) + C6 * rho * omega_ab;
    double gamma_sum = C2 + C3;
    
    double term1 = -C1 / rho - gamma_sum * C7 / 8.0;
    double term2 = (C2 * rho + C3) / (rho * rho * rho * rho * rho);
    double term3 = -LYP_C8 / 8.0 * gamma_sum / rho;
    
    double epsilon = term1 + term2 * d + term3;
    
    // More accurate LYP formula
    double omega = rhoa * rhob;
    double omega2 = omega * omega;
    double rho_pow = std::pow(rho, -1.0/3.0);
    double rho_pow2 = rho_pow * rho_pow;
    double rho_43 = std::pow(rho, 4.0/3.0);
    
    double delta = gamma_aa * rhoa * rhoa + 2.0 * gamma_ab * omega + gamma_bb * rhob * rhob;
    
    double a = 1.0 / (2.0 * std::pow(PI, 2.0/3.0));
    double b = LYP_C3;
    double c = LYP_C6;
    double d_param = LYP_C8;
    
    double term_a = -omega / rho * (a * rho_pow - c / 8.0 * delta / rho_43 - d_param / 8.0 * gamma_sum / rho);
    double term_b = -a / 18.0 * rhoa * rhob * std::pow(rhoa * rhob, -2.0/3.0);
    double term_c = c / 18.0 * omega * delta / rho_pow2;
    double term_d = d_param / 72.0 * omega * gamma_sum;
    
    epsilon = term_a + b * term_b + c * term_c + d_param * term_d;
    
    return epsilon * rho;
}

B3LYPFunctional::B3LYPFunctional() = default;

XCResult evaluate_b3lyp_xc(double rho, double gamma) {
    XCResult result = {0.0, 0.0, 0.0};
    
    if (rho <= 1e-10) return result;
    
    // Slater exchange (80%)
    double exc_slater = 0.8 * slater_exchange(rho);
    
    // Becke88 correction (80% of 0.42)
    double exc_becke88 = 0.8 * becke88_exchange(rho, gamma, B3LYP_GAMMA);
    
    // VWN correlation (19%)
    double exc_vwn = 0.19 * lda_c_vwn_energy(rho, 0.0);
    
    result.exc = exc_slater + exc_becke88 + exc_vwn;
    
    // Simplified potential (for GGA)
    result.vxc = 0.0;
    result.vsigma = 0.8 * B3LYP_GAMMA;  // Approximation
    
    return result;
}

GGAResult evaluate_b3lyp_gga(double rhoa, double rhob,
                               double gamma_aa, double gamma_ab, double gamma_bb) {
    GGAResult result = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    if (rhoa <= 1e-10 || rhob <= 1e-10) return result;
    
    double rho = rhoa + rhob;
    
    // Exchange part (Slater + Becke88)
    double exc_x = 0.8 * slater_exchange(rho);
    double exc_becke = 0.8 * becke88_exchange(rho, gamma_aa + 2.0 * gamma_ab + gamma_bb, B3LYP_GAMMA);
    
    // Correlation part (LYP + VWN)
    double exc_lyp = 0.81 * lyp_correlation(rhoa, rhob, gamma_aa, gamma_ab, gamma_bb);
    double exc_vwn = 0.19 * rho * lda_c_vwn_energy(rho, 0.0);
    
    result.exc = exc_x + exc_becke + exc_lyp + exc_vwn;
    
    // Potential derivatives (simplified)
    result.vrho_a = 0.8 * (-4.0/3.0) * std::pow(3.0 * rho / PI, 1.0/3.0);
    result.vrho_b = result.vrho_a;
    
    // sigma derivatives
    result.vsigma_aa = 0.8 * B3LYP_GAMMA;
    result.vsigma_ab = 1.6 * B3LYP_GAMMA;
    result.vsigma_bb = 0.8 * B3LYP_GAMMA;
    
    return result;
}

void B3LYPFunctional::compute_xc_energy(const std::vector<double>& rho,
                                       const std::vector<double>& gamma,
                                       double& exc, std::vector<double>& vxc) {
    exc = 0.0;
    vxc.resize(rho.size(), 0.0);
    
    for (size_t i = 0; i < rho.size(); ++i) {
        auto result = evaluate_b3lyp_xc(rho[i], gamma[i]);
        exc += result.exc;
        vxc[i] = result.vxc;
    }
}

void B3LYPFunctional::compute_gga_potential(const std::vector<double>& rhoa,
                                            const std::vector<double>& rhob,
                                            const std::vector<double>& gamma_aa,
                                            const std::vector<double>& gamma_ab,
                                            const std::vector<double>& gamma_bb,
                                            std::vector<double>& vrhoa,
                                            std::vector<double>& vrhob,
                                            std::vector<double>& vsigma_aa,
                                            std::vector<double>& vsigma_ab,
                                            std::vector<double>& vsigma_bb) {
    size_t n = rhoa.size();
    vrhoa.resize(n);
    vrhob.resize(n);
    vsigma_aa.resize(n);
    vsigma_ab.resize(n);
    vsigma_bb.resize(n);
    
    for (size_t i = 0; i < n; ++i) {
        auto result = evaluate_b3lyp_gga(rhoa[i], rhob[i], 
                                         gamma_aa[i], gamma_ab[i], gamma_bb[i]);
        vrhoa[i] = result.vrho_a;
        vrhob[i] = result.vrho_b;
        vsigma_aa[i] = result.vsigma_aa;
        vsigma_ab[i] = result.vsigma_ab;
        vsigma_bb[i] = result.vsigma_bb;
    }
}

void compute_density(const std::vector<double>& ao_values,
                    const std::vector<double>& mo_coefficients,
                    int n_occupied,
                    std::vector<double>& rho) {
    // Simplified: rho = sum_{occ} |phi_i|^2
    // where phi_i = sum_mu C_mu_i * chi_mu
    int n_ao = mo_coefficients.size() / n_occupied;
    int n_grid = ao_values.size() / n_ao;
    
    rho.assign(n_grid, 0.0);
    
    // For each occupied orbital
    for (int i = 0; i < n_occupied; ++i) {
        // For each grid point
        for (int j = 0; j < n_grid; ++j) {
            double rho_contrib = 0.0;
            for (int mu = 0; mu < n_ao; ++mu) {
                rho_contrib += mo_coefficients[i * n_ao + mu] * 
                              ao_values[j * n_ao + mu];
            }
            rho[j] += 2.0 * rho_contrib * rho_contrib;  // Factor 2 for spin
        }
    }
}

void compute_density_gradient(const std::vector<double>& ao_values,
                              const std::vector<double>& ao_grad_x,
                              const std::vector<double>& ao_grad_y,
                              const std::vector<double>& ao_grad_z,
                              const std::vector<double>& mo_coefficients,
                              int n_occupied,
                              std::vector<double>& gamma) {
    int n_ao = mo_coefficients.size() / n_occupied;
    int n_grid = ao_values.size() / n_ao;
    
    gamma.assign(n_grid, 0.0);
    
    for (int i = 0; i < n_occupied; ++i) {
        for (int j = 0; j < n_grid; ++j) {
            double gx = 0.0, gy = 0.0, gz = 0.0;
            for (int mu = 0; mu < n_ao; ++mu) {
                double c = mo_coefficients[i * n_ao + mu];
                gx += c * ao_grad_x[j * n_ao + mu];
                gy += c * ao_grad_y[j * n_ao + mu];
                gz += c * ao_grad_z[j * n_ao + mu];
            }
            gamma[j] += 4.0 * (gx*gx + gy*gy + gz*gz);
        }
    }
}

} // namespace dft
} // namespace pyscf
