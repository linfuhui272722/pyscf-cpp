#ifndef PYSCF_CPP_OPTIMIZER_H
#define PYSCF_CPP_OPTIMIZER_H

#include "scf.h"
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>

namespace pyscf {
namespace geom {

// Optimization convergence criteria (same as geomeTRIC defaults)
struct OptConvergence {
    double gradient_max = 4.5e-4;    // Eh/Bohr
    double gradient_rms = 3.0e-4;    // Eh/Bohr  
    double step_max = 1.8e-3;       // Bohr
    double step_rms = 1.2e-3;       // Bohr
    double energy_tol = 1.0e-6;     // Eh
    int max_iterations = 100;
    bool verbose = true;
};

// Optimization result
struct OptimizationResult {
    std::vector<double> coordinates;  // Final coordinates in Bohr
    double final_energy;
    std::vector<double> final_gradient;
    int iterations;
    bool converged;
    std::vector<double> energy_history;
    std::vector<double> gradient_norms;
    std::string message;
};

// BFGS optimizer for geometry optimization
// Reference: Nocedal & Wright, Numerical Optimization, Chapter 6
class BFGSOptimizer {
public:
    BFGSOptimizer() : conv_(), verbose_(true) {}
    ~BFGSOptimizer() = default;

    void set_convergence(const OptConvergence& conv) { conv_ = conv; }
    void set_verbose(bool v) { verbose_ = v; }
    const OptConvergence& convergence() const { return conv_; }

    // Main optimization interface
    OptimizationResult optimize(
        std::shared_ptr<dft::RKS> scf,
        std::shared_ptr<Molecule> mol_init
    );

private:
    OptConvergence conv_;
    bool verbose_;

    // Hessian approximation (inverse Hessian for BFGS)
    std::vector<double> inv_hessian_;
    int n_variables_;

    // Initialize Hessian (identity matrix for first iteration)
    void initialize_hessian(int n);
    
    // BFGS Hessian update formula:
    // B_new = B + (y*y^T)/(y^T*s) - (B*s*s^T*B)/(s^T*B*s)
    // where s = dx, y = dg
    void update_inverse_hessian(
        const std::vector<double>& dx,
        const std::vector<double>& dg
    );
    
    // Compute search direction: p = -H * g
    std::vector<double> compute_search_direction(const std::vector<double>& gradient);
    
    // Strong Wolfe condition line search
    // Returns step length alpha
    double line_search(
        std::shared_ptr<dft::RKS> scf,
        std::shared_ptr<Molecule> mol,
        const std::vector<double>& position,
        const std::vector<double>& direction,
        double energy_init,
        const std::vector<double>& gradient_init
    );
    
    // Check convergence criteria
    bool check_convergence(
        const std::vector<double>& gradient,
        const std::vector<double>& step,
        double energy,
        double energy_prev
    );
    
    // Helper: vector operations
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    static double norm(const std::vector<double>& v);
    static void axpy(double a, const std::vector<double>& x, std::vector<double>& y);
    
    // Matrix-vector product: y = H * x
    void hessian_times_vector(const std::vector<double>& x, std::vector<double>& y);
};

// Scanner to evaluate energy and gradient at given coordinates
class Scanner {
public:
    Scanner(std::shared_ptr<dft::RKS> scf) : scf_(scf), mol_(nullptr) {}
    
    void set_molecule(std::shared_ptr<Molecule> mol) { mol_ = mol; }
    
    // Evaluate energy and gradient at given coordinates
    // Returns: {energy, gradient}
    std::pair<double, std::vector<double>> operator()(const std::vector<double>& coords);
    
private:
    std::shared_ptr<dft::RKS> scf_;
    std::shared_ptr<Molecule> mol_;
    
    // Numerical gradient computation (central difference)
    std::vector<double> compute_numerical_gradient(
        std::shared_ptr<dft::RKS> scf,
        std::shared_ptr<Molecule> mol,
        const std::vector<double>& coords,
        double h = 1e-5
    );
};

// Compute nuclear gradient from electronic structure
// Reference: PySCF grad/rhf.py
std::vector<double> compute_nuclear_gradient(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol
);

} // namespace geom
} // namespace pyscf

#endif // PYSCF_CPP_OPTIMIZER_H
