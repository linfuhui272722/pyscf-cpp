/**
 * Geometry Optimizer Implementation
 * 
 * Implements BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm for
 * molecular geometry optimization.
 * 
 * Reference:
 * - Nocedal & Wright, Numerical Optimization, Chapter 6
 * - PySCF geomopt/geometric_solver.py
 * - PySCF geomopt/berny_solver.py
 */

#include "optimizer.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

extern "C" {
// BLAS/LAPACK interfaces
void dgemm_(char* transa, char* transb, int* m, int* n, int* k,
           double* alpha, double* a, int* lda, double* b, int* ldb,
           double* beta, double* c, int* ldc);
}

namespace pyscf {
namespace geom {

// ============================================================================
// BFGS Optimizer Implementation
// ============================================================================

void BFGSOptimizer::initialize_hessian(int n) {
    n_variables_ = n;
    inv_hessian_.assign(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        inv_hessian_[i * n + i] = 1.0;  // Start with identity matrix
    }
}

void BFGSOptimizer::update_inverse_hessian(
    const std::vector<double>& dx,
    const std::vector<double>& dg
) {
    int n = n_variables_;
    
    // Compute s^T * y (scalar)
    double s_dot_y = dot(dx, dg);
    if (std::abs(s_dot_y) < 1e-12) {
        // Skip update if curvature condition not satisfied
        if (verbose_) {
            std::cerr << "Warning: Skipping BFGS update (s^T*y too small)\n";
        }
        return;
    }
    
    // Temporary vectors for matrix operations
    std::vector<double> temp(n);
    std::vector<double> h_dg(n);
    
    // h_dg = H * dg (matrix-vector product)
    hessian_times_vector(dg, h_dg);
    
    // dg^T * H * dg (scalar)
    double dg_h_dg = dot(dg, h_dg);
    
    // BFGS update formula for inverse Hessian:
    // H_new = H + (1 + dg^T*H*dg / s^T*dg) * s*s^T / s^T*dg - s*dg^T*H / s^T*dg - H*dg*s^T / s^T*dg
    double coeff1 = 1.0 + dg_h_dg / s_dot_y;
    double coeff2 = 1.0 / s_dot_y;
    
    // Update inv_hessian_ in-place
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inv_hessian_[i * n + j] += coeff1 * dx[i] * dx[j] * coeff2
                                      - (dx[i] * h_dg[j] + h_dg[i] * dx[j]) / s_dot_y;
        }
    }
}

void BFGSOptimizer::hessian_times_vector(
    const std::vector<double>& x,
    std::vector<double>& y
) {
    int n = n_variables_;
    int m = 1;  // result is n x 1
    int k = n;  // H is n x n
    char trans = 'N';
    double alpha = 1.0, beta = 0.0;
    
    // y = H * x (using inv_hessian_ as approximation)
    // y (n x 1) = H (n x n) * x (n x 1)
    // Note: dgemm requires non-const pointers, so we copy the input
    std::vector<double> x_copy = x;
    dgemm_(&trans, &trans, &n, &m, &k,
           &alpha, inv_hessian_.data(), &n, x_copy.data(), &k,
           &beta, y.data(), &n);
}

std::vector<double> BFGSOptimizer::compute_search_direction(
    const std::vector<double>& gradient
) {
    int n = n_variables_;
    std::vector<double> direction(n);
    
    // direction = -H * gradient (using inv_hessian_)
    hessian_times_vector(gradient, direction);
    for (int i = 0; i < n; ++i) {
        direction[i] = -direction[i];
    }
    
    return direction;
}

double BFGSOptimizer::line_search(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol,
    const std::vector<double>& position,
    const std::vector<double>& direction,
    double energy_init,
    const std::vector<double>& gradient_init
) {
    // Parameters for strong Wolfe conditions
    double c1 = 1e-4;   // Sufficient decrease
    double c2 = 0.9;    // Curvature condition (high for Newton-like methods)
    double alpha_max = 1.0;
    double alpha = 1.0;
    double alpha_lo = 0.0, alpha_hi = 0.0;
    
    // Initial directional derivative
    double dir_deriv_init = dot(gradient_init, direction);
    if (dir_deriv_init >= 0) {
        // Not a descent direction, just return unit step
        return 1.0;
    }
    
    // Create temporary molecule for testing
    auto mol_test = std::make_shared<Molecule>(*mol);
    std::vector<double> coords_test(3 * mol->num_atoms());
    std::vector<double> coords_new(3 * mol->num_atoms());
    
    // Evaluate at alpha = 1.0 first
    for (size_t i = 0; i < position.size(); ++i) {
        coords_new[i] = position[i] + alpha * direction[i];
    }
    mol_test->set_coordinates(coords_new);
    scf->compute();
    double energy_hi = scf->get_total_energy();
    double dir_deriv_hi = 0.0;  // Would need to compute gradient
    
    // Simple backtracking line search
    int max_iter = 20;
    for (int i = 0; i < max_iter; ++i) {
        // Check Armijo condition: f(x + alpha*p) <= f(x) + c1 * alpha * f'(x)^T * p
        if (energy_hi <= energy_init + c1 * alpha * dir_deriv_init) {
            // Sufficient decrease satisfied
            return alpha;
        }
        
        // Zoom to find acceptable step
        if (alpha < 0.5) {
            alpha *= 2.0;
        } else {
            alpha *= 0.5;
        }
        
        for (size_t j = 0; j < position.size(); ++j) {
            coords_new[j] = position[j] + alpha * direction[j];
        }
        mol_test->set_coordinates(coords_new);
        scf->compute();
        energy_hi = scf->get_total_energy();
    }
    
    // If line search fails, return small step
    return alpha;
}

bool BFGSOptimizer::check_convergence(
    const std::vector<double>& gradient,
    const std::vector<double>& step,
    double energy,
    double energy_prev
) {
    // Gradient norms
    double grad_max = 0.0, grad_rms = 0.0;
    for (size_t i = 0; i < gradient.size(); ++i) {
        double abs_g = std::abs(gradient[i]);
        grad_max = std::max(grad_max, abs_g);
        grad_rms += gradient[i] * gradient[i];
    }
    grad_rms = std::sqrt(grad_rms / gradient.size());
    
    // Step norms
    double step_max = 0.0, step_rms = 0.0;
    for (size_t i = 0; i < step.size(); ++i) {
        double abs_s = std::abs(step[i]);
        step_max = std::max(step_max, abs_s);
        step_rms += step[i] * step[i];
    }
    step_rms = std::sqrt(step_rms / step.size());
    
    // Energy change
    double dE = std::abs(energy - energy_prev);
    
    if (verbose_) {
        std::cout << "  Convergence check:\n";
        std::cout << "    |g|_max = " << grad_max << " (tol = " << conv_.gradient_max << ")\n";
        std::cout << "    |g|_rms = " << grad_rms << " (tol = " << conv_.gradient_rms << ")\n";
        std::cout << "    |s|_max = " << step_max << " (tol = " << conv_.step_max << ")\n";
        std::cout << "    |s|_rms = " << step_rms << " (tol = " << conv_.step_rms << ")\n";
        std::cout << "    dE      = " << dE << " (tol = " << conv_.energy_tol << ")\n";
    }
    
    // Check all criteria
    return (grad_max < conv_.gradient_max &&
            grad_rms < conv_.gradient_rms &&
            step_max < conv_.step_max &&
            step_rms < conv_.step_rms &&
            dE < conv_.energy_tol);
}

double BFGSOptimizer::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double BFGSOptimizer::norm(const std::vector<double>& v) {
    return std::sqrt(dot(v, v));
}

void BFGSOptimizer::axpy(double a, const std::vector<double>& x, std::vector<double>& y) {
    for (size_t i = 0; i < x.size(); ++i) {
        y[i] += a * x[i];
    }
}

OptimizationResult BFGSOptimizer::optimize(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol_init
) {
    OptimizationResult result;
    
    // Use the same molecule object, don't create a copy
    auto mol = mol_init;
    int n_atoms = mol->num_atoms();
    n_variables_ = 3 * n_atoms;
    
    if (verbose_) {
        std::cout << "\n=== Geometry Optimization with BFGS ===\n\n";
        std::cout << "Initial geometry:\n";
        for (int i = 0; i < n_atoms; ++i) {
            auto atom = mol->get_atom(i);
            std::cout << "  " << element_symbol(atom.atomic_number);
            std::cout << "  " << std::fixed << std::setprecision(6);
            std::cout << "  " << atom.x << "  " << atom.y << "  " << atom.z << "\n";
        }
        std::cout << "\n";
    }
    
    // Initialize Hessian to identity
    initialize_hessian(n_variables_);
    
    // Get initial coordinates
    auto coords = mol->get_coordinates();
    std::vector<double> gradient_prev(n_variables_, 0.0);
    std::vector<double> coords_prev = coords;
    
    // Scanner for energy and gradient evaluation
    Scanner scanner(scf);
    scanner.set_molecule(mol);
    
    // Initial evaluation
    auto scan_result = scanner(coords);
    double energy_prev = scan_result.first;
    std::vector<double> gradient = scan_result.second;
    
    result.energy_history.push_back(energy_prev);
    result.gradient_norms.push_back(norm(gradient));
    
    if (verbose_) {
        std::cout << "Cycle 0: E = " << std::scientific << std::setprecision(8) << energy_prev;
        std::cout << "  |g| = " << norm(gradient) << "\n";
    }
    
    // Main optimization loop
    for (int iteration = 1; iteration <= conv_.max_iterations; ++iteration) {
        // Compute search direction: p = -H * g
        std::vector<double> search_dir = compute_search_direction(gradient);
        
        // Line search
        double step_length = line_search(scf, mol, coords, search_dir, energy_prev, gradient);
        
        // Update coordinates
        std::vector<double> step(n_variables_);
        for (int i = 0; i < n_variables_; ++i) {
            step[i] = step_length * search_dir[i];
            coords[i] = coords_prev[i] + step[i];
        }
        mol->set_coordinates(coords);
        
        // Evaluate at new point
        auto scan_result_new = scanner(coords);
        double energy = scan_result_new.first;
        std::vector<double> gradient_new = scan_result_new.second;
        
        result.energy_history.push_back(energy);
        result.gradient_norms.push_back(norm(gradient_new));
        
        if (verbose_) {
            std::cout << "Cycle " << iteration << ": E = " << std::scientific << std::setprecision(8) << energy;
            std::cout << "  dE = " << (energy - energy_prev);
            std::cout << "  |g| = " << norm(gradient_new);
            std::cout << "  step = " << norm(step) << "\n";
        }
        
        // Check convergence
        if (check_convergence(gradient_new, step, energy, energy_prev)) {
            result.converged = true;
            result.iterations = iteration;
            result.final_energy = energy;
            result.final_gradient = gradient_new;
            result.coordinates = coords;
            result.message = "Optimization converged successfully";
            
            if (verbose_) {
                std::cout << "\n*** Optimization converged in " << iteration << " cycles ***\n\n";
                std::cout << "Final geometry:\n";
                for (int i = 0; i < n_atoms; ++i) {
                    auto atom = mol->get_atom(i);
                    std::cout << "  " << element_symbol(atom.atomic_number);
                    std::cout << "  " << std::fixed << std::setprecision(6);
                    std::cout << "  " << atom.x << "  " << atom.y << "  " << atom.z << "\n";
                }
            }
            return result;
        }
        
        // BFGS update: only if step is reasonable
        double step_norm = norm(step);
        if (step_norm > 1e-6 && step_norm < 10.0) {
            std::vector<double> dg(n_variables_);
            for (int i = 0; i < n_variables_; ++i) {
                dg[i] = gradient_new[i] - gradient[i];
            }
            update_inverse_hessian(step, dg);
        }
        
        // Prepare for next iteration
        coords_prev = coords;
        gradient_prev = gradient;
        gradient = gradient_new;
        energy_prev = energy;
    }
    
    // Optimization failed to converge
    result.converged = false;
    result.iterations = conv_.max_iterations;
    result.final_energy = energy_prev;
    result.final_gradient = gradient;
    result.coordinates = coords;
    result.message = "Maximum iterations reached without convergence";
    
    if (verbose_) {
        std::cout << "\n*** Optimization failed to converge in " << conv_.max_iterations << " cycles ***\n";
    }
    
    return result;
}

// ============================================================================
// Scanner Implementation
// ============================================================================

std::pair<double, std::vector<double>> Scanner::operator()(
    const std::vector<double>& coords
) {
    if (!mol_) {
        throw std::runtime_error("Scanner: molecule not set");
    }
    
    mol_->set_coordinates(coords);
    scf_->compute();
    
    double energy = scf_->get_total_energy();
    std::vector<double> gradient = compute_numerical_gradient(scf_, mol_, coords);
    
    return {energy, gradient};
}

std::vector<double> Scanner::compute_numerical_gradient(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol,
    const std::vector<double>& coords,
    double h
) {
    int n = coords.size();
    std::vector<double> gradient(n, 0.0);
    
    // Numerical differentiation using central difference
    // g_i = (E(x+h) - E(x-h)) / (2*h)
    for (int i = 0; i < n; ++i) {
        // Evaluate at x + h
        std::vector<double> coords_plus = coords;
        coords_plus[i] += h;
        mol->set_coordinates(coords_plus);
        double energy_plus = scf->compute().energy;
        
        // Evaluate at x - h
        std::vector<double> coords_minus = coords;
        coords_minus[i] -= h;
        mol->set_coordinates(coords_minus);
        double energy_minus = scf->compute().energy;
        
        // Restore original coordinate
        std::vector<double> coords_restore = coords;
        mol->set_coordinates(coords_restore);
        
        gradient[i] = (energy_plus - energy_minus) / (2.0 * h);
    }
    
    return gradient;
}

// ============================================================================
// Nuclear Gradient (stub - requires derivative integrals)
// ============================================================================

std::vector<double> compute_nuclear_gradient(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol
) {
    // This is a placeholder - full implementation would require
    // analytical gradients from derivative integrals (int1e_ipnuc, etc.)
    // For now, we use numerical gradients from Scanner
    int n_atoms = mol->num_atoms();
    return std::vector<double>(3 * n_atoms, 0.0);
}

} // namespace geom
} // namespace pyscf
