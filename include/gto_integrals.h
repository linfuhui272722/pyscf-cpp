#ifndef PYSCF_CPP_GTO_INTEGRALS_H
#define PYSCF_CPP_GTO_INTEGRALS_H

#include <vector>
#include <cmath>
#include <memory>

namespace pyscf {

// Forward declaration for Molecule
class Molecule;

namespace gto {

// IntegralEngine uses libcint for accurate GTO integrals
class IntegralEngine {
public:
    IntegralEngine();
    ~IntegralEngine();
    
    // Initialize with molecule data
    void initialize(const Molecule& mol);
    
    // Compute overlap matrix S
    void compute_overlap(const Molecule& mol, double* S);
    
    // Compute kinetic matrix T
    void compute_kinetic(const Molecule& mol, double* T);
    
    // Compute nuclear attraction matrix V
    void compute_nuclear(const Molecule& mol, double* V);
    
    // Compute two-electron repulsion integrals (ij|kl)
    // Returns pointer to allocated buffer, caller responsible for deletion
    double* compute_eri(const Molecule& mol, size_t& size);
    
    // Compute AO density at grid points
    void eval_ao(const Molecule& mol, const double* grid, int n_points, double* ao_values);
    
private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

// Helper functions
inline double factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; ++i) result *= i;
    return result;
}

inline double double_factorial(int n) {
    if (n <= 0) return 1.0;
    double result = 1.0;
    for (int i = n; i > 0; i -= 2) result *= i;
    return result;
}

} // namespace gto
} // namespace pyscf

#endif // PYSCF_CPP_GTO_INTEGRALS_H
