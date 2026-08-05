#ifndef PYSCF_CPP_LEBEDEV_GRID_H
#define PYSCF_CPP_LEBEDEV_GRID_H

#include <vector>
#include <cmath>
#include <memory>

namespace pyscf {

// Forward declaration for Molecule
class Molecule;

namespace dft {

// Grid point structure
struct GridPoint {
    double x, y, z;  // Cartesian coordinates
    double w;        // weight
};

// Radial grid types
enum class RadialGrid { MURA_KNOWLES, GAUSS_CHEBYSHEV, TREUTLER };

// Lebedev grid generator
class LebedevGrid {
public:
    LebedevGrid() = default;
    
    // Generate Lebedev grid for a given number of points
    static std::vector<GridPoint> generate(int n_points);
    
    // Get all available Lebedev grid sizes
    static std::vector<int> available_sizes();
    
private:
    // Internal Lebedev grid generation using published formulas
    static void generate_lebedev_003(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_005(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_007(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_011(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_013(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_015(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_017(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_019(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_023(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_029(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_035(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_041(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_047(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_053(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_059(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_065(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_071(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_077(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_083(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_089(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_095(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_101(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_107(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_113(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_119(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_125(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_131(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_137(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_143(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_149(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_155(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_161(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_167(int n, std::vector<GridPoint>& grid);
    static void generate_lebedev_173(int n, std::vector<GridPoint>& grid);
};

// Radial grid generator
std::vector<double> generate_radial_grid_mura_knowles(int n_points, double r_max = 30.0);
std::vector<double> generate_radial_grid_gauss_chebyshev(int n_points, double r_max = 30.0);

// Pruning schemes
enum class PruneScheme { NONE, TREUTLER, STRATMANN };

// Atomic grid generator
std::vector<GridPoint> generate_atomic_grid(
    int atomic_number,
    int n_radial,
    int n_angular,
    double r_scale = 1.0,
    PruneScheme prune = PruneScheme::NONE
);

// Molecular grid generator
// Note: Molecule is defined in pyscf namespace, accessed via forward declaration
std::vector<GridPoint> generate_molecular_grid(
    const Molecule& mol,
    int n_radial = 50,
    int n_angular = 110,
    PruneScheme prune = PruneScheme::TREUTLER
);

} // namespace dft
} // namespace pyscf

#endif // PYSCF_CPP_LEBEDEV_GRID_H
