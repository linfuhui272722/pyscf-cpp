#ifndef PYSCF_CPP_MOLECULE_H
#define PYSCF_CPP_MOLECULE_H

#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <map>

namespace pyscf {

// Physical constants (atomic units)
constexpr double BOHR_TO_ANGSTROM = 0.529177210903;
constexpr double ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM;
constexpr double PI = 3.14159265358979323846;
constexpr double FPI = 3.14159265358979323846;

// Element data
struct ElementData {
    int atomic_number;
    double atomic_mass;
    std::string symbol;
};

// Basis function types
enum class GTOType { S, P, D, F };

struct GTOContracted {
    GTOType type;
    std::vector<double> exponents;
    std::vector<double> coefficients;
};

// Shell information
struct Shell {
    int atom_index;
    GTOType type;
    int l, m, n;  // angular momentum components
    std::vector<double> exponents;
    std::vector<double> contractions;  // contraction coefficients
};

// Atom structure
struct Atom {
    int atomic_number;
    double x, y, z;  // coordinates in Bohr
};

// Molecule class
class Molecule {
public:
    Molecule();
    ~Molecule();
    
    void add_atom(int atomic_number, double x, double y, double z);
    void set_basis(const std::string& basis_name);
    void load_basis(const std::string& filename);
    
    int num_atoms() const { return atoms_.size(); }
    int num_electrons() const;
    int num_shells() const { return shells_.size(); }
    int num_basis_functions() const;
    
    const Atom& get_atom(int i) const { return atoms_[i]; }
    Atom& get_atom(int i) { return atoms_[i]; }
    const Shell& get_shell(int i) const { return shells_[i]; }
    const std::vector<Shell>& get_shells() const { return shells_; }
    std::vector<Shell>& get_shells() { return shells_; }

    // Set coordinates from a flat array [x1, y1, z1, x2, y2, z2, ...]
    void set_coordinates(const std::vector<double>& coords) {
        for (size_t i = 0; i < atoms_.size() && i * 3 + 2 < coords.size(); ++i) {
            atoms_[i].x = coords[i * 3];
            atoms_[i].y = coords[i * 3 + 1];
            atoms_[i].z = coords[i * 3 + 2];
        }
    }
    
    double total_nuclear_charge() const;
    double nuclear_repulsion_energy() const;
    
    std::vector<double> get_coordinates() const;
    
private:
    std::vector<Atom> atoms_;
    std::vector<Shell> shells_;
    std::string basis_name_;
    
    void parse_basis_file(const std::string& content);
    void parse_polarization_file(const std::string& content, bool add_to_all_atoms);
    int shell_to_ao_offset(int shell) const;
};

// Get element symbol from atomic number
std::string element_symbol(int Z);

// Get atomic number from symbol
int element_atomic_number(const std::string& symbol);

// Common element data
extern const std::map<int, ElementData> ELEMENT_DATA;

} // namespace pyscf

#endif // PYSCF_CPP_MOLECULE_H
