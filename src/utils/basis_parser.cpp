#include "basis_parser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace pyscf {
namespace utils {

std::map<int, std::vector<Shell>> load_basis_from_file(const std::string& filename) {
    std::map<int, std::vector<Shell>> shells_by_atom;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        return shells_by_atom;
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();
    
    std::istringstream iss(content);
    std::string line;
    std::string current_element;
    int current_atomic = 0;
    GTOType current_type = GTOType::S;
    
    int line_num = 0;
    
    while (std::getline(iss, line)) {
        line_num++;
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        // Skip header
        if (line.find("BASIS") != std::string::npos || 
            line.find("END") != std::string::npos) {
            continue;
        }
        
        std::istringstream line_iss(line);
        std::string token;
        line_iss >> token;
        
        // Skip empty tokens
        if (token.empty()) continue;
        
        // Check if it's an element symbol
        static const std::vector<std::string> elements = {
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
            "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd"
        };
        
        bool is_element = false;
        for (const auto& elem : elements) {
            if (token == elem) {
                is_element = true;
                current_element = token;
                
                // Get atomic number
                if (elem == "H") current_atomic = 1;
                else if (elem == "He") current_atomic = 2;
                else if (elem == "Li") current_atomic = 3;
                else if (elem == "Be") current_atomic = 4;
                else if (elem == "B") current_atomic = 5;
                else if (elem == "C") current_atomic = 6;
                else if (elem == "N") current_atomic = 7;
                else if (elem == "O") current_atomic = 8;
                else if (elem == "F") current_atomic = 9;
                else if (elem == "Ne") current_atomic = 10;
                else if (elem == "Na") current_atomic = 11;
                else if (elem == "Mg") current_atomic = 12;
                else if (elem == "Al") current_atomic = 13;
                else if (elem == "Si") current_atomic = 14;
                else if (elem == "P") current_atomic = 15;
                else if (elem == "S") current_atomic = 16;
                else if (elem == "Cl") current_atomic = 17;
                else if (elem == "Ar") current_atomic = 18;
                else if (elem == "K") current_atomic = 19;
                else if (elem == "Ca") current_atomic = 20;
                else if (elem == "Zn") current_atomic = 30;
                else if (elem == "Br") current_atomic = 35;
                else if (elem == "Kr") current_atomic = 36;
                else if (elem == "I") current_atomic = 53;
                break;
            }
        }
        
        if (is_element) {
            // Read shell type
            std::string type_str;
            line_iss >> type_str;
            
            if (type_str == "SP") {
                current_type = GTOType::S;  // Will handle both
            } else if (type_str == "S") {
                current_type = GTOType::S;
            } else if (type_str == "P") {
                current_type = GTOType::P;
            } else if (type_str == "D") {
                current_type = GTOType::D;
            } else if (type_str == "F") {
                current_type = GTOType::F;
            }
            
            shells_by_atom[current_atomic] = std::vector<Shell>();
            continue;
        }
        
        // Set l, m, n based on shell type
        int l_val = 0, m_val = 0, n_val = 0;
        if (current_type == GTOType::S) {
            l_val = m_val = n_val = 0;
        }
        
        // Parse exponents and coefficients
        double exponent, coeff1, coeff2 = 0.0;
        if (line_iss >> exponent >> coeff1) {
            line_iss >> coeff2;  // May be zero for non-SP shells
            
            Shell shell;
            shell.atom_index = -1;
            shell.type = current_type;
            shell.l = l_val;
            shell.m = m_val;
            shell.n = n_val;
            shell.exponents.push_back(exponent);
            shell.contractions.push_back(coeff1);
            
            if (current_type == GTOType::S && coeff2 != 0.0) {
                // SP shell - add both S and P contractions
                shells_by_atom[current_atomic].push_back(shell);
                shell.type = GTOType::P;
                shell.l = 1; shell.m = 0; shell.n = 0;  // P shell
                shell.contractions.clear();
                shell.contractions.push_back(coeff2);
            }
            shells_by_atom[current_atomic].push_back(shell);
        }
    }
    
    for (const auto& kv : shells_by_atom) {
    }
    
    file.close();
    return shells_by_atom;
}

BasisSet parse_basis_file(const std::string& filename) {
    BasisSet basis;
    basis.name = filename;
    basis.shells_by_atom = load_basis_from_file(filename);
    return basis;
}

BasisSet load_basis_set(const std::string& name) {
    
    // Try different paths and name variations
    std::vector<std::string> name_variations = {
        name,
        // Convert lowercase to various formats
        "sto-3g"  // special case
    };
    
    // Build path list
    std::vector<std::string> paths = {
        "basis/" + name + ".dat",
        "../basis/" + name + ".dat",
        "/workspace/project/pyscf_cpp/basis/" + name + ".dat",
    };
    
    // Handle special basis set names
    std::string filename;
    if (name == "sto-3g") {
        filename = "sto-3g.dat";
        paths.push_back("/workspace/project/pyscf_cpp/basis/" + filename);
    } else if (name == "6-31g" || name == "6-31G") {
        filename = "6-31Gss.dat";
        paths.push_back("/workspace/project/pyscf_cpp/basis/" + filename);
    } else if (name == "6-31g*" || name == "6-31G*") {
        filename = "6-31G-polarization-d.dat";
        paths.push_back("/workspace/project/pyscf_cpp/basis/" + filename);
    } else {
        filename = name + ".dat";
        paths.push_back("/workspace/project/pyscf_cpp/basis/" + filename);
    }
    
    // Also check PySCF's basis directory
    paths.push_back("/workspace/project/pyscf/pyscf/gto/basis/" + filename);
    
    for (const auto& path : paths) {
        std::ifstream file(path);
        if (file.is_open()) {
            file.close();
            return parse_basis_file(path);
        }
    }
    
    return BasisSet();
}

std::vector<std::string> available_basis_sets() {
    return {
        "sto-3g",
        "3-21g",
        "6-31g",
        "6-31g*",
        "6-31+g*",
        "cc-pvdz",
        "cc-pvtz",
        "cc-pvqz",
        "cc-pv5z",
        "aug-cc-pvdz",
        "aug-cc-pvtz",
        "aug-cc-pvqz",
        "def2-svp",
        "def2-tzvp",
        "def2-qzvp"
    };
}

bool has_basis_set(const std::string& name) {
    auto available = available_basis_sets();
    return std::find(available.begin(), available.end(), name) != available.end();
}

std::vector<Shell> get_basis_for_element(const BasisSet& basis, int atomic_number) {
    auto it = basis.shells_by_atom.find(atomic_number);
    if (it != basis.shells_by_atom.end()) {
        return it->second;
    }
    return std::vector<Shell>();
}

// Compute Gaussian type orbital value at a point
double gto_value(double alpha, int l, int m, int n, 
                 double dx, double dy, double dz) {
    double r2 = dx*dx + dy*dy + dz*dz;
    double pre = std::pow(dx, l) * std::pow(dy, m) * std::pow(dz, n);
    return pre * std::exp(-alpha * r2);
}

// Evaluate contracted GTO at a point
double contracted_gto_value(const Shell& shell, int l, int m, int n,
                          double dx, double dy, double dz) {
    double sum = 0.0;
    for (size_t i = 0; i < shell.exponents.size(); ++i) {
        sum += shell.contractions[i] * 
               gto_value(shell.exponents[i], l, m, n, dx, dy, dz);
    }
    return sum;
}

std::vector<double> eval_ao_at_point(const Molecule& mol, 
                                      double x, double y, double z) {
    std::vector<double> ao_values;
    
    for (int shell_idx = 0; shell_idx < mol.num_shells(); ++shell_idx) {
        const auto& shell = mol.get_shell(shell_idx);
        const auto& atom = mol.get_atom(shell.atom_index);
        
        double dx = x - atom.x;
        double dy = y - atom.y;
        double dz = z - atom.z;
        
        switch (shell.type) {
            case GTOType::S:
                ao_values.push_back(contracted_gto_value(shell, 0, 0, 0, dx, dy, dz));
                break;
            case GTOType::P:
                ao_values.push_back(contracted_gto_value(shell, 1, 0, 0, dx, dy, dz));  // px
                ao_values.push_back(contracted_gto_value(shell, 0, 1, 0, dx, dy, dz));  // py
                ao_values.push_back(contracted_gto_value(shell, 0, 0, 1, dx, dy, dz));  // pz
                break;
            case GTOType::D:
                // 6 Cartesian d functions
                ao_values.push_back(contracted_gto_value(shell, 2, 0, 0, dx, dy, dz));  // dxx
                ao_values.push_back(contracted_gto_value(shell, 1, 1, 0, dx, dy, dz));  // dxy
                ao_values.push_back(contracted_gto_value(shell, 1, 0, 1, dx, dy, dz));  // dxz
                ao_values.push_back(contracted_gto_value(shell, 0, 2, 0, dx, dy, dz));  // dyy
                ao_values.push_back(contracted_gto_value(shell, 0, 1, 1, dx, dy, dz));  // dyz
                ao_values.push_back(contracted_gto_value(shell, 0, 0, 2, dx, dy, dz));  // dzz
                break;
            case GTOType::F:
                // 10 Cartesian f functions
                for (int i = 0; i < 10; ++i) {
                    ao_values.push_back(0.0);  // Placeholder
                }
                break;
        }
    }
    
    return ao_values;
}

void eval_ao_at_grid(const Molecule& mol,
                    const std::vector<double>& coords,
                    std::vector<double>& ao_values) {
    int n_grid = coords.size() / 3;
    int n_ao = mol.num_basis_functions();
    
    ao_values.assign(n_grid * n_ao, 0.0);
    
    for (int i = 0; i < n_grid; ++i) {
        double x = coords[3 * i];
        double y = coords[3 * i + 1];
        double z = coords[3 * i + 2];
        
        auto ao_point = eval_ao_at_point(mol, x, y, z);
        for (size_t j = 0; j < ao_point.size() && (size_t)(3 * i + j) < ao_values.size(); ++j) {
            ao_values[i * n_ao + j] = ao_point[j];
        }
    }
}

// Simplified overlap matrix
std::vector<double> compute_overlap_matrix(const Molecule& mol) {
    int nbf = mol.num_basis_functions();
    std::vector<double> S(nbf * nbf, 0.0);
    
    // Simplified: identity matrix (real code uses integral evaluation)
    for (int i = 0; i < nbf; ++i) {
        S[i * nbf + i] = 1.0;
    }
    
    return S;
}

// Simplified kinetic matrix
std::vector<double> compute_kinetic_matrix(const Molecule& mol) {
    int nbf = mol.num_basis_functions();
    std::vector<double> T(nbf * nbf, 0.0);
    
    // Simplified
    for (int i = 0; i < nbf; ++i) {
        T[i * nbf + i] = 1.0;
    }
    
    return T;
}

// Simplified nuclear matrix
std::vector<double> compute_nuclear_matrix(const Molecule& mol) {
    int nbf = mol.num_basis_functions();
    std::vector<double> V(nbf * nbf, 0.0);
    
    // Simplified: negative of nuclear charge
    double Z = mol.total_nuclear_charge();
    for (int i = 0; i < nbf; ++i) {
        V[i * nbf + i] = -Z;
    }
    
    return V;
}

// H_core = T + V
std::vector<double> compute_h_core(const Molecule& mol) {
    int nbf = mol.num_basis_functions();
    auto T = compute_kinetic_matrix(mol);
    auto V = compute_nuclear_matrix(mol);
    
    std::vector<double> H(nbf * nbf, 0.0);
    for (int i = 0; i < nbf * nbf; ++i) {
        H[i] = T[i] + V[i];
    }
    
    return H;
}

} // namespace utils
} // namespace pyscf
