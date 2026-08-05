#include "molecule.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <iostream>

namespace pyscf {

// Element data table
const std::map<int, ElementData> ELEMENT_DATA = {
    {1, {1, 1.00794, "H"}},
    {2, {2, 4.002602, "He"}},
    {3, {3, 6.941, "Li"}},
    {4, {4, 9.012182, "Be"}},
    {5, {5, 10.811, "B"}},
    {6, {6, 12.0107, "C"}},
    {7, {7, 14.0067, "N"}},
    {8, {8, 15.9994, "O"}},
    {9, {9, 18.9984032, "F"}},
    {10, {10, 20.1797, "Ne"}},
    {11, {11, 22.98976928, "Na"}},
    {12, {12, 24.3050, "Mg"}},
    {13, {13, 26.9815386, "Al"}},
    {14, {14, 28.0855, "Si"}},
    {15, {15, 30.973762, "P"}},
    {16, {16, 32.065, "S"}},
    {17, {17, 35.453, "Cl"}},
    {18, {18, 39.948, "Ar"}},
    {19, {19, 39.0983, "K"}},
    {20, {20, 40.078, "Ca"}},
    {36, {36, 83.798, "Kr"}},
    {53, {53, 126.90447, "I"}},
};

std::string element_symbol(int Z) {
    auto it = ELEMENT_DATA.find(Z);
    if (it != ELEMENT_DATA.end()) {
        return it->second.symbol;
    }
    return "X";
}

int element_atomic_number(const std::string& symbol) {
    for (const auto& pair : ELEMENT_DATA) {
        if (pair.second.symbol == symbol) {
            return pair.first;
        }
    }
    return 0;
}

Molecule::Molecule() {}

Molecule::~Molecule() {}

void Molecule::add_atom(int atomic_number, double x, double y, double z) {
    Atom atom;
    atom.atomic_number = atomic_number;
    atom.x = x * ANGSTROM_TO_BOHR;  // Convert from Angstrom to Bohr
    atom.y = y * ANGSTROM_TO_BOHR;
    atom.z = z * ANGSTROM_TO_BOHR;
    atoms_.push_back(atom);
}

int Molecule::num_electrons() const {
    int total_charge = 0;
    for (const auto& atom : atoms_) {
        total_charge += atom.atomic_number;
    }
    // For neutral molecule, electrons = nuclear charge
    return total_charge;
}

int Molecule::num_basis_functions() const {
    int nbf = 0;
    for (const auto& shell : shells_) {
        switch (shell.type) {
            case GTOType::S: nbf += 1; break;
            case GTOType::P: nbf += 3; break;
            case GTOType::D: nbf += 6; break;
            case GTOType::F: nbf += 10; break;
        }
    }
    return nbf;
}

double Molecule::total_nuclear_charge() const {
    double charge = 0.0;
    for (const auto& atom : atoms_) {
        charge += atom.atomic_number;
    }
    return charge;
}

double Molecule::nuclear_repulsion_energy() const {
    double energy = 0.0;
    for (size_t i = 0; i < atoms_.size(); ++i) {
        for (size_t j = i + 1; j < atoms_.size(); ++j) {
            double dx = atoms_[i].x - atoms_[j].x;
            double dy = atoms_[i].y - atoms_[j].y;
            double dz = atoms_[i].z - atoms_[j].z;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            energy += atoms_[i].atomic_number * atoms_[j].atomic_number / r;
        }
    }
    return energy;
}

std::vector<double> Molecule::get_coordinates() const {
    std::vector<double> coords(3 * atoms_.size());
    for (size_t i = 0; i < atoms_.size(); ++i) {
        coords[3*i] = atoms_[i].x;
        coords[3*i+1] = atoms_[i].y;
        coords[3*i+2] = atoms_[i].z;
    }
    return coords;
}

void Molecule::set_basis(const std::string& basis_name) {
    basis_name_ = basis_name;
    shells_.clear();  // Clear any existing shells
    
    // Parse the basis name to extract base and polarization
    // Examples: "6-31G*", "6-31G(d)", "6-31G**", "6-31G(d,p)"
    std::string base_basis = basis_name;
    bool has_d_pol = false;
    bool has_p_pol = false;
    bool has_diffuse = false;
    
    // Check for polarization markers (use uppercase for comparison)
    std::string upper_name = basis_name;
    std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(), ::toupper);
    
    // Handle "6-31G*" and "6-31G**" notation
    if (upper_name.find("*") != std::string::npos) {
        has_d_pol = true;
        base_basis = upper_name;
        // Remove asterisks
        base_basis.erase(std::remove(base_basis.begin(), base_basis.end(), '*'), base_basis.end());
        // Check for P polarization (double asterisk ** means P polarization on H)
        if (upper_name.find("**") != std::string::npos) {
            has_p_pol = true;
        }
    } else {
        base_basis = upper_name;
    }
    
    // Handle "(d)" notation
    if (upper_name.find("(D)") != std::string::npos || upper_name.find("(D,P)") != std::string::npos) {
        has_d_pol = true;
        base_basis = upper_name;
        size_t pos = base_basis.find("(D)");
        if (pos != std::string::npos) base_basis.erase(pos, 3);
        pos = base_basis.find("(D,P)");
        if (pos != std::string::npos) {
            base_basis.erase(pos, 5);
            has_p_pol = true;
        }
    }
    
    // Handle "(p)" notation  
    if (upper_name.find("(P)") != std::string::npos) {
        has_p_pol = true;
        base_basis = upper_name;
        size_t pos = base_basis.find("(P)");
        if (pos != std::string::npos) base_basis.erase(pos, 3);
    }
    
    // Handle diffuse functions
    if (upper_name.find("++") != std::string::npos) {
        has_diffuse = true;
    }
    
    // Remove + signs for base basis
    base_basis.erase(std::remove(base_basis.begin(), base_basis.end(), '+'), base_basis.end());
    
    // Determine actual file names for Pople basis sets
    // 6-31G.dat has only S shells, 6-31Gs.dat has S+P shells
    std::string base_filename = base_basis;
    if (base_basis == "6-31G" || base_basis == "631G") {
        base_filename = "6-31G";
    } else if (base_basis == "3-21G" || base_basis == "321G") {
        base_filename = "3-21G";
    } else if (base_basis == "STO-3G" || base_basis == "STO3G") {
        base_filename = "sto-3g";
    } else if (base_basis == "6-311G") {
        base_filename = "6-311G";
    }
    
    
    // Find paths for base basis - prioritize absolute paths
    std::vector<std::string> base_paths = {
        "/workspace/project/pyscf_cpp/basis/" + base_filename + ".dat",
        "basis/" + base_filename + ".dat",
        "../basis/" + base_filename + ".dat",
        "/workspace/project/pyscf/pyscf/gto/basis/" + base_filename + ".dat",
        "/workspace/project/pyscf/pyscf/gto/basis/pople-basis/" + base_filename + ".dat",
    };
    
    // Load base basis
    bool loaded = false;
    for (const auto& path : base_paths) {
        std::ifstream file(path);
        if (file.is_open()) {
            std::stringstream buffer;
            buffer << file.rdbuf();
            std::string content = buffer.str();
            parse_basis_file(content);
            file.close();
            if (!shells_.empty()) {
                loaded = true;
                break;
            }
        }
    }
    
    if (!loaded) {
        std::cerr << "WARNING: Failed to load base basis: " << base_filename << "\n";
    }
    
    // Add polarization functions if needed
    // For Pople basis sets:
    // - 6-31G*: D polarization on heavy atoms (Z >= 3), no P polarization
    // - 6-31G**: D polarization on heavy atoms + P polarization on all atoms (including H)
    
    // For 6-31G*: add d only to non-H atoms (Z > 2)
    // For 6-31G**: add d to all atoms (including H)
    bool add_to_all_atoms = false;  // true for 6-31G**, false for 6-31G*
    
    if (has_d_pol && (base_basis.find("6-31") != std::string::npos || base_basis.find("6-311") != std::string::npos)) {
        // Find the polarization-d file
        std::string pol_file_name = (base_basis == "6-311G") ? "6-311G-polarization-d.dat" : "6-31G-polarization-d.dat";
        std::vector<std::string> pol_paths = {
            "/workspace/project/pyscf_cpp/basis/" + pol_file_name,
            "basis/" + pol_file_name,
            "/workspace/project/pyscf/pyscf/gto/basis/pople-basis/" + pol_file_name,
        };
        
        for (const auto& path : pol_paths) {
            std::ifstream pol_file(path);
            if (pol_file.is_open()) {
                std::stringstream buffer;
                buffer << pol_file.rdbuf();
                parse_polarization_file(buffer.str(), add_to_all_atoms);
                pol_file.close();
                break;
            }
        }
    }
    
    // Add p polarization if needed (for 6-31G** style)
    if (has_p_pol && (base_basis.find("6-31") != std::string::npos || base_basis.find("6-311") != std::string::npos)) {
        std::string pol_file_name = (base_basis == "6-311G") ? "6-311G-polarization-p.dat" : "6-31G-polarization-p.dat";
        std::vector<std::string> pol_paths = {
            "/workspace/project/pyscf_cpp/basis/" + pol_file_name,
            "basis/" + pol_file_name,
            "/workspace/project/pyscf/pyscf/gto/basis/pople-basis/" + pol_file_name,
        };
        
        for (const auto& path : pol_paths) {
            std::ifstream pol_file(path);
            if (pol_file.is_open()) {
                std::stringstream buffer;
                buffer << pol_file.rdbuf();
                parse_polarization_file(buffer.str(), true);  // Add p to all atoms
                pol_file.close();
                break;
            }
        }
    }
}

void Molecule::parse_basis_file(const std::string& content) {
    std::istringstream iss(content);
    std::string line;
    std::string current_element;
    int current_atomic = 0;
    GTOType current_type = GTOType::S;
    int current_l = 0, current_m = 0, current_n = 0;
    Shell current_shell, p_shell_for_sp;
    bool in_sp_shell = false;
    bool need_this_atom = false;
    bool prev_need_this_atom = false;
    
    while (std::getline(iss, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Trim leading spaces
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start != std::string::npos) {
            line = line.substr(start);
        }
        if (line.empty()) continue;
        
        // Skip header/footer
        if (line.substr(0, 5) == "BASIS" || line.substr(0, 3) == "END") {
            continue;
        }
        
        std::istringstream line_iss(line);
        std::string token;
        line_iss >> token;
        
        // Check if it's an element symbol
        if (token == "H" || token == "He" || token == "Li" || token == "Be" ||
            token == "B" || token == "C" || token == "N" || token == "O" ||
            token == "F" || token == "Ne" || token == "Na" || token == "Mg" ||
            token == "Al" || token == "Si" || token == "P" || token == "S" ||
            token == "Cl" || token == "Ar" || token == "K" || token == "Ca" ||
            token == "Kr" || token == "I") {
            
            // Check if current atom had shells to add (before resetting)
            if (need_this_atom && current_shell.exponents.size() > 0) {
                // Create a shell for each atom of this type
                for (int i = 0; i < num_atoms(); ++i) {
                    if (atoms_[i].atomic_number == current_atomic) {
                        Shell s = current_shell;
                        s.atom_index = i;
                        shells_.push_back(s);
                        if (in_sp_shell && p_shell_for_sp.exponents.size() > 0) {
                            Shell p = p_shell_for_sp;
                            p.atom_index = i;
                            shells_.push_back(p);
                        }
                    }
                }
            }
            
            current_element = token;
            current_atomic = element_atomic_number(token);
            current_shell = Shell();
            p_shell_for_sp = Shell();
            in_sp_shell = false;
            
            // Check if this atom is in our molecule
            need_this_atom = false;
            for (const auto& atom : atoms_) {
                if (atom.atomic_number == current_atomic) {
                    need_this_atom = true;
                    break;
                }
            }
            
            // DEBUG
            
            // Read type (S, SP, P, D, etc.)
            std::string type_str;
            line_iss >> type_str;
            
            if (type_str == "SP") {
                in_sp_shell = true;
                current_type = GTOType::S;
                current_l = current_m = current_n = 0;
                current_shell.type = GTOType::S;
                current_shell.l = 0; current_shell.m = 0; current_shell.n = 0;
                p_shell_for_sp.type = GTOType::P;
                p_shell_for_sp.l = 1; p_shell_for_sp.m = 0; p_shell_for_sp.n = 0;
            } else if (type_str == "S") {
                current_type = GTOType::S;
                current_l = current_m = current_n = 0;
                current_shell.type = GTOType::S;
                current_shell.l = 0; current_shell.m = 0; current_shell.n = 0;
            } else if (type_str == "P") {
                current_type = GTOType::P;
                current_l = 1; current_m = 0; current_n = 0;
                current_shell.type = GTOType::P;
                current_shell.l = 1; current_shell.m = 0; current_shell.n = 0;
            } else if (type_str == "D") {
                current_type = GTOType::D;
                current_l = 2; current_m = 0; current_n = 0;
                current_shell.type = GTOType::D;
                current_shell.l = 2; current_shell.m = 0; current_shell.n = 0;
            } else if (type_str == "F") {
                current_type = GTOType::F;
                current_l = 3; current_m = 0; current_n = 0;
                current_shell.type = GTOType::F;
                current_shell.l = 3; current_shell.m = 0; current_shell.n = 0;
            }
            continue;
        }
        
        // Skip if this atom is not needed
        if (!need_this_atom) {
            continue;
        }
        
        // Parse exponents and coefficients - create fresh stream
        std::istringstream data_iss(line);
        double exponent, coeff1, coeff2 = 0.0;
        if (data_iss >> exponent >> coeff1) {
            data_iss >> coeff2;  // May be zero for non-SP shells
            
            // Add primitive GTO to current shell
            current_shell.exponents.push_back(exponent);
            current_shell.contractions.push_back(coeff1);
            
            // For SP shells, also add to P shell
            if (in_sp_shell) {
                p_shell_for_sp.exponents.push_back(exponent);
                p_shell_for_sp.contractions.push_back(coeff2);
            }
        }
    }
    
    // Don't forget the last atom's shells
    if (need_this_atom && current_shell.exponents.size() > 0) {
        for (int i = 0; i < num_atoms(); ++i) {
            if (atoms_[i].atomic_number == current_atomic) {
                Shell s = current_shell;
                s.atom_index = i;
                shells_.push_back(s);
                if (in_sp_shell && p_shell_for_sp.exponents.size() > 0) {
                    Shell p = p_shell_for_sp;
                    p.atom_index = i;
                    shells_.push_back(p);
                }
            }
        }
    }
}

void Molecule::parse_polarization_file(const std::string& content, bool add_to_all_atoms) {
    // Parse polarization basis set and add shells to appropriate atoms
    // For 6-31G*: add d to non-hydrogen atoms only (Z >= 3 for Li-Ca, or Z > 2 for heavier)
    // For 6-31G**: add d to all atoms
    
    std::istringstream iss(content);
    std::string line;
    int current_atomic = 0;
    Shell current_shell;
    bool in_shell = false;
    
    while (std::getline(iss, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        // Trim leading spaces
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start != std::string::npos) {
            line = line.substr(start);
        }
        if (line.empty()) continue;
        
        // Skip header/footer
        if (line.substr(0, 5) == "BASIS" || line.substr(0, 3) == "END") {
            continue;
        }
        
        std::istringstream line_iss(line);
        std::string token;
        line_iss >> token;
        
        // Check if it's an element symbol
        if (token == "H" || token == "He" || token == "Li" || token == "Be" ||
            token == "B" || token == "C" || token == "N" || token == "O" ||
            token == "F" || token == "Ne" || token == "Na" || token == "Mg" ||
            token == "Al" || token == "Si" || token == "P" || token == "S" ||
            token == "Cl" || token == "Ar" || token == "K" || token == "Ca") {
            
            // Add previous shell if needed
            if (in_shell && current_shell.exponents.size() > 0) {
                // Determine if we should add to this atom type
                bool should_add = false;
                
                if (add_to_all_atoms) {
                    should_add = true;
                } else {
                    // For 6-31G* style: only add to non-hydrogen atoms (Z > 2)
                    // Note: H and He don't get d polarization in 6-31G*
                    if (current_atomic > 2) {
                        should_add = true;
                    }
                }
                
                if (should_add) {
                    for (int i = 0; i < num_atoms(); ++i) {
                        if (atoms_[i].atomic_number == current_atomic) {
                            Shell s = current_shell;
                            s.atom_index = i;
                            shells_.push_back(s);
                        }
                    }
                }
            }
            
            // Get atomic number
            current_atomic = element_atomic_number(token);
            
            // Read shell type (should be D for polarization files)
            std::string type_str;
            line_iss >> type_str;
            
            current_shell = Shell();
            current_shell.atom_index = -1;
            in_shell = true;
            
            if (type_str == "D") {
                current_shell.type = GTOType::D;
                current_shell.l = 2; current_shell.m = 0; current_shell.n = 0;
            } else if (type_str == "P") {
                current_shell.type = GTOType::P;
                current_shell.l = 1; current_shell.m = 0; current_shell.n = 0;
            } else {
                in_shell = false;
            }
            
            continue;
        }
        
        // Skip if not in a shell
        if (!in_shell) continue;
        
        // Parse exponents and coefficients
        std::istringstream data_iss(line);
        double exponent, coeff = 0.0;
        if (data_iss >> exponent >> coeff) {
            current_shell.exponents.push_back(exponent);
            current_shell.contractions.push_back(coeff);
        }
    }
    
    // Don't forget the last shell
    if (in_shell && current_shell.exponents.size() > 0) {
        bool should_add = false;
        
        if (add_to_all_atoms) {
            should_add = true;
        } else {
            // For 6-31G* style: only add to non-hydrogen atoms
            if (current_atomic > 2) {
                should_add = true;
            }
        }
        
        if (should_add) {
            for (int i = 0; i < num_atoms(); ++i) {
                if (atoms_[i].atomic_number == current_atomic) {
                    Shell s = current_shell;
                    s.atom_index = i;
                    shells_.push_back(s);
                }
            }
        }
    }
}

} // namespace pyscf

// Debug helper
namespace {
void debug_shells(const pyscf::Molecule& mol) {
}
}
