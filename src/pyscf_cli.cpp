/**
 * PySCF C++ - Universal Command-line Interface
 * Supports any molecule input (XYZ, PDB, coordinate format)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <algorithm>
#include "pyscf_cpp.h"

using namespace pyscf;

// Parse XYZ format
std::shared_ptr<Molecule> parse_xyz(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        // Try as direct coordinate input
        return nullptr;
    }
    
    std::string line;
    int natoms = 0;
    
    // Read number of atoms
    std::getline(file, line);
    natoms = std::stoi(line);
    
    // Skip comment line
    std::getline(file, line);
    
    auto mol = std::make_shared<Molecule>();
    
    for (int i = 0; i < natoms; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        std::string element;
        double x, y, z;
        iss >> element >> x >> y >> z;
        
        int atomic_num = element_atomic_number(element);
        if (atomic_num > 0) {
            // Convert Angstrom to Bohr
            mol->add_atom(atomic_num, x * 1.889726, y * 1.889726, z * 1.889726);
        }
    }
    
    return mol;
}

// Parse coordinate string
std::shared_ptr<Molecule> parse_coordinates(const std::string& input) {
    std::istringstream iss(input);
    std::string line;
    
    auto mol = std::make_shared<Molecule>();
    bool found_coords = false;
    
    while (std::getline(iss, line)) {
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t") == std::string::npos) continue;
        
        // Try to parse as: Element X Y Z
        std::istringstream line_ss(line);
        std::string element;
        double x, y, z;
        
        if (line_ss >> element >> x >> y >> z) {
            // Check if first token is an element symbol
            int atomic_num = element_atomic_number(element);
            if (atomic_num > 0) {
                // Convert Angstrom to Bohr
                mol->add_atom(atomic_num, x * 1.889726, y * 1.889726, z * 1.889726);
                found_coords = true;
            }
        }
    }
    
    return found_coords ? mol : nullptr;
}

// Print usage
void print_usage(const char* program) {
    std::cout << "PySCF C++ - Universal Quantum Chemistry Calculator\n";
    std::cout << "==================================================\n\n";
    std::cout << "Usage: " << program << " [options] [input]\n\n";
    std::cout << "Input options:\n";
    std::cout << "  -f <file>     Input file (XYZ format)\n";
    std::cout << "  -c <coords>   Inline coordinates (Element X Y Z per line)\n";
    std::cout << "  -i            Interactive mode\n\n";
    std::cout << "Calculation options:\n";
    std::cout << "  -m <method>   Method: rks, uks (default: rks)\n";
    std::cout << "  -x <xc>       XC functional (default: b3lyp)\n";
    std::cout << "  -b <basis>    Basis set (default: sto-3g)\n";
    std::cout << "  -t <tol>      SCF convergence tolerance (default: 1e-6)\n";
    std::cout << "  -n <iter>     Max SCF iterations (default: 100)\n\n";
    std::cout << "Analysis options:\n";
    std::cout << "  -o            Run geometry optimization\n";
    std::cout << "  -B            Compute Mayer bond orders\n";
    std::cout << "  -v            Verbose output\n";
    std::cout << "  -h            Show this help\n\n";
    std::cout << "Examples:\n";
    std::cout << "  " << program << " -f water.xyz -b cc-pvdz\n";
    std::cout << "  " << program << " -c \"O 0 0 0\\nH 0.96 0 0\\nH 0 0.96 70.5\"\n\n";
}

// Interactive input
std::shared_ptr<Molecule> interactive_input() {
    std::cout << "\n=== Interactive Molecule Input ===\n";
    std::cout << "Enter atoms as: Element X Y Z (in Angstrom)\n";
    std::cout << "Enter empty line when done, 'quit' to exit.\n\n";
    
    auto mol = std::make_shared<Molecule>();
    std::string line;
    
    while (true) {
        std::cout << "> ";
        std::getline(std::cin, line);
        
        if (line == "quit" || line == "q") return nullptr;
        if (line.empty()) break;
        
        std::istringstream iss(line);
        std::string element;
        double x, y, z;
        
        if (iss >> element >> x >> y >> z) {
            int atomic_num = element_atomic_number(element);
            if (atomic_num > 0) {
                mol->add_atom(atomic_num, x * 1.889726, y * 1.889726, z * 1.889726);
                std::cout << "  Added " << element << " at (" << x << ", " << y << ", " << z << ")\n";
            } else {
                std::cerr << "  Unknown element: " << element << "\n";
            }
        } else {
            std::cerr << "  Invalid format. Use: Element X Y Z\n";
        }
    }
    
    return mol->num_atoms() > 0 ? mol : nullptr;
}

int main(int argc, char** argv) {
    std::string input_file, coordinates, method = "rks", xc = "b3lyp", basis = "sto-3g";
    double conv_tol = 1e-6;
    int max_iter = 100;
    bool do_optimize = false, do_bond = false, verbose = false, interactive = false;
    
    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-f" && i + 1 < argc) {
            input_file = argv[++i];
        } else if (arg == "-c" && i + 1 < argc) {
            coordinates = argv[++i];
        } else if (arg == "-i") {
            interactive = true;
        } else if (arg == "-m" && i + 1 < argc) {
            method = argv[++i];
        } else if (arg == "-x" && i + 1 < argc) {
            xc = argv[++i];
        } else if (arg == "-b" && i + 1 < argc) {
            basis = argv[++i];
        } else if (arg == "-t" && i + 1 < argc) {
            conv_tol = std::stod(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            max_iter = std::stoi(argv[++i]);
        } else if (arg == "-o") {
            do_optimize = true;
        } else if (arg == "-B") {
            do_bond = true;
        } else if (arg == "-v") {
            verbose = true;
        } else if (arg == "-h") {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    std::cout << "==================================================\n";
    std::cout << "   PySCF C++ - Universal Quantum Chemistry\n";
    std::cout << "==================================================\n\n";
    
    std::shared_ptr<Molecule> mol;
    
    // Get molecule
    if (interactive) {
        mol = interactive_input();
        if (!mol) {
            std::cout << "No molecule entered. Exiting.\n";
            return 0;
        }
    } else if (!input_file.empty()) {
        mol = parse_xyz(input_file);
        if (!mol) {
            std::cerr << "Error: Could not read file: " << input_file << "\n";
            return 1;
        }
    } else if (!coordinates.empty()) {
        mol = parse_coordinates(coordinates);
        if (!mol || mol->num_atoms() == 0) {
            std::cerr << "Error: Could not parse coordinates.\n";
            return 1;
        }
    } else {
        std::cerr << "Error: No input specified. Use -h for help.\n";
        return 1;
    }
    
    // Load basis set
    std::cout << "Loading basis set: " << basis << "\n";
    mol->set_basis(basis);
    
    std::cout << "\n--- Molecule Summary ---\n";
    std::cout << "Atoms: " << mol->num_atoms() << "\n";
    std::cout << "Electrons: " << mol->num_electrons() << "\n";
    std::cout << "Basis functions: " << mol->num_basis_functions() << "\n";
    std::cout << "Nuclear repulsion: " << mol->nuclear_repulsion_energy() << " Eh\n";
    
    if (mol->num_basis_functions() == 0) {
        std::cerr << "\nError: No basis functions loaded!\n";
        return 1;
    }
    
    // Run SCF calculation
    std::cout << "\n--- SCF Calculation ---\n";
    std::cout << "Method: " << method << "\n";
    std::cout << "XC Functional: " << xc << "\n";
    
    auto scf = std::make_shared<dft::RKS>(mol);
    scf->set_xc_functional(xc);
    
    dft::SCFOptions options;
    options.conv_tol = conv_tol;
    options.max_cycle = max_iter;
    options.verbose = verbose;
    scf->set_options(options);
    
    auto result = scf->compute();
    
    std::cout << "\n--- Results ---\n";
    std::cout << "Converged: " << (result.converged ? "Yes" : "No") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Total energy: " << result.energy << " Eh\n";
    std::cout << "             " << result.energy * 27.2114 << " eV\n";
    std::cout << "Electronic energy: " << scf->get_electronic_energy() << " Eh\n";
    std::cout << "Nuclear repulsion: " << scf->get_nuclear_repulsion() << " Eh\n";
    
    // Geometry optimization
    if (do_optimize) {
        std::cout << "\n--- Geometry Optimization ---\n";
        geom::BFGSOptimizer opt;
        geom::OptConvergence conv;
        conv.max_iterations = 50;
        opt.set_convergence(conv);
        auto opt_result = opt.optimize(scf, mol);
        
        std::cout << "Optimized energy: " << opt_result.final_energy << " Eh\n";
        std::cout << "Converged: " << (opt_result.converged ? "Yes" : "No") << "\n";
        std::cout << "Iterations: " << opt_result.iterations << "\n\n";
        std::cout << "Optimized coordinates (Angstrom):\n";
        for (int i = 0; i < mol->num_atoms(); ++i) {
            const auto& atom = mol->get_atom(i);
            std::cout << "  " << element_symbol(atom.atomic_number)
                      << "  " << (atom.x / 1.889726) << "  "
                      << (atom.y / 1.889726) << "  "
                      << (atom.z / 1.889726) << "\n";
        }
    }
    
    // Mayer bond order analysis
    if (do_bond) {
        std::cout << "\n--- Mayer Bond Order Analysis ---\n";
        auto bond_result = bond::analyze_bonds(scf, mol);
        
        for (const auto& b : bond_result.bonds) {
            int i, j;
            double order;
            std::tie(i, j, order) = b;
            std::cout << "  " << element_symbol(mol->get_atom(i).atomic_number) << "-" 
                      << element_symbol(mol->get_atom(j).atomic_number)
                      << " (" << i << "-" << j << "): "
                      << order << "\n";
        }
    }
    
    std::cout << "\n==================================================\n";
    std::cout << "Calculation complete!\n";
    std::cout << "==================================================\n";
    
    return result.converged ? 0 : 1;
}
