# PySCF iOS Executable

Universal quantum chemistry calculator for iOS (arm64).

## Installation

1. Download `pyscf` and `basis/` directory
2. Put them in the **same directory**
3. Give execute permission: `chmod +x pyscf`

## Usage

```bash
# Hydrogen molecule (H2)
./pyscf -c "H 0 0 0" -c "H 0.74 0 0"

# Water molecule (H2O)
./pyscf -c "O 0 0 0" -c "H 0.96 0 0" -c "H 0 0.96 70.5"

# Geometry optimization
./pyscf -c "H 0 0 0" -c "H 0.74 0 0" -o

# Show help
./pyscf -h
```

## Basis Sets

**Important**: If basis set name contains `*` or `+`, use quotes:

```bash
# WRONG (zsh will expand *)
./pyscf -c "H 0 0 0" -c "H 0.74 0 0" -b 6-31G*

# CORRECT (quote the name)
./pyscf -c "H 0 0 0" -c "H 0.74 0 0" -b "6-31G*"
./pyscf -c "H 0 0 0" -c "H 0.74 0 0" -b "6-31G**"
```

Available basis sets:
- `sto-3g` (default)
- `6-31g`
- `6-31g*` (use quotes!)
- `6-31g**` (use quotes!)
- `6-31gs`
- `cc-pvdz`
- `cc-pvtz`

## Options

```bash
-b <basis>    Basis set (default: sto-3g)
-o            Run geometry optimization
-B            Compute Mayer bond orders
-x <func>     XC functional (b3lyp, hf, lda)
-t <tol>      SCF convergence tolerance
-n <iter>     Max SCF iterations
-v            Verbose output
-h            Show help
```

## Directory Structure

```
your_folder/
├── pyscf          # Main executable
└── basis/         # Basis set files (required!)
    ├── sto-3g.dat
    ├── 6-31G.dat
    └── ...
```

## Requirements

- iOS 16.5+ or later
- arm64 device (Apple Silicon)
