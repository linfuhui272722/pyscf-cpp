# PySCF-C++

> **⚠️ 重要声明 | Important Disclaimer**
> 
> 本项目由AI编写，不保证安全性、稳定性或科学准确性。
> 本代码仅供学习研究使用，请勿用于生产环境或科学发表。
> 
> **This project was written by AI and does not guarantee safety, stability, or scientific accuracy. This code is for learning and research purposes only. Do not use it for production environments or scientific publications.**

[English](README_en.md) | 中文

---

A C++ implementation of PySCF-like density functional theory (DFT) calculations, featuring:

- **B3LYP Hybrid Functional** - The most widely used exchange-correlation functional
- **Multiple Basis Sets** - STO-3G, 6-31G, 6-31G*, 6-31G**, 3-21G, and more
- **Geometry Optimization** - BFGS-based molecular structure optimization
- **Bond Order Analysis** - Mayer bond orders and Wiberg indices
- **Gaussian-Type Orbitals** - Powered by libcint library

## 项目状态 | Project Status

⚠️ **本项目为实验性项目，不适合用于生产环境或科学研究。**

⚠️ **This is an experimental project and is NOT suitable for production use or scientific research.**

## 快速开始 | Quick Start

### 安装依赖 | Install Dependencies

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y cmake g++ libopenblas-dev liblapack-dev make

# Install libcint (required for integrals)
git clone https://github.com/sunqm/libcint.git
cd libcint && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make -j4 && sudo make install
```

### 编译 | Build

```bash
cd pyscf-cpp
mkdir build && cd build
cmake ..
make -j4
```

### 运行示例 | Run Examples

```bash
# Energy calculation
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH ./pyscf_demo 6-31g

# Bond order analysis
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH ./pyscf_bond_demo

# Geometry optimization
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH ./pyscf_optimize_demo
```

## 功能特性 | Features

### 支持的基组 | Supported Basis Sets

| 基组 | 类型 | 描述 |
|------|------|------|
| STO-3G | 最小基 | Minimal basis set |
| 3-21G | 分裂价 | Split valence |
| 6-31G | 分裂价 | Split valence |
| 6-31G* | 极化 | With d polarization |
| 6-31G** | 双极化 | Double polarization |

### 计算功能 | Calculation Features

- ✅ 能量计算 (Energy Calculation)
- ✅ SCF收敛 (SCF Convergence)
- ✅ B3LYP泛函 (B3LYP Functional)
- ⚠️ 几何优化 (Geometry Optimization) - 实验性
- ⚠️ Mayer键级 (Mayer Bond Orders) - 实验性

## 示例代码 | Example Code

### 基础能量计算 | Basic Energy Calculation

```cpp
#include <iostream>
#include <memory>
#include "molecule.h"
#include "scf.h"

int main() {
    // 创建分子
    auto mol = std::make_shared<pyscf::Molecule>();
    mol->add_atom(8, 0.0, 0.0, 0.0);   // O
    mol->add_atom(1, 1.4, 1.2, 0.0);   // H1
    mol->add_atom(1, -1.4, 1.2, 0.0);  // H2
    mol->set_basis("6-31g");
    
    // 创建SCF计算器
    auto scf = std::make_shared<pyscf::dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    
    // 运行计算
    scf->compute();
    
    // 获取结果
    std::cout << "总能量: " << scf->get_total_energy() << " Eh\n";
    
    return 0;
}
```

### Mayer键级分析 | Mayer Bond Order Analysis

```cpp
#include "mayer_bond.h"

// ...
auto result = pyscf::bond::analyze_bonds(scf, mol);
pyscf::bond::print_bond_analysis(result, *mol);
```

## 项目结构 | Project Structure

```
pyscf-cpp/
├── include/                    # 头文件
│   ├── molecule.h              # 分子结构
│   ├── scf.h                  # SCF计算
│   ├── optimizer.h             # 几何优化
│   ├── mayer_bond.h           # 键级分析
│   └── ...
├── src/                        # 源文件
│   ├── molecule.cpp
│   ├── scf.cpp
│   ├── gto/                    # GTO积分
│   ├── dft/                    # DFT计算
│   ├── geom/                   # 几何优化
│   └── bond/                   # 键级分析
├── basis/                      # 基组数据
├── examples/                   # 示例程序
└── CMakeLists.txt
```

## 文档 | Documentation

详细教程请查看：
- [中文教程](TUTORIAL_zh.md) - 详细的C++编程指南
- [English Tutorial](TUTORIAL_en.md) - Comprehensive C++ guide

## 已知问题 | Known Issues

⚠️ **SCF Convergence**: SCF may produce NaN for certain molecular geometries
⚠️ **Numerical Stability**: Some calculations may be numerically unstable
⚠️ **Limited Basis Sets**: Only a few Pople basis sets are supported
⚠️ **Limited Functionals**: Only B3LYP is implemented

## 限制和免责声明 | Limitations and Disclaimer

### 科学准确性 | Scientific Accuracy

本项目：
- ❌ 不保证计算结果的科学准确性
- ❌ 不保证与商业软件（如Gaussian、ORCA）的一致性
- ❌ 不适合用于任何形式的科学研究或发表
- ❌ 不适合用于任何生产环境

This project:
- ❌ Does NOT guarantee scientific accuracy of results
- ❌ Does NOT guarantee consistency with commercial software
- ❌ Is NOT suitable for any scientific research or publication
- ❌ Is NOT suitable for any production environment

### 用途 | Purpose

本项目仅用于：
- 学习量子化学和DFT的基本原理
- 理解C++编程和数值计算
- 教学演示目的

This project is only for:
- Learning quantum chemistry and DFT basics
- Understanding C++ programming and numerical computing
- Teaching demonstration purposes

## 许可证 | License

本项目采用 **MIT许可证**。

This project is licensed under the **MIT License**.

请在使用前阅读 [LICENSE](LICENSE) 文件。

Please read the [LICENSE](LICENSE) file before using.

## 致谢 | Acknowledgments

本项目参考了以下开源项目：
- [PySCF](https://github.com/pyscf/pyscf) - Python-based quantum chemistry package
- [libcint](https://github.com/sunqm/libcint) - Analytical integral library for Gaussian-type functions
- [OpenBLAS](https://github.com/xianyi/OpenBLAS) - Optimized BLAS library

## 贡献者 | Contributors

*本项目由AI (OpenHands) 编写*

*This project was written by AI (OpenHands)*

## 联系方式 | Contact

本项目不接受任何形式的贡献或支持请求。

This project does not accept any form of contributions or support requests.

---

**最后更新 | Last Updated**: 2026年8月 | August 2026
