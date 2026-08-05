# PySCF-C++ 分子计算C++程序编写详细教程

> **⚠️ 免责声明**: 本教程由AI编写，不保证代码的安全性、稳定性或科学准确性。本代码仅供学习研究使用，请勿将其用于生产环境或科学发表。如需进行严肃的量子化学计算，请使用经过验证的成熟软件如Gaussian、ORCA、PySCF等。

## 目录

1. [简介](#1-简介)
2. [开发环境配置](#2-开发环境配置)
3. [项目结构](#3-项目结构)
4. [核心概念](#4-核心概念)
5. [基础示例：水分子的能量计算](#5-基础示例水分子的能量计算)
6. [进阶示例：几何优化](#6-进阶示例几何优化)
7. [进阶示例：Mayer键级分析](#7-进阶示例mayer键级分析)
8. [完整示例程序](#8-完整示例程序)
9. [基组处理](#9-基组处理)
10. [常见问题与调试](#10-常见问题与调试)

---

## 1. 简介

本项目是一个C++实现的密度泛函理论(DFT)计算框架，参考了PySCF的设计思想。主要功能包括：

- **B3LYP泛函**：最常用的混合泛函
- **多种基组**：STO-3G、6-31G、6-31G*、6-31G**、3-21G等
- **几何优化**：基于BFGS算法的分子构型优化
- **键级分析**：Mayer键级和Wiberg指数计算

### 1.1 系统架构

```
┌─────────────────────────────────────────────────────────┐
│                    用户程序 (main.cpp)                   │
├─────────────────────────────────────────────────────────┤
│  Molecule │ RKS (SCF) │ Optimizer │ BondAnalyzer       │
├─────────────────────────────────────────────────────────┤
│                    核心计算层                             │
├─────────────────────────────────────────────────────────┤
│  GTO积分 │ Lebedev积分 │ B3LYP泛函 │ BLAS/LAPACK      │
├─────────────────────────────────────────────────────────┤
│              libcint (积分库) / OpenBLAS                │
└─────────────────────────────────────────────────────────┘
```

---

## 2. 开发环境配置

### 2.1 依赖项

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    cmake \
    g++ \
    libopenblas-dev \
    liblapack-dev \
    libcint-dev \
    make

# 或者使用 apt-get 安装 libcint
sudo apt-get install -y libcint-dev
```

### 2.2 从源码编译libcint

如果系统没有libcint，需要从源码编译：

```bash
git clone https://github.com/sunqm/libcint.git
cd libcint
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make -j4
sudo make install
```

### 2.3 编译项目

```bash
cd pyscf-cpp
mkdir build && cd build
cmake ..
make -j4
```

---

## 3. 项目结构

```
pyscf-cpp/
├── include/                    # 头文件
│   ├── molecule.h              # 分子结构定义
│   ├── scf.h                   # SCF计算接口
│   ├── optimizer.h              # 几何优化器
│   ├── mayer_bond.h            # 键级分析
│   ├── b3lyp_xc.h             # B3LYP泛函
│   └── lebedev_grid.h          # 数值积分网格
├── src/                        # 源文件
│   ├── molecule.cpp
│   ├── scf.cpp
│   ├── gto/                    # GTO积分
│   ├── dft/                    # DFT计算
│   ├── geom/                   # 几何优化
│   └── bond/                   # 键级分析
├── basis/                      # 基组数据文件
├── examples/                   # 示例程序
└── CMakeLists.txt
```

---

## 4. 核心概念

### 4.1 原子单位制

量子化学计算通常使用原子单位(a.u.)：
- 长度：1 Bohr = 0.529177 Å
- 能量：1 Eh (Hartree) = 27.2114 eV

### 4.2 高斯型轨道(GTO)

基函数表示为高斯函数的线性组合：

```
χ(r) = Σ c_i × g_i(r)
g_i(r) = (x-ax)^l (y-ay)^m (z-az)^n × exp(-α_i × |r-R|^2)
```

### 4.3 密度矩阵

对于闭壳层体系，密度矩阵：

```
P_μν = 2 × Σ_i(occ) C_μi × C_νi
```

### 4.4 B3LYP泛函

B3LYP = 20% HF + 72% LDA/GGA + 8% GGA交换

```
E_XC^B3LYP = (1-a) × E_X^LDA + a × E_X^HF 
            + c × ΔE_X^B88 + (1-c) × E_C^VWN 
            + d × ΔE_C^LYP
其中 a=0.20, c=0.72, d=0.81
```

---

## 5. 基础示例：水分子的能量计算

### 5.1 完整代码

```cpp
/**
 * 基础示例：水分子能量计算
 * 
 * 本程序演示如何使用PySCF-C++进行基本的DFT能量计算
 */

#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

#include "molecule.h"
#include "scf.h"

using namespace pyscf;

int main() {
    std::cout << "========================================\n";
    std::cout << "  水分子 B3LYP/6-31G 能量计算\n";
    std::cout << "========================================\n\n";
    
    // ========== 第一步：创建分子对象 ==========
    // 1. 定义分子结构
    // 坐标单位是 Bohr (1 Bohr = 0.529 Å)
    // 水分子：O在原点，H在~104.5度键角
    
    auto mol = std::make_shared<Molecule>();
    
    // 添加原子：原子序数, x, y, z (单位：Bohr)
    mol->add_atom(8,  0.0,  0.0,  0.0);       // O at origin
    mol->add_atom(1,  1.4,  1.2,  0.0);      // H1
    mol->add_atom(1, -1.4,  1.2,  0.0);      // H2
    
    // 2. 设置基组
    mol->set_basis("6-31g");
    
    // 打印分子信息
    std::cout << "分子信息:\n";
    std::cout << "  原子数: " << mol->num_atoms() << "\n";
    std::cout << "  电子数: " << mol->num_electrons() << "\n";
    std::cout << "  基函数数: " << mol->num_basis_functions() << "\n";
    std::cout << "  核排斥能: " << mol->nuclear_repulsion_energy() << " Eh\n\n";
    
    // ========== 第二步：创建SCF对象 ==========
    // 创建限制性Kohn-Sham (RKS) 计算器
    
    auto scf = std::make_shared<dft::RKS>(mol);
    
    // 设置交换-相关泛函
    scf->set_xc_functional("b3lyp");
    
    // 设置SCF选项（可选）
    dft::SCFOptions opts;
    opts.conv_tol = 1e-8;      // 收敛精度
    opts.max_cycle = 100;       // 最大迭代次数
    opts.verbose = true;        // 打印详细信息
    scf->set_options(opts);
    
    // ========== 第三步：运行计算 ==========
    std::cout << "开始SCF迭代...\n\n";
    
    scf->compute();
    
    // ========== 第四步：获取结果 ==========
    std::cout << "\n========================================\n";
    std::cout << "            计算结果\n";
    std::cout << "========================================\n\n";
    
    std::cout << "总能量: " << std::scientific << std::setprecision(8)
              << scf->get_total_energy() << " Eh\n";
    std::cout << "电子能量: " << scf->get_electronic_energy() << " Eh\n";
    std::cout << "核排斥能: " << scf->get_nuclear_repulsion() << " Eh\n\n";
    
    // 获取分子轨道信息
    const auto& mo_coeffs = scf->get_mo_coefficients();
    std::cout << "分子轨道系数矩阵大小: " << mo_coeffs.size() << "\n";
    
    return 0;
}
```

### 5.2 代码详解

#### 5.2.1 包含头文件

```cpp
#include "molecule.h"   // 分子定义
#include "scf.h"       // SCF计算
```

#### 5.2.2 创建分子

```cpp
auto mol = std::make_shared<Molecule>();
mol->add_atom(8, 0.0, 0.0, 0.0);  // 原子序数8=氧
```

**原子序数参考表：**
| 元素 | 原子序数 | 元素 | 原子序数 |
|------|----------|------|----------|
| H | 1 | S | 16 |
| C | 6 | Cl | 17 |
| N | 7 | Br | 35 |
| O | 8 | I | 53 |

#### 5.2.3 设置基组

```cpp
mol->set_basis("sto-3g");    // 最简单的基组
mol->set_basis("6-31g");     // Pople基组
mol->set_basis("6-31g*");    // 带极化
mol->set_basis("6-31g**");   // 双极化
```

#### 5.2.4 创建SCF对象并计算

```cpp
auto scf = std::make_shared<dft::RKS>(mol);
scf->set_xc_functional("b3lyp");
scf->compute();
```

---

## 6. 进阶示例：几何优化

### 6.1 完整代码

```cpp
/**
 * 几何优化示例
 * 
 * 使用BFGS算法优化分子几何结构
 */

#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

#include "molecule.h"
#include "scf.h"
#include "optimizer.h"

using namespace pyscf;

int main() {
    std::cout << "========================================\n";
    std::cout << "  H2分子几何优化 (BFGS)\n";
    std::cout << "========================================\n\n";
    
    // ========== 创建分子 ==========
    auto mol = std::make_shared<Molecule>();
    
    // 初始猜测：H-H距离1.4 Bohr (实验值~1.4 Bohr)
    mol->add_atom(1, 0.0, 0.0, 0.0);      // H1 at origin
    mol->add_atom(1, 1.4, 0.0, 0.0);     // H2 at 1.4 Bohr
    mol->set_basis("sto-3g");
    
    std::cout << "初始几何构型:\n";
    std::cout << "  H1: (0.0, 0.0, 0.0)\n";
    std::cout << "  H2: (1.4, 0.0, 0.0)\n";
    std::cout << "  距离: 1.4 Bohr\n\n";
    
    // ========== 创建SCF对象 ==========
    auto scf = std::make_shared<dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    
    // ========== 创建优化器 ==========
    geom::BFGSOptimizer optimizer;
    
    // 设置收敛参数
    geom::OptConvergence conv;
    conv.gradient_max = 4.5e-4;   // 最大梯度 (Eh/Bohr)
    conv.gradient_rms = 3.0e-4;    // RMS梯度
    conv.step_max = 1.8e-3;       // 最大步长 (Bohr)
    conv.step_rms = 1.2e-3;       // RMS步长
    conv.energy_tol = 1.0e-6;     // 能量收敛
    conv.max_iterations = 50;     // 最大迭代次数
    conv.verbose = true;          // 打印详细信息
    optimizer.set_convergence(conv);
    
    // ========== 运行优化 ==========
    std::cout << "开始几何优化...\n\n";
    
    auto result = optimizer.optimize(scf, mol);
    
    // ========== 输出结果 ==========
    std::cout << "\n========================================\n";
    std::cout << "           优化结果\n";
    std::cout << "========================================\n\n";
    
    std::cout << "收敛: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "迭代次数: " << result.iterations << "\n";
    std::cout << "最终能量: " << std::scientific << std::setprecision(8)
              << result.final_energy << " Eh\n";
    
    // 计算最终H-H距离
    auto final_mol = std::make_shared<Molecule>(*mol);
    final_mol->set_coordinates(result.coordinates);
    auto atom1 = final_mol->get_atom(0);
    auto atom2 = final_mol->get_atom(1);
    double dx = atom2.x - atom1.x;
    double dy = atom2.y - atom1.y;
    double dz = atom2.z - atom1.z;
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "最终H-H距离: " << std::fixed << std::setprecision(4)
              << distance << " Bohr (" << distance * 0.529177 << " Å)\n";
    
    return 0;
}
```

### 6.2 BFGS优化器原理

BFGS是一种拟牛顿法，通过近似Hessian矩阵来加速收敛：

```
搜索方向: p_k = -H_k × g_k
步长: α_k 通过线搜索确定
更新: H_{k+1} = H_k + (y_k × y_k^T) / (y_k^T × s_k) 
                - (H_k × s_k × s_k^T × H_k) / (s_k^T × H_k × s_k)
其中: s_k = x_{k+1} - x_k
      y_k = g_{k+1} - g_k
```

### 6.3 收敛判据

| 参数 | 默认值 | 说明 |
|------|--------|------|
| gradient_max | 4.5e-4 | 最大梯度分量 (Eh/Bohr) |
| gradient_rms | 3.0e-4 | 梯度均方根 |
| step_max | 1.8e-3 | 最大位移 (Bohr) |
| step_rms | 1.2e-3 | 位移均方根 |
| energy_tol | 1.0e-6 | 能量变化 |

---

## 7. 进阶示例：Mayer键级分析

### 7.1 完整代码

```cpp
/**
 * Mayer键级分析示例
 * 
 * 计算分子中的键级和键类型
 */

#include <iostream>
#include <iomanip>
#include <memory>

#include "molecule.h"
#include "scf.h"
#include "mayer_bond.h"

using namespace pyscf;

int main() {
    std::cout << "========================================\n";
    std::cout << "  Mayer键级分析示例\n";
    std::cout << "========================================\n\n";
    
    // ========== 水分子 ==========
    std::cout << "--- 水分子 (H2O) ---\n\n";
    
    auto mol = std::make_shared<Molecule>();
    mol->add_atom(8,  0.0,  0.0,  0.0);   // O
    mol->add_atom(1,  1.4,  1.2,  0.0);   // H1
    mol->add_atom(1, -1.4,  1.2,  0.0);   // H2
    mol->set_basis("sto-3g");
    
    // SCF计算
    auto scf = std::make_shared<dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    scf->compute();
    
    std::cout << "总能量: " << scf->get_total_energy() << " Eh\n\n";
    
    // 键级分析
    auto result = bond::analyze_bonds(scf, mol);
    bond::print_bond_analysis(result, *mol);
    
    // ========== 氮气分子 (三键) ==========
    std::cout << "\n--- 氮气分子 (N2) ---\n\n";
    
    auto mol_n2 = std::make_shared<Molecule>();
    mol_n2->add_atom(7, 0.0, 0.0, 0.0);   // N
    mol_n2->add_atom(7, 2.0, 0.0, 0.0);   // N
    mol_n2->set_basis("sto-3g");
    
    auto scf_n2 = std::make_shared<dft::RKS>(mol_n2);
    scf_n2->set_xc_functional("b3lyp");
    scf_n2->compute();
    
    auto result_n2 = bond::analyze_bonds(scf_n2, mol_n2);
    bond::print_bond_analysis(result_n2, *mol_n2);
    
    return 0;
}
```

### 7.2 Mayer键级公式

Mayer键级定义为：

```
W_AB = Σ_μ∈A Σ_ν∈B P_μν × S_μν

其中：
- P_μν: 密度矩阵元
- S_μν: 重叠矩阵元
```

### 7.3 键级参考值

| 分子 | 键 | Mayer键级 | 键类型 |
|------|-----|-----------|--------|
| H2 | H-H | ~1.0 | 单键 |
| H2O | O-H | ~0.8 | 单键 |
| N2 | N≡N | ~3.0 | 三键 |
| CO | C≡O | ~2.5 | 三键 |
| C2H4 | C=C | ~2.0 | 双键 |

---

## 8. 完整示例程序

### 8.1 能量计算模板

```cpp
#include <iostream>
#include <memory>
#include "molecule.h"
#include "scf.h"

int main() {
    // 1. 创建分子
    auto mol = std::make_shared<pyscf::Molecule>();
    // mol->add_atom(原子序数, x, y, z);
    mol->set_basis("基组名称");
    
    // 2. 创建SCF计算器
    auto scf = std::make_shared<pyscf::dft::RKS>(mol);
    scf->set_xc_functional("b3lyp");
    
    // 3. 运行计算
    scf->compute();
    
    // 4. 获取结果
    std::cout << "总能量: " << scf->get_total_energy() << " Eh\n";
    
    return 0;
}
```

### 8.2 常用分子几何构型

```cpp
// 水分子 (实验值)
O  0.0000  0.0000  0.0000
H  0.7570  0.5860  0.0000
H -0.7570  0.5860  0.0000
(单位: Angstrom)

// 甲烷 (四面体)
C  0.0000  0.0000  0.0000
H  0.6290  0.6290  0.6290
H -0.6290 -0.6290  0.6290
H -0.6290  0.6290 -0.6290
H  0.6290 -0.6290 -0.6290

// 氨气 (三角锥)
N  0.0000  0.0000  0.0000
H  0.0000  0.9390  0.8110
H  0.8130 -0.4700  0.8110
H -0.8130 -0.4700  0.8110
```

### 8.3 编译链接

```bash
# 编译
g++ -std=c++17 -O3 -I./include \
    main.cpp \
    src/molecule.cpp \
    src/dft/scf.cpp \
    src/dft/b3lyp_xc.cpp \
    src/dft/lebedev_grid.cpp \
    src/gto/gto_integrals.cpp \
    -o my_program \
    -L/usr/local/lib -lcint \
    -lopenblas -llapack -lm -lpthread

# 运行
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH ./my_program
```

---

## 9. 基组处理

### 9.1 支持的基组

| 基组 | 类型 | 说明 |
|------|------|------|
| STO-3G | 最小基 | 最早也是最简单的基组 |
| 3-21G | 分裂价 | 每个内层两个高斯，每个价层一个 |
| 6-31G | 分裂价 | 更精确的内层和价层表示 |
| 6-31G* | 极化 | 添加d极化函数 |
| 6-31G** | 双极化 | 添加p极化到H |
| 6-311G | 三分裂 | 价层三分裂 |

### 9.2 基组文件格式

```text
# STO-3G 基组示例
H S
  3
  18.7311370  0.03349460
   2.8253937  0.23472695
   0.6401217  0.81375733
  1
   0.1612778  1.00000000
```

---

## 10. 常见问题与调试

### 10.1 SCF不收敛

**问题**: 计算产生NaN或SCF不收敛

**可能原因**:
1. 初始几何不合理（原子靠得太近）
2. 基组不适合该分子
3. SCF参数需要调整

**解决方案**:
```cpp
// 降低收敛要求
opts.conv_tol = 1e-6;
opts.max_cycle = 200;
scf->set_options(opts);

// 或使用更好的初始猜测
```

### 10.2 内存不足

**问题**: 基组太大导致内存溢出

**解决方案**:
- 使用更小的基组（如sto-3g代替cc-pVTZ）
- 减少积分精度
- 使用密度拟合近似

### 10.3 编译错误

**问题**: 找不到libcint

**解决方案**:
```bash
# 检查安装位置
ldconfig -p | grep cint
# 或手动指定
cmake -DCINT_LIBRARY=/path/to/libcint.so ..
```

### 10.4 能量值异常

**问题**: 计算的能量与预期差异很大

**检查项**:
- [ ] 坐标单位是否正确（Bohr）
- [ ] 基组是否正确加载
- [ ] 核电荷数是否正确
- [ ] 电子数是否正确

---

## 附录：完整API参考

### A.1 Molecule类

```cpp
class Molecule {
public:
    void add_atom(int atomic_number, double x, double y, double z);
    void set_basis(const std::string& basis_name);
    int num_atoms() const;
    int num_electrons() const;
    int num_basis_functions() const;
    double nuclear_repulsion_energy() const;
    std::vector<double> get_coordinates() const;
    void set_coordinates(const std::vector<double>& coords);
};
```

### A.2 RKS类

```cpp
class RKS {
public:
    RKS(std::shared_ptr<Molecule> mol);
    void set_xc_functional(const std::string& xc);
    SCFResult compute();
    double get_total_energy() const;
    double get_electronic_energy() const;
    double get_nuclear_repulsion() const;
    const std::vector<double>& get_density_matrix() const;
    const std::vector<double>& get_overlap_matrix() const;
    const std::vector<double>& get_mo_coefficients() const;
};
```

### A.3 BFGSOptimizer类

```cpp
class BFGSOptimizer {
public:
    OptimizationResult optimize(
        std::shared_ptr<dft::RKS> scf,
        std::shared_ptr<Molecule> mol
    );
    void set_convergence(const OptConvergence& conv);
    void set_verbose(bool v);
};
```

### A.4 键级分析函数

```cpp
BondAnalysisResult analyze_bonds(
    std::shared_ptr<dft::RKS> scf,
    std::shared_ptr<Molecule> mol
);
void print_bond_analysis(
    const BondAnalysisResult& result,
    const Molecule& mol
);
```

---

## 参考文献

1. Becke, A. D. (1993). "A new mixing of Hartree-Fock and local density-functional theories". The Journal of Chemical Physics, 98(2), 1372.
2. Lee, C., Yang, W., & Parr, R. G. (1988). "Development of the Colle-Salvetti correlation-energy formula into a functional of the electron density". Physical Review B, 37(2), 785.
3. Mayer, I. (1986). "Bond order and valence indices: a personal account". Journal of Computational Chemistry, 28(1), 204-221.
4. Nocedal, J., & Wright, S. J. (2006). Numerical Optimization. Springer.

---

**⚠️ 免责声明**: 本教程由AI (OpenHands) 编写，仅供学习参考。代码不保证：
- 科学计算的准确性
- 数值稳定性
- 生产环境安全性

如需进行严肃的科学研究，请使用经过验证的商业或开源软件。

*最后更新: 2026年8月*
