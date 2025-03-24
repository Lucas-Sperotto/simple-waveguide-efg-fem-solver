# Waveguide Solver: FEM vs EFGM

This project implements the solution to a **homogeneous waveguide problem** using two numerical methods:
- **FEM (Finite Element Method)**
- **EFGM (Element-Free Galerkin Method)**

We solve the **2D Helmholtz eigenvalue problem** and compare the fundamental propagation modes (k²) using both approaches.

## 🧠 Problem Description

We solve the equation:

```
∇²φ + k²φ = 0  in Ω  
φ = 0          on ∂Ω
```

Where `Ω` is a rectangular waveguide domain. This leads to a generalized eigenvalue problem:

```
K x = λ M x    where λ = k²
```

- `K`: stiffness matrix
- `M`: mass matrix

---

## 🚀 How to Run

### 🐧 Dependencies (Ubuntu / WSL)

```bash
sudo apt update
sudo apt install build-essential liblapacke-dev libblas-dev liblapack-dev
sudo apt install python3 python3-pip
pip3 install numpy matplotlib
```

### 🔧 Compilation

```bash
g++ main.cpp fem_solver.cpp -o fem_solver -llapacke -llapack -lblas
g++ main_efgm.cpp efgm_solver.cpp mls_shape.cpp fem_solver.cpp -o efgm_solver -llapacke -llapack -lblas -I /usr/include/eigen3
```

> Make sure you have the [Eigen](https://eigen.tuxfamily.org/) headers installed (`libeigen3-dev`).

### ▶️ Execution

```bash
./fem_solver
./efgm_solver
python3 plot_eigenvalues.py
```

This will generate:
- `autovalores_fem.txt`
- `autovalores_efgm.txt`
- `comparacao_autovalores.png`

---

## 📊 Output Example

![Comparison](comparacao_autovalores.png)

---

## 📚 Method Explanation

### 🔷 FEM (Finite Element Method)

- Structured rectangular mesh
- Quadrilateral bilinear elements (Q4)
- Stiffness and mass matrices assembled using analytical formulas
- Dirichlet boundary conditions enforced directly

### 🔶 EFGM (Element-Free Galerkin Method)

- Node-based cloud (no mesh)
- MLS (Moving Least Squares) shape functions
- Gauss integration over regular background cells
- Boundary conditions applied via penalty method

---

## 📁 Files Overview

| File                 | Description                         |
|----------------------|-------------------------------------|
| `main.cpp`           | FEM solver                          |
| `main_efgm.cpp`      | EFGM solver                         |
| `fem_solver.cpp/.h`  | FEM mesh generation and assembly    |
| `efgm_solver.cpp/.h` | Point cloud and EFGM matrix assembly |
| `mls_shape.cpp/.h`   | MLS shape function implementation   |
| `utils.cpp/.h`       | Support functions (e.g. weights)    |
| `plot_eigenvalues.py`| Python visualization script         |

---

## 📄 License

This project is licensed under the MIT License.
