// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include "fem_solver.h"

int main() {
    // Dimensões do guia de onda
    double a = 1.0; // largura
    double b = 0.5; // altura
    int nx = 10;    // número de elementos em x
    int ny = 5;     // número de elementos em y

    // Geração da malha e das matrizes FEM
    std::vector<std::vector<double>> K, M;
    generate_fem_matrices(a, b, nx, ny, K, M);

    // Resolução do problema de autovalores: K * x = lambda * M * x
    std::vector<double> eigenvalues;
    solve_generalized_eigenproblem(K, M, eigenvalues);

    // Exibe os primeiros autovalores (k^2)
    std::cout << "Primeiros autovalores (k^2):\n";
    for (size_t i = 0; i < std::min(eigenvalues.size(), size_t(5)); ++i) {
        std::cout << "λ" << i+1 << " = " << eigenvalues[i] << "\n";
    }

    // Salvar autovalores
std::ofstream out("autovalores_fem.txt");
for (double val : eigenvalues)
    out << val << "\n";

    return 0;
}
