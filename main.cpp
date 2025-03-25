// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "fem_solver.h"

int main() {
    // Dimensões do guia de onda
    double a = 1.0; // largura
    double b = 0.5; // altura
    int nx = 81;    // número de elementos em x
    int ny = 41;     // número de elementos em y

    // Geração da malha e das matrizes FEM
    std::vector<std::vector<double>> K, M;
    generate_fem_matrices(a, b, nx, ny, K, M);

    // Resolução do problema de autovalores: K * x = lambda * M * x
    std::vector<double> eigenvalues;
    solve_generalized_eigenproblem(K, M, eigenvalues);
    std::vector<double> filtered;
    for (double val : eigenvalues) {
        if (val > 1.1) // ignora autovalores triviais = 1
            filtered.push_back(val);
    }
    
    // Exibir os primeiros valores relevantes
    std::cout << "Primeiros autovalores (k^2):\n";
    for (size_t i = 0; i < std::min(filtered.size(), size_t(10)); ++i) {
        std::cout << "λ[" << i + 1 << "] = " << filtered[i] << std::endl;
    }
    // Exibe os primeiros autovalores (k^2)
//    std::cout << "Primeiros autovalores (k^2):\n";
//    for (size_t i = 0; i < std::min(eigenvalues.size(), size_t(5)); ++i) {
 //       std::cout << "λ" << i + 1 << " = " << eigenvalues[i] << "\n";
 //   }

    // Salvar autovalores
std::ofstream out("autovalores_fem.txt");
for (double val : eigenvalues)
    out << (sqrt(val) / 2.0)  << "\n";

    return 0;
}
