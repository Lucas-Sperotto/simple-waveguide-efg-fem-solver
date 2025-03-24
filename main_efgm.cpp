#include <iostream>
#include <fstream>
#include <vector>
#include "efgm_solver.h"

int main() {
    // Dimensões do guia de onda
    double a = 1.0;
    double b = 0.5;
    int nx = 10; // pontos em x
    int ny = 5;  // pontos em y

    // Gerar nuvem de pontos
    std::vector<Node> cloud;
    generate_point_cloud(a, b, nx, ny, cloud);

    // Construir e resolver sistema EFGM
    std::vector<std::vector<double>> K, M;
    assemble_efgm_matrices(cloud, a, b, K, M);

    // Resolver problema de autovalores
    std::vector<double> eigenvalues;
    solve_generalized_eigenproblem(K, M, eigenvalues);

    std::cout << "Autovalores (k²) via EFGM:\n";
    for (size_t i = 0; i < std::min(eigenvalues.size(), size_t(5)); ++i)
        std::cout << "λ" << i+1 << " = " << eigenvalues[i] << "\n";


    // Salvar autovalores
std::ofstream out("autovalores_fem.txt");
for (double val : eigenvalues)
    out << val << "\n";    
    return 0;
}
