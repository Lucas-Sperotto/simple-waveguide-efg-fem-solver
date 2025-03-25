#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "efgm_solver.h"

int main() {
    // Dimensões do guia de onda
    double a = 2.0;
    double b = 1.0;
    int nx = 80; // pontos em x
    int ny = 40;  // pontos em y

    // Gerar nuvem de pontos
    std::vector<Node> cloud;
    generate_point_cloud(a, b, nx, ny, cloud);

    // Construir e resolver sistema EFGM
    std::vector<std::vector<double>> K, M;
    assemble_efgm_matrices(cloud, a, b, K, M);

    // Resolver problema de autovalores
    std::vector<double> eigenvalues;
    solve_generalized_eigenproblem(K, M, eigenvalues);
    std::vector<double> filtered;
    for (double val : eigenvalues) {
        if (val > 1.00001) // ignora autovalores triviais = 1
            filtered.push_back(val);
    }
    
    // Exibir os primeiros valores relevantes
    std::cout << "Primeiros autovalores (k^2):\n";
    for (size_t i = 0; i < std::min(filtered.size(), size_t(10)); ++i) {
        std::cout << "λ[" << i + 1 << "] = " << filtered[i] << "\t";
        std::cout << "Kc[" << i + 1 << "] = " << sqrt(filtered[i]) << std::endl;
    }
   // std::cout << "Autovalores (k²) via EFGM:\n";
  //  for (size_t i = 0; i < std::min(eigenvalues.size(), size_t(5)); ++i)
   //     std::cout << "λ" << i+1 << " = " << eigenvalues[i] << "\n";


    // Salvar autovalores
std::ofstream out("autovalores_efgm.txt");
for (double val : eigenvalues)
    out << val << "\n";    
    return 0;
}
