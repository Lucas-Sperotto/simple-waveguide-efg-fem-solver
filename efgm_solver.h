#ifndef EFGM_SOLVER_H
#define EFGM_SOLVER_H

#include <vector>

// Estrutura de ponto
struct Node {
    double x, y;
    int id;
    bool is_boundary;
};

// Geração da nuvem de pontos
void generate_point_cloud(double a, double b, int nx, int ny,
                          std::vector<Node>& cloud);

// Montagem das matrizes via EFGM
void assemble_efgm_matrices(const std::vector<Node>& cloud,
                            double a, double b,
                            std::vector<std::vector<double>>& K,
                            std::vector<std::vector<double>>& M);

// Resolução via LAPACK
void solve_generalized_eigenproblem(const std::vector<std::vector<double>>& K,
                                     const std::vector<std::vector<double>>& M,
                                     std::vector<double>& eigenvalues);

#endif
