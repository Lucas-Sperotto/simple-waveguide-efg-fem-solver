// fem_solver.h
#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <vector>

// Geração das matrizes K e M para uma malha retangular estruturada
void generate_fem_matrices(double a, double b, int nx, int ny,
                           std::vector<std::vector<double>>& K,
                           std::vector<std::vector<double>>& M);

// Resolução do problema generalizado de autovalores K x = λ M x
void solve_generalized_eigenproblem(const std::vector<std::vector<double>>& K,
                                     const std::vector<std::vector<double>>& M,
                                     std::vector<double>& eigenvalues);

#endif
