// fem_solver.cpp
#include "fem_solver.h"
#include <iostream>
#include <cmath>
#include <cstring> // memcpy
#include <lapacke.h>

using std::vector;

// Matrizes locais para elemento Q4 (simplificadas)
void compute_element_matrices(double hx, double hy, double Ke[4][4], double Me[4][4]) {
    double A = hx * hy; // área do elemento

    // Matriz de rigidez local simplificada (para Q4 e malha regular)
    double valK = 1.0 / (hx * hy);
    double Ke_template[4][4] = {
        { 2, -1, -1,  0},
        {-1,  2,  0, -1},
        {-1,  0,  2, -1},
        { 0, -1, -1,  2}
    };

    // Matriz de massa local (consistente)
    double valM = A / 36.0;
    double Me_template[4][4] = {
        {4, 2, 2, 1},
        {2, 4, 1, 2},
        {2, 1, 4, 2},
        {1, 2, 2, 4}
    };

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            Ke[i][j] = valK * Ke_template[i][j];
            Me[i][j] = valM * Me_template[i][j];
        }
}

void generate_fem_matrices(double a, double b, int nx, int ny,
                           vector<vector<double>>& K, vector<vector<double>>& M) {
    int npx = nx + 1;
    int npy = ny + 1;
    int ndof = npx * npy;

    K.assign(ndof, vector<double>(ndof, 0.0));
    M.assign(ndof, vector<double>(ndof, 0.0));

    double hx = a / nx;
    double hy = b / ny;

    double Ke[4][4], Me[4][4];
    compute_element_matrices(hx, hy, Ke, Me);

    // Mapeamento local -> global
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int n0 = j * npx + i;
            int n1 = n0 + 1;
            int n2 = n0 + npx;
            int n3 = n2 + 1;

            int nodes[4] = {n0, n1, n2, n3};

            for (int m = 0; m < 4; ++m)
                for (int n = 0; n < 4; ++n) {
                    K[nodes[m]][nodes[n]] += Ke[m][n];
                    M[nodes[m]][nodes[n]] += Me[m][n];
                }
        }
    }

    // Aplicar condições de contorno de Dirichlet homogêneas (φ = 0)
    for (int j = 0; j < npy; ++j) {
        for (int i = 0; i < npx; ++i) {
            if (i == 0 || i == npx - 1 || j == 0 || j == npy - 1) {
                int idx = j * npx + i;
                for (int k = 0; k < ndof; ++k) {
                    K[idx][k] = K[k][idx] = 0.0;
                    M[idx][k] = M[k][idx] = 0.0;
                }
                K[idx][idx] = 1.0;
            }
        }
    }
}

void solve_generalized_eigenproblem(const vector<vector<double>>& K_mat,
    const vector<vector<double>>& M_mat,
    vector<double>& eigenvalues) {
int n = K_mat.size();
if (n == 0 || K_mat[0].size() != n) {
std::cerr << "Matrizes inválidas!\n";
return;
}

// Alocar e copiar matrizes para formato compatível com LAPACK
double* K = new double[n * n];
double* M = new double[n * n];

for (int i = 0; i < n; ++i)
for (int j = 0; j < n; ++j) {
K[i * n + j] = K_mat[i][j];
M[i * n + j] = M_mat[i][j];
}

eigenvalues.resize(n);

int info = LAPACKE_dsygv(LAPACK_ROW_MAJOR, 1, 'V', 'U', n,
K, n, M, n, eigenvalues.data());

if (info > 0) {
std::cerr << "A matriz não pôde ser diagonalizada.\n";
} else if (info < 0) {
std::cerr << "Erro no argumento " << -info << " da chamada LAPACKE_dsygv.\n";
}

delete[] K;
delete[] M;
}
