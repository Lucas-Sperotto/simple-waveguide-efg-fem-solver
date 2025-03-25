// fem_solver.cpp
#include "fem_solver.h"
#include <iostream>
#include <cmath>
#include <cstring> // memcpy
#include <lapacke.h>

using std::vector;
/*
// Matrizes locais para elemento Q4 (simplificadas)
void compute_element_matrices(double hx, double hy, double Ke[4][4], double Me[4][4]) {
    double A = hx * hy; // área do elemento

    // Matriz de rigidez local simplificada (para Q4 e malha regular)
    double valK = 1.0 / A;
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
*/
#include <array>
#include <cmath>

// Funções de forma e derivadas no espaço natural (ξ, η)
void shape_functions(double xi, double eta,
                     std::array<double, 4>& N,
                     std::array<double, 4>& dN_dxi,
                     std::array<double, 4>& dN_deta) {
    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);

    dN_dxi[0]  = -0.25 * (1 - eta);
    dN_dxi[1]  =  0.25 * (1 - eta);
    dN_dxi[2]  =  0.25 * (1 + eta);
    dN_dxi[3]  = -0.25 * (1 + eta);

    dN_deta[0] = -0.25 * (1 - xi);
    dN_deta[1] = -0.25 * (1 + xi);
    dN_deta[2] =  0.25 * (1 + xi);
    dN_deta[3] =  0.25 * (1 - xi);
}

void compute_element_matrices(double hx, double hy, double Ke[4][4], double Me[4][4]) {
    const double gauss_pts[2] = { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };
    const double weights[2] = { 1.0, 1.0 };

    std::array<double, 4> N, dN_dxi, dN_deta;

    // Inicializa Ke e Me com zero
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            Ke[i][j] = Me[i][j] = 0.0;

    // Jacobiano constante para elemento retangular alinhado aos eixos
    double J[2][2] = {
        { hx / 2.0, 0 },
        { 0, hy / 2.0 }
    };
    double detJ = (hx * hy) / 4.0;

    // Inverso do Jacobiano
    double invJ[2][2] = {
        { 2.0 / hx, 0 },
        { 0, 2.0 / hy }
    };

    for (int i = 0; i < 2; ++i) {
        double xi = gauss_pts[i];
        double wi = weights[i];
        for (int j = 0; j < 2; ++j) {
            double eta = gauss_pts[j];
            double wj = weights[j];

            // Avalia funções de forma e derivadas
            shape_functions(xi, eta, N, dN_dxi, dN_deta);

            // Derivadas reais (em x, y) via J⁻¹
            std::array<std::array<double, 2>, 4> dNdx;
            for (int a = 0; a < 4; ++a) {
                dNdx[a][0] = invJ[0][0] * dN_dxi[a] + invJ[0][1] * dN_deta[a]; // ∂N/∂x
                dNdx[a][1] = invJ[1][0] * dN_dxi[a] + invJ[1][1] * dN_deta[a]; // ∂N/∂y
            }

            // Peso total do ponto de integração
            double weight = wi * wj * detJ;

            for (int a = 0; a < 4; ++a) {
                for (int b = 0; b < 4; ++b) {
                    double grad_a_dot_grad_b =
                        dNdx[a][0] * dNdx[b][0] + dNdx[a][1] * dNdx[b][1];

                    Ke[a][b] += grad_a_dot_grad_b * weight;
                    Me[a][b] += N[a] * N[b] * weight;
                }
            }
        }
    }
}


/*
void compute_element_matrices(double hx, double hy, double Ke[4][4], double Me[4][4]) {
    double hx2 = hx * hx;
    double hy2 = hy * hy;
    double hxhy = hx * hy;

    // Matriz de Rigidez (K_e)
    double factorK = 1.0 / (6.0 * hxhy);

    double Ktemplate[4][4] = {
        { 2 * (hx2 + hy2),   -2 * hx2 + hy2,   -hx2 - hy2,      hx2 - 2 * hy2 },
        { -2 * hx2 + hy2,    2 * (hx2 + hy2),   hx2 - 2 * hy2,  -hx2 - hy2    },
        { -hx2 - hy2,         hx2 - 2 * hy2,    2 * (hx2 + hy2), -2 * hx2 + hy2 },
        { hx2 - 2 * hy2,     -hx2 - hy2,       -2 * hx2 + hy2,   2 * (hx2 + hy2) }
    };

    // Matriz de Massa (M_e)
    double factorM = hxhy / 36.0;

    double Mtemplate[4][4] = {
        {4, 2, 1, 2},
        {2, 4, 2, 1},
        {1, 2, 4, 2},
        {2, 1, 2, 4}
    };

    // Preenche Ke e Me com os valores multiplicados pelos fatores
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            Ke[i][j] = factorK * Ktemplate[i][j];
            Me[i][j] = factorM * Mtemplate[i][j];
        }
}
*/
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
                M[idx][idx] = 1.0;
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

double traceK = 0.0, traceM = 0.0;
for (int i = 0; i < n; ++i) {
    traceK += K[i * n + i];
    traceM += M[i * n + i];
}
std::cout << "Trace K: " << traceK << ", Trace M: " << traceM << "\n";


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
