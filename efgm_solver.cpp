#include "efgm_solver.h"
#include "mls_shape.h"
#include "utils.h"
#include <cmath>
#include <lapacke.h>
#include <iostream>

void generate_point_cloud(double a, double b, int nx, int ny,
                          std::vector<Node>& cloud) {
    double hx = a / (nx - 1);
    double hy = b / (ny - 1);
    int id = 0;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            Node n;
            n.x = i * hx;
            n.y = j * hy;
            n.id = id++;
            n.is_boundary = (i == 0 || j == 0 || i == nx-1 || j == ny-1);
            cloud.push_back(n);
        }
    }
}

void assemble_efgm_matrices(const std::vector<Node>& cloud,
    double a, double b,
    std::vector<std::vector<double>>& K,
    std::vector<std::vector<double>>& M) {
int nx = std::round(std::sqrt(cloud.size() * a / b));
int ny = std::round(cloud.size() / nx);
int ndof = cloud.size();
K.assign(ndof, std::vector<double>(ndof, 0.0));
M.assign(ndof, std::vector<double>(ndof, 0.0));

double dx = a / (nx - 1);
double dy = b / (ny - 1);
double dmax = 2.1 * std::max(dx, dy); // raio de influência MLS

// Pontos de Gauss (2x2)
const double gp[2] = {-1.0 / std::sqrt(3), 1.0 / std::sqrt(3)};
const double w[2] = {1.0, 1.0};

for (int j = 0; j < ny - 1; ++j) {
for (int i = 0; i < nx - 1; ++i) {
double x0 = i * dx;
double y0 = j * dy;

for (int gy = 0; gy < 2; ++gy) {
for (int gx = 0; gx < 2; ++gx) {
double xi = gp[gx];
double eta = gp[gy];
double wx = w[gx];
double wy = w[gy];

double x = x0 + (1 + xi) * dx / 2.0;
double y = y0 + (1 + eta) * dy / 2.0;
double detJ = (dx / 2.0) * (dy / 2.0);
double weight = wx * wy * detJ;

ShapeFunction sf = compute_mls_shape(x, y, cloud, dmax);

for (int a = 0; a < ndof; ++a) {
if (sf.phi[a] == 0.0 && sf.dphidx[a] == 0.0) continue;
for (int b = 0; b < ndof; ++b) {
    if (sf.phi[b] == 0.0 && sf.dphidx[b] == 0.0) continue;

    double k_ij = (sf.dphidx[a] * sf.dphidx[b] +
                   sf.dphidy[a] * sf.dphidy[b]) * weight;
    double m_ij = sf.phi[a] * sf.phi[b] * weight;

    K[a][b] += k_ij;
    M[a][b] += m_ij;
}
}
}
}
}
}

// Penalização das condições de contorno
double penalty = 1e10;
for (int i = 0; i < ndof; ++i) {
if (cloud[i].is_boundary) {
K[i][i] += penalty;
M[i][i] += penalty;
}
}

std::cout << "Matrizes EFGM montadas.\n";
}