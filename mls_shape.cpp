#include "mls_shape.h"
#include "utils.h"
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

ShapeFunction compute_mls_shape(double x, double y,
                                const std::vector<Node>& cloud,
                                double dmax) {
    const int m = 3; // base linear: {1, x, y}
    int np = cloud.size();

    std::vector<int> neighbors;
    for (int i = 0; i < np; ++i) {
        double dx = cloud[i].x - x;
        double dy = cloud[i].y - y;
        double dist = std::sqrt(dx * dx + dy * dy);
        if (dist < dmax)
            neighbors.push_back(i);
    }

    int nnb = neighbors.size();
    ShapeFunction result;
    result.phi.resize(np, 0.0);
    result.dphidx.resize(np, 0.0);
    result.dphidy.resize(np, 0.0);

    if (nnb < m) return result; // MLS requer pelo menos m pontos

    MatrixXd P(m, nnb);
    MatrixXd W = MatrixXd::Zero(nnb, nnb);
    MatrixXd A(m, m);
    VectorXd px(m), dpxdx(m), dpxdy(m);

    for (int j = 0; j < nnb; ++j) {
        int id = neighbors[j];
        double xi = cloud[id].x;
        double yi = cloud[id].y;

        double dx = x - xi;
        double dy = y - yi;
        double w = weight_gaussian(dx, dy, dmax);
        W(j, j) = w;

        P(0, j) = 1.0;
        P(1, j) = xi;
        P(2, j) = yi;
    }

    px << 1.0, x, y;
    dpxdx << 0.0, 1.0, 0.0;
    dpxdy << 0.0, 0.0, 1.0;

    A = P * W * P.transpose();
    if (A.determinant() < 1e-10) return result;

    MatrixXd Ainv = A.inverse();
    MatrixXd coeff = Ainv * P * W;

    for (int j = 0; j < nnb; ++j) {
        int id = neighbors[j];
        VectorXd pj(m);
        pj << 1.0, cloud[id].x, cloud[id].y;

        result.phi[id]    = (px.transpose()    * coeff.col(j))(0);
        result.dphidx[id] = (dpxdx.transpose() * coeff.col(j))(0);
        result.dphidy[id] = (dpxdy.transpose() * coeff.col(j))(0);
    }

    return result;
}
