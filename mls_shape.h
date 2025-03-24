#ifndef MLS_SHAPE_H
#define MLS_SHAPE_H

#include <vector>
#include "efgm_solver.h"

// Função de forma MLS e suas derivadas
struct ShapeFunction {
    std::vector<double> phi;     // valores das funções de forma
    std::vector<double> dphidx;  // derivadas em x
    std::vector<double> dphidy;  // derivadas em y
};

// Calcula função de forma no ponto (x, y)
ShapeFunction compute_mls_shape(double x, double y,
                                const std::vector<Node>& cloud,
                                double dmax); // raio de influência

#endif
