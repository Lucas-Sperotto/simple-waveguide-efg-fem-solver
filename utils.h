#ifndef UTILS_H
#define UTILS_H

#include <cmath>

// Função peso Gaussiana compacta
inline double weight_gaussian(double dx, double dy, double dmax) {
    double r = std::sqrt(dx*dx + dy*dy) / dmax;
    return (r >= 1.0) ? 0.0 : std::exp(-6.25 * r * r);
}

#endif
