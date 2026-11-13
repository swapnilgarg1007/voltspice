#include "solver.hpp"
#include <algorithm>
#include <cmath>

namespace vs {

bool luSolve(std::vector<double>& A, std::vector<double>& b, int n) {
    double scale = 0;
    for (double v : A) scale = std::max(scale, std::fabs(v));
    const double tiny = std::max(scale, 1e-30) * 1e-14;

    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(A[i * n + k]) > std::fabs(A[piv * n + k])) piv = i;
        if (std::fabs(A[piv * n + k]) < tiny) return false;
        if (piv != k) {
            for (int j = 0; j < n; ++j) std::swap(A[k * n + j], A[piv * n + j]);
            std::swap(b[k], b[piv]);
        }
        for (int i = k + 1; i < n; ++i) {
            double f = A[i * n + k] / A[k * n + k];
            if (f == 0) continue;
            for (int j = k; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
            b[i] -= f * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[i * n + j] * b[j];
        b[i] = s / A[i * n + i];
    }
    return true;
}

}  // namespace vs
