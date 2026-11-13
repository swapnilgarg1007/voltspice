#include "sparse.hpp"
#include <algorithm>
#include "solver.hpp"
#include <cstdio>
#include <random>
using namespace vs;
int main() {
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> U(-1, 1);
    double worst = 0;
    for (int trial = 0; trial < 200; ++trial) {
        int n = 5 + trial % 60;
        std::vector<double> A(size_t(n) * n, 0.0);
        std::vector<Triplet<double>> t;
        auto put = [&](int r, int c, double v) { A[size_t(r) * n + c] += v; t.push_back({r, c, v}); };
        std::vector<int> perm(n); for (int i = 0; i < n; ++i) perm[i] = i;
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int i = 0; i < n; ++i) put(i, perm[i], 1 + std::fabs(U(rng)));   // nonsingular backbone, zero diagonal mostly
        for (int k = 0; k < 3 * n; ++k) put(rng() % n, rng() % n, U(rng) * 0.3);
        std::vector<double> b(n), x;
        for (auto& v : b) v = U(rng);
        std::vector<double> Ad = A, xd = b;
        if (!luSolve(Ad, xd, n)) continue;
        SparseLU<double> lu;
        if (!lu.factor(n, t)) { std::printf("sparse failed n=%d\n", n); return 1; }
        x = b; lu.solve(x);
        double err = 0; for (int i = 0; i < n; ++i) err = std::max(err, std::fabs(x[i] - xd[i]));
        // second RHS reusing the factorisation
        std::vector<double> b2(n), x2, xd2; for (auto& v : b2) v = U(rng);
        xd2 = b2; Ad = A; luSolve(Ad, xd2, n); x2 = b2; lu.solve(x2);
        for (int i = 0; i < n; ++i) err = std::max(err, std::fabs(x2[i] - xd2[i]));
        worst = std::max(worst, err);
    }
    std::printf("worst |sparse - dense| over 200 random systems: %.2e\n", worst);
    // complex 2x2: (1+j)x + 2y = 3 ; x - j y = 1
    using C = std::complex<double>;
    SparseLU<C> lc; lc.factor(2, {{0,0,C(1,1)},{0,1,C(2,0)},{1,0,C(1,0)},{1,1,C(0,-1)}});
    std::vector<C> b{C(3,0), C(1,0)}; lc.solve(b);
    C r0 = C(1,1)*b[0] + C(2,0)*b[1], r1 = b[0] - C(0,1)*b[1];
    std::printf("complex residual: %.2e %.2e\n", std::abs(r0 - C(3,0)), std::abs(r1 - C(1,0)));
    return worst < 1e-8 ? 0 : 1;
}
