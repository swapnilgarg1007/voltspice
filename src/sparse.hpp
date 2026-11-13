#pragma once
// Sparse LU for MNA matrices (real or complex), header-only.
//
// Right-looking Gaussian elimination on row-sparse storage with Markowitz-style pivot
// selection: at each step the active column with the fewest entries is chosen, and within it
// the row with the fewest entries among those passing a threshold-pivoting test
// (|a| >= 0.1 * column max). This keeps fill-in low on circuit matrices while staying stable
// for MNA systems, which have structural zeros on the diagonal (voltage sources, inductors).
//
// The factorisation stores the elimination multipliers (L) and pivot rows (U), so one
// factorisation can be reused for many right-hand sides: linear transient analysis and
// DC sweeps of linear circuits factor once and then only run solve().
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

namespace vs {

template <class T>
struct Triplet {
    int r, c;
    T v;
};

template <class T>
class SparseLU {
public:
    // Factorises the n x n matrix given as triplets (duplicate positions are summed).
    // Returns false if the matrix is numerically singular.
    bool factor(int n, std::vector<Triplet<T>> trip) {
        n_ = n;
        p_.clear(); q_.clear(); urow_.clear(); l_.clear();
        std::sort(trip.begin(), trip.end(), [](const Triplet<T>& a, const Triplet<T>& b) {
            return a.r < b.r || (a.r == b.r && a.c < b.c);
        });
        std::vector<Row> rows(n);
        double scale = 0;
        for (size_t k = 0; k < trip.size();) {
            size_t j = k;
            T s{};
            while (j < trip.size() && trip[j].r == trip[k].r && trip[j].c == trip[k].c) s += trip[j++].v;
            rows[trip[k].r].emplace_back(trip[k].c, s);
            scale = std::max(scale, std::abs(s));
            k = j;
        }
        nnzA_ = 0;
        const double tiny = std::max(scale, 1e-30) * 1e-14;

        std::vector<int> colCount(n, 0);
        std::vector<std::vector<int>> colRows(n);
        for (int i = 0; i < n; ++i)
            for (auto& e : rows[i]) { ++colCount[e.first]; colRows[e.first].push_back(i); ++nnzA_; }
        std::vector<char> rowActive(n, 1);
        std::set<std::pair<int, int>> byCount;  // (active entries, column)
        for (int j = 0; j < n; ++j) byCount.insert({colCount[j], j});
        auto bump = [&](int j, int delta) {
            byCount.erase({colCount[j], j});
            colCount[j] += delta;
            byCount.insert({colCount[j], j});
        };

        p_.reserve(n); q_.reserve(n); urow_.reserve(n); l_.reserve(n);
        nnzLU_ = 0;
        for (int step = 0; step < n; ++step) {
            const int q = byCount.begin()->second;
            byCount.erase(byCount.begin());

            std::vector<int> cand;
            double cmax = 0;
            for (int i : colRows[q])
                if (rowActive[i]) { cand.push_back(i); cmax = std::max(cmax, std::abs(get(rows[i], q))); }
            if (cand.empty() || cmax < tiny) return false;

            int p = -1;
            size_t bestSize = static_cast<size_t>(-1);
            double bestAbs = 0;
            for (int i : cand) {
                double a = std::abs(get(rows[i], q));
                if (a < 0.1 * cmax) continue;
                size_t sz = rows[i].size();
                if (sz < bestSize || (sz == bestSize && a > bestAbs)) { bestSize = sz; bestAbs = a; p = i; }
            }

            Row prow = std::move(rows[p]);
            rowActive[p] = 0;
            const T piv = get(prow, q);
            std::vector<std::pair<int, T>> mult;
            for (int i : cand) {
                if (i == p) continue;
                const T f = get(rows[i], q) / piv;
                mult.emplace_back(i, f);
                Row out;
                out.reserve(rows[i].size() + prow.size());
                auto a = rows[i].begin(), ae = rows[i].end();
                auto b = prow.begin(), be = prow.end();
                while (a != ae || b != be) {
                    if (a != ae && a->first == q) { ++a; continue; }
                    if (b != be && b->first == q) { ++b; continue; }
                    if (b == be || (a != ae && a->first < b->first)) { out.push_back(*a); ++a; }
                    else if (a == ae || b->first < a->first) {  // fill-in
                        out.emplace_back(b->first, -f * b->second);
                        bump(b->first, +1);
                        colRows[b->first].push_back(i);
                        ++b;
                    } else { out.emplace_back(a->first, a->second - f * b->second); ++a; ++b; }
                }
                rows[i] = std::move(out);
            }
            for (auto& e : prow) if (e.first != q) bump(e.first, -1);  // pivot row leaves the active set

            nnzLU_ += prow.size() + mult.size();
            p_.push_back(p); q_.push_back(q);
            urow_.push_back(std::move(prow));
            l_.push_back(std::move(mult));
        }
        return true;
    }

    // Overwrites b with the solution of A x = b.
    void solve(std::vector<T>& b) const {
        for (int k = 0; k < n_; ++k) {
            const T bp = b[p_[k]];
            for (auto& m : l_[k]) b[m.first] -= m.second * bp;
        }
        std::vector<T> x(n_);
        for (int k = n_ - 1; k >= 0; --k) {
            T s = b[p_[k]], piv{};
            for (auto& e : urow_[k]) {
                if (e.first == q_[k]) piv = e.second;
                else s -= e.second * x[e.first];
            }
            x[q_[k]] = s / piv;
        }
        b = std::move(x);
    }

    size_t nnzA() const { return nnzA_; }    // non-zeros in the input matrix
    size_t nnzLU() const { return nnzLU_; }  // non-zeros stored in L and U (fill-in = nnzLU - nnzA)

private:
    using Row = std::vector<std::pair<int, T>>;
    static T get(const Row& r, int col) {
        auto it = std::lower_bound(r.begin(), r.end(), col, [](const std::pair<int, T>& e, int c) { return e.first < c; });
        return (it != r.end() && it->first == col) ? it->second : T{};
    }
    int n_ = 0;
    size_t nnzA_ = 0, nnzLU_ = 0;
    std::vector<int> p_, q_;
    std::vector<Row> urow_;
    std::vector<std::vector<std::pair<int, T>>> l_;
};

}  // namespace vs
