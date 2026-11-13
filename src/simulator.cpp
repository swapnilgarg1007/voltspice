#include "simulator.hpp"
#include "solver.hpp"
#include "sparse.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>

namespace vs {
namespace {

using cplx = std::complex<double>;
constexpr double kPi = 3.14159265358979323846;
constexpr double kGmin = 1e-12;     // conductance from every node to ground
constexpr double kVt = 0.025852;    // thermal voltage at 300 K
constexpr double kRelTol = 1e-3;
constexpr double kAbsTolV = 1e-6;
constexpr double kAbsTolI = 1e-9;
constexpr double kStepLimit = 1.0;  // max node-voltage change per Newton iteration
constexpr double kMaxCells = 8e6;   // refuse to build outputs larger than this

struct Ctx {
    bool tran = false;
    double t = 0, h = 0;
    bool trap = false;       // trapezoidal vs backward Euler for this step
    double srcScale = 1.0;   // for source stepping
    bool ac = false;
    double omega = 0;
};

// Sparse system: real part and imaginary part (only AC uses `im`) plus right-hand side.
struct Sys {
    std::vector<Triplet<double>> re, im;
    std::vector<double> bre, bim;
};

struct Key {
    bool tran, trap;
    double h;
    bool operator==(const Key& o) const { return tran == o.tran && trap == o.trap && h == o.h; }
};

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

class Engine {
public:
    explicit Engine(Circuit& c) : c_(c) {
        nV_ = c_.numNodes() - 1;
        int nb = 0;
        branch_.assign(c_.elements.size(), -1);
        for (size_t i = 0; i < c_.elements.size(); ++i) {
            ElemType t = c_.elements[i].type;
            if (t == ElemType::V || t == ElemType::L) branch_[i] = nb++;
            if (t == ElemType::D || t == ElemType::M) nonlinear_ = true;
        }
        n_ = nV_ + nb;
        if (n_ == 0) throw SimError("circuit is empty");
        vPrev_.assign(c_.elements.size(), 0);
        iPrev_.assign(c_.elements.size(), 0);
    }

    int size() const { return n_; }
    void fillStats(Result& r) const {
        r.size = n_; r.nnzA = nnzA_; r.nnzLU = nnzLU_;
        r.factorizations = nFactor_; r.solves = nSolve_;
        r.factorSeconds = tFactor_; r.solveSeconds = tSolve_;
    }

    std::vector<std::string> columns() const {
        std::vector<std::string> cols;
        for (int i = 1; i <= nV_; ++i) cols.push_back("v(" + c_.nodeName(i) + ")");
        for (size_t i = 0; i < c_.elements.size(); ++i)
            if (branch_[i] >= 0) cols.push_back("i(" + c_.elements[i].name + ")");
        return cols;
    }

    // DC operating point. Falls back to source stepping if plain Newton fails.
    void operatingPoint(std::vector<double>& x) {
        Ctx ctx;
        if (solve(x, ctx)) return;
        std::fill(x.begin(), x.end(), 0.0);
        const int steps = 25;
        for (int s = 1; s <= steps; ++s) {
            ctx.srcScale = double(s) / steps;
            if (!solve(x, ctx)) throw SimError("DC operating point did not converge");
        }
    }

    void initTransient(const std::vector<double>& x) {
        for (size_t i = 0; i < c_.elements.size(); ++i) {
            const Element& e = c_.elements[i];
            if (e.type == ElemType::C) { vPrev_[i] = volt(x, e.nodes[0]) - volt(x, e.nodes[1]); iPrev_[i] = 0; }
            if (e.type == ElemType::L) { vPrev_[i] = 0; iPrev_[i] = x[nV_ + branch_[i]]; }
        }
    }

    bool step(std::vector<double>& x, double t, double h, bool trap) {
        Ctx ctx;
        ctx.tran = true; ctx.t = t; ctx.h = h; ctx.trap = trap;
        if (!solve(x, ctx)) return false;
        for (size_t i = 0; i < c_.elements.size(); ++i) {
            const Element& e = c_.elements[i];
            if (e.type != ElemType::C && e.type != ElemType::L) continue;
            double v = volt(x, e.nodes[0]) - volt(x, e.nodes[1]);
            if (e.type == ElemType::C) {
                double g = (trap ? 2 : 1) * e.value / h;
                iPrev_[i] = g * v - g * vPrev_[i] - (trap ? iPrev_[i] : 0.0);
            } else {
                iPrev_[i] = x[nV_ + branch_[i]];
            }
            vPrev_[i] = v;
        }
        return true;
    }

    // Small-signal solution at one frequency, linearised around the operating point `xop`.
    void acSolve(const std::vector<double>& xop, double freq, std::vector<cplx>& out) {
        Ctx ctx;
        ctx.ac = true;
        ctx.omega = 2 * kPi * freq;
        Sys s;
        s.bre.assign(n_, 0.0); s.bim.assign(n_, 0.0);
        stamp(s, xop, ctx);
        std::vector<Triplet<cplx>> t;
        t.reserve(s.re.size() + s.im.size());
        for (auto& e : s.re) t.push_back({e.r, e.c, cplx(e.v, 0)});
        for (auto& e : s.im) t.push_back({e.r, e.c, cplx(0, e.v)});
        SparseLU<cplx> lu;
        auto t0 = Clock::now();
        if (!lu.factor(n_, std::move(t))) throw SimError("singular AC matrix");
        tFactor_ += since(t0); ++nFactor_;
        nnzA_ = lu.nnzA(); nnzLU_ = lu.nnzLU();
        out.resize(n_);
        for (int i = 0; i < n_; ++i) out[i] = cplx(s.bre[i], s.bim[i]);
        lu.solve(out);
        ++nSolve_;
    }

private:
    Circuit& c_;
    int nV_ = 0, n_ = 0;
    bool nonlinear_ = false;
    std::vector<int> branch_;
    std::vector<double> vPrev_, iPrev_;
    std::vector<std::pair<Key, SparseLU<double>>> cache_;  // linear circuits: reuse factorisations
    size_t nnzA_ = 0, nnzLU_ = 0;
    long nFactor_ = 0, nSolve_ = 0;
    double tFactor_ = 0, tSolve_ = 0;

    static double volt(const std::vector<double>& x, int node) { return node == 0 ? 0.0 : x[node - 1]; }

    void linearSolve(Sys& s, std::vector<double>& rhs, const Ctx& ctx) {
        if (c_.denseSolver) {
            auto t0 = Clock::now();
            std::vector<double> A(size_t(n_) * n_, 0.0);
            for (auto& e : s.re) A[size_t(e.r) * n_ + e.c] += e.v;
            if (!luSolve(A, rhs, n_)) throw SimError("singular matrix: check for a loop of voltage sources/inductors or a short");
            tFactor_ += since(t0); ++nFactor_; ++nSolve_;
            return;
        }
        if (!nonlinear_) {  // matrix is constant for a given (analysis, method, step): factor once
            Key k{ctx.tran, ctx.trap, ctx.h};
            auto it = std::find_if(cache_.begin(), cache_.end(), [&](auto& p) { return p.first == k; });
            if (it == cache_.end()) {
                if (cache_.size() >= 4) cache_.clear();
                cache_.emplace_back(k, SparseLU<double>());
                it = cache_.end() - 1;
                auto t0 = Clock::now();
                if (!it->second.factor(n_, s.re)) {
                    cache_.pop_back();
                    throw SimError("singular matrix: check for a loop of voltage sources/inductors or a short");
                }
                tFactor_ += since(t0); ++nFactor_;
                nnzA_ = it->second.nnzA(); nnzLU_ = it->second.nnzLU();
            }
            auto t1 = Clock::now();
            it->second.solve(rhs);
            tSolve_ += since(t1); ++nSolve_;
            return;
        }
        SparseLU<double> lu;
        auto t0 = Clock::now();
        if (!lu.factor(n_, s.re)) throw SimError("singular matrix: check for a loop of voltage sources/inductors or a short");
        tFactor_ += since(t0); ++nFactor_;
        nnzA_ = lu.nnzA(); nnzLU_ = lu.nnzLU();
        auto t1 = Clock::now();
        lu.solve(rhs);
        tSolve_ += since(t1); ++nSolve_;
    }

    bool solve(std::vector<double>& x, const Ctx& ctx) {
        const int maxIter = nonlinear_ ? 150 : 1;
        for (int it = 0; it < maxIter; ++it) {
            Sys s;
            s.bre.assign(n_, 0.0); s.bim.assign(n_, 0.0);
            stamp(s, x, ctx);
            std::vector<double> xn = s.bre;
            linearSolve(s, xn, ctx);
            if (!nonlinear_) { x = xn; return true; }
            bool converged = true;
            for (int i = 0; i < n_; ++i) {
                double d = xn[i] - x[i];
                double tol = kRelTol * std::max(std::fabs(xn[i]), std::fabs(x[i])) + (i < nV_ ? kAbsTolV : kAbsTolI);
                if (std::fabs(d) > tol) converged = false;
                if (i < nV_ && std::fabs(d) > kStepLimit) xn[i] = x[i] + std::copysign(kStepLimit, d);
            }
            x = xn;
            if (converged) return true;
        }
        return false;
    }

    void stamp(Sys& s, const std::vector<double>& x, const Ctx& ctx) {
        auto add = [&](int r, int c, double v) { if (r >= 0 && c >= 0) s.re.push_back({r, c, v}); };
        auto addi = [&](int r, int c, double v) { if (r >= 0 && c >= 0) s.im.push_back({r, c, v}); };
        auto rhs = [&](int r, double v) { if (r >= 0) s.bre[r] += v; };
        auto conductance = [&](int na, int nb, double g) {
            add(na - 1, na - 1, g); add(nb - 1, nb - 1, g);
            add(na - 1, nb - 1, -g); add(nb - 1, na - 1, -g);
        };
        // constant current i flowing from na through the element to nb (no effect in AC analysis)
        auto current = [&](int na, int nb, double i) { if (ctx.ac) return; rhs(na - 1, -i); rhs(nb - 1, i); };

        for (int i = 0; i < nV_; ++i) add(i, i, kGmin);

        for (size_t k = 0; k < c_.elements.size(); ++k) {
            const Element& e = c_.elements[k];
            switch (e.type) {
                case ElemType::R:
                    conductance(e.nodes[0], e.nodes[1], 1.0 / e.value);
                    break;
                case ElemType::C:
                    if (ctx.ac) {
                        double b = ctx.omega * e.value;
                        addi(e.nodes[0] - 1, e.nodes[0] - 1, b); addi(e.nodes[1] - 1, e.nodes[1] - 1, b);
                        addi(e.nodes[0] - 1, e.nodes[1] - 1, -b); addi(e.nodes[1] - 1, e.nodes[0] - 1, -b);
                    } else if (ctx.tran) {  // companion model: i = g*v + ieq
                        double g = (ctx.trap ? 2 : 1) * e.value / ctx.h;
                        double ieq = -g * vPrev_[k] - (ctx.trap ? iPrev_[k] : 0.0);
                        conductance(e.nodes[0], e.nodes[1], g);
                        current(e.nodes[0], e.nodes[1], ieq);
                    }
                    break;  // open circuit in DC
                case ElemType::L: {
                    int ib = nV_ + branch_[k];
                    add(e.nodes[0] - 1, ib, 1); add(e.nodes[1] - 1, ib, -1);
                    add(ib, e.nodes[0] - 1, 1); add(ib, e.nodes[1] - 1, -1);
                    if (ctx.ac) {
                        addi(ib, ib, -ctx.omega * e.value);
                    } else if (ctx.tran) {
                        double req = (ctx.trap ? 2 : 1) * e.value / ctx.h;
                        add(ib, ib, -req);
                        rhs(ib, -req * iPrev_[k] - (ctx.trap ? vPrev_[k] : 0.0));
                    }  // DC: short circuit (v+ - v- = 0)
                    break;
                }
                case ElemType::V: {
                    int ib = nV_ + branch_[k];
                    add(e.nodes[0] - 1, ib, 1); add(e.nodes[1] - 1, ib, -1);
                    add(ib, e.nodes[0] - 1, 1); add(ib, e.nodes[1] - 1, -1);
                    if (ctx.ac) {
                        double ph = e.wf.acPhase * kPi / 180;
                        s.bre[ib] += e.wf.acMag * std::cos(ph);
                        s.bim[ib] += e.wf.acMag * std::sin(ph);
                    } else {
                        rhs(ib, ctx.srcScale * (ctx.tran ? e.wf.eval(ctx.t) : e.wf.dcValue()));
                    }
                    break;
                }
                case ElemType::I:
                    if (ctx.ac) {
                        double ph = e.wf.acPhase * kPi / 180, re = e.wf.acMag * std::cos(ph), im = e.wf.acMag * std::sin(ph);
                        if (e.nodes[0] > 0) { s.bre[e.nodes[0] - 1] -= re; s.bim[e.nodes[0] - 1] -= im; }
                        if (e.nodes[1] > 0) { s.bre[e.nodes[1] - 1] += re; s.bim[e.nodes[1] - 1] += im; }
                    } else {
                        current(e.nodes[0], e.nodes[1], ctx.srcScale * (ctx.tran ? e.wf.eval(ctx.t) : e.wf.dcValue()));
                    }
                    break;
                case ElemType::D: {
                    auto it = c_.models.find(e.model);
                    double is = 1e-14, nn = 1;
                    if (it != c_.models.end()) { is = it->second.get("is", 1e-14); nn = it->second.get("n", 1); }
                    double vd = volt(x, e.nodes[0]) - volt(x, e.nodes[1]);
                    double nvt = nn * kVt, arg = vd / nvt, id, gd;
                    if (arg > 40) {  // linear extension keeps exp() finite
                        double ex = std::exp(40.0);
                        gd = is * ex / nvt;
                        id = is * (ex - 1) + gd * (vd - 40 * nvt);
                    } else {
                        double ex = std::exp(arg);
                        id = is * (ex - 1);
                        gd = is * ex / nvt;
                    }
                    gd += kGmin;
                    conductance(e.nodes[0], e.nodes[1], gd);
                    current(e.nodes[0], e.nodes[1], id - gd * vd);
                    break;
                }
                case ElemType::M: stampMos(e, s, x, ctx.ac); break;
            }
        }
    }

    // Square-law MOSFET with channel-length modulation, linearised for Newton-Raphson / AC.
    void stampMos(const Element& e, Sys& s, const std::vector<double>& x, bool ac) {
        const Model& m = c_.models.at(e.model);
        const double p = (m.type == "pmos") ? -1.0 : 1.0;
        const double vto = m.get("vto", p > 0 ? 0.5 : -0.5);
        const double beta = m.get("kp", 2e-5) * e.w / e.l;
        const double lambda = m.get("lambda", 0.0);

        int nd = e.nodes[0], ng = e.nodes[1], ns = e.nodes[2];
        double vd = volt(x, nd), vg = volt(x, ng), vs = volt(x, ns);
        if (p * (vd - vs) < 0) { std::swap(nd, ns); std::swap(vd, vs); }  // device is symmetric

        double vgs = p * (vg - vs), vds = p * (vd - vs), vov = vgs - p * vto;
        double id = 0, gm = 0, gds = 0;
        if (vov > 0) {
            double clm = 1 + lambda * vds;
            if (vds < vov) {  // triode
                double core = vov * vds - 0.5 * vds * vds;
                id = beta * core * clm;
                gm = beta * vds * clm;
                gds = beta * (vov - vds) * clm + beta * core * lambda;
            } else {          // saturation
                id = 0.5 * beta * vov * vov * clm;
                gm = beta * vov * clm;
                gds = 0.5 * beta * vov * vov * lambda;
            }
        }
        gds += kGmin;
        const double ieq = p * id - gm * vg - gds * vd + (gm + gds) * vs;

        auto add = [&](int r, int c, double v) { if (r >= 0 && c >= 0) s.re.push_back({r, c, v}); };
        add(nd - 1, ng - 1, gm);  add(nd - 1, nd - 1, gds);  add(nd - 1, ns - 1, -(gm + gds));
        add(ns - 1, ng - 1, -gm); add(ns - 1, nd - 1, -gds); add(ns - 1, ns - 1, gm + gds);
        if (!ac) {
            if (nd > 0) s.bre[nd - 1] -= ieq;
            if (ns > 0) s.bre[ns - 1] += ieq;
        }
    }
};

std::vector<double> acFrequencies(const Analysis& a) {
    std::vector<double> f;
    if (a.acSweep == "lin") {
        if (a.acPoints == 1) f.push_back(a.fstart);
        else for (int k = 0; k < a.acPoints; ++k) f.push_back(a.fstart + k * (a.fstop - a.fstart) / (a.acPoints - 1));
    } else {
        const double base = a.acSweep == "dec" ? 10.0 : 2.0;
        const long n = static_cast<long>(std::floor(std::log(a.fstop / a.fstart) / std::log(base) * a.acPoints + 1e-9)) + 1;
        for (long k = 0; k < n; ++k) f.push_back(a.fstart * std::pow(base, double(k) / a.acPoints));
    }
    return f;
}

}  // namespace

Result simulate(Circuit circuit) {
    Engine eng(circuit);
    Result res;
    const std::vector<std::string> allCols = eng.columns();
    const int unknowns = static_cast<int>(allCols.size());

    // Which unknowns to output (.probe), default all.
    std::vector<int> keep;
    if (circuit.probes.empty()) {
        for (int i = 0; i < unknowns; ++i) keep.push_back(i);
    } else {
        for (const auto& p : circuit.probes) {
            auto it = std::find(allCols.begin(), allCols.end(), p);
            if (it == allCols.end()) throw SimError("unknown probe '" + p + "'");
            int idx = static_cast<int>(it - allCols.begin());
            if (std::find(keep.begin(), keep.end(), idx) == keep.end()) keep.push_back(idx);
        }
    }
    auto guard = [&](double rows, double perSignal) {
        if (rows * (keep.size() * perSignal + 1) > kMaxCells)
            throw SimError("output too large (" + std::to_string(static_cast<long long>(rows)) + " rows x " +
                           std::to_string(keep.size()) + " signals); select signals with .probe v(node) ...");
    };
    auto project = [&](double lead, bool hasLead, const std::vector<double>& x) {
        std::vector<double> row;
        row.reserve(keep.size() + 1);
        if (hasLead) row.push_back(lead);
        for (int k : keep) row.push_back(x[k]);
        return row;
    };
    auto header = [&](const char* lead) {
        if (lead) res.columns.push_back(lead);
        for (int k : keep) res.columns.push_back(allCols[k]);
    };

    std::vector<double> x(unknowns, 0.0);
    Analysis::Kind kind = circuit.analysis.kind == Analysis::NONE ? Analysis::OP : circuit.analysis.kind;
    const Analysis& a = circuit.analysis;

    if (kind == Analysis::OP) {
        res.analysis = "op";
        header(nullptr);
        guard(1, 1);
        eng.operatingPoint(x);
        res.rows.push_back(project(0, false, x));
    } else if (kind == Analysis::DC) {
        res.analysis = "dc";
        auto it = std::find_if(circuit.elements.begin(), circuit.elements.end(),
                               [&](const Element& e) { return e.name == a.sweepSource; });
        if (it == circuit.elements.end() || (it->type != ElemType::V && it->type != ElemType::I))
            throw SimError("sweep source '" + a.sweepSource + "' not found");
        double span = a.stop - a.start;
        if (span * a.step < 0) throw SimError(".dc step has the wrong sign");
        long n = static_cast<long>(std::floor(span / a.step + 1e-9)) + 1;
        guard(double(n), 1);
        header(a.sweepSource.c_str());
        // The engine holds a reference to `circuit`, so editing *it updates the sweep source.
        for (long k = 0; k < n; ++k) {
            double v = a.start + k * a.step;
            it->wf.kind = Waveform::DC;
            it->wf.dc = v;
            try { eng.operatingPoint(x); }
            catch (const SimError&) { std::fill(x.begin(), x.end(), 0.0); eng.operatingPoint(x); }
            res.rows.push_back(project(v, true, x));
        }
    } else if (kind == Analysis::TRAN) {
        res.analysis = "tran";
        long n = static_cast<long>(std::llround(a.tstop / a.tstep));
        if (n > 2000000) throw SimError("too many time steps; increase tstep");
        guard(double(n + 1), 1);
        header("time");
        eng.operatingPoint(x);
        eng.initTransient(x);
        res.rows.push_back(project(0.0, true, x));
        for (long k = 1; k <= n; ++k) {
            double t = k * a.tstep;
            // First step uses backward Euler to avoid trapezoidal ringing at source discontinuities.
            if (!eng.step(x, t, a.tstep, circuit.trapezoidal && k > 1))
                throw SimError("transient step did not converge at t=" + std::to_string(t));
            res.rows.push_back(project(t, true, x));
        }
    } else {
        res.analysis = "ac";
        bool anyAc = false;
        for (const auto& e : circuit.elements)
            if ((e.type == ElemType::V || e.type == ElemType::I) && e.wf.acMag != 0) anyAc = true;
        if (!anyAc) throw SimError("no AC stimulus: add 'AC 1' to a source, e.g. V1 in 0 DC 0 AC 1");
        std::vector<double> freqs = acFrequencies(a);
        if (freqs.size() > 200000) throw SimError("too many frequency points");
        guard(double(freqs.size()), 2);
        res.columns.push_back("freq");
        for (int k : keep) {  // v(out) -> vdb(out), vp(out)
            const std::string& nm = allCols[k];
            res.columns.push_back(nm.substr(0, 1) + "db" + nm.substr(1));
            res.columns.push_back(nm.substr(0, 1) + "p" + nm.substr(1));
        }
        eng.operatingPoint(x);  // bias point for nonlinear devices
        std::vector<cplx> xc;
        for (double f : freqs) {
            eng.acSolve(x, f, xc);
            std::vector<double> row{f};
            for (int k : keep) {
                double mag = std::abs(xc[k]);
                row.push_back(20 * std::log10(std::max(mag, 1e-20)));
                row.push_back(std::arg(xc[k]) * 180 / kPi);
            }
            res.rows.push_back(std::move(row));
        }
    }
    eng.fillStats(res);
    return res;
}

}  // namespace vs
