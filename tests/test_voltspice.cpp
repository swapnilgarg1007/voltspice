#include "parser.hpp"
#include "simulator.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace vs;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { ++failures; std::printf("FAIL: %s (line %d)\n", msg, __LINE__); } } while (0)

static int col(const Result& r, const std::string& name) {
    for (size_t i = 0; i < r.columns.size(); ++i) if (r.columns[i] == name) return static_cast<int>(i);
    std::printf("missing column %s\n", name.c_str());
    ++failures;
    return 0;
}
static Result run(const char* net) { return simulate(parseNetlist(net)); }
static bool near(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

static void testParseValues() {
    CHECK(near(parseValue("1k"), 1e3, 1e-9), "1k");
    CHECK(near(parseValue("2.2u"), 2.2e-6, 1e-15), "2.2u");
    CHECK(near(parseValue("1meg"), 1e6, 1e-6), "1meg");
    CHECK(near(parseValue("10n"), 1e-8, 1e-18), "10n");
    CHECK(near(parseValue("4.7kohm"), 4700, 1e-9), "4.7kohm");
    bool threw = false;
    try { parseValue("abc"); } catch (const ParseError&) { threw = true; }
    CHECK(threw, "bad number throws");
}

static void testDivider() {
    Result r = run("V1 in 0 DC 10\nR1 in out 1k\nR2 out 0 3k\n.op\n");
    CHECK(near(r.rows[0][col(r, "v(out)")], 7.5, 1e-6), "divider v(out)=7.5");
    CHECK(near(r.rows[0][col(r, "i(v1)")], -2.5e-3, 1e-9), "source current -2.5mA (SPICE sign)");
}

static void testRcStep(bool trap) {
    // tau = 1 ms; v(t) = 5 (1 - exp(-t/tau))
    Result r = run(trap ? "V1 in 0 PULSE(0 5 0 0 0 1 2)\nR1 in out 1k\nC1 out 0 1u\n.tran 5u 5m\n"
                        : "V1 in 0 PULSE(0 5 0 0 0 1 2)\nR1 in out 1k\nC1 out 0 1u\n.options method=be\n.tran 5u 5m\n");
    int t = col(r, "time"), v = col(r, "v(out)");
    double maxErr = 0;
    for (auto& row : r.rows) maxErr = std::max(maxErr, std::fabs(row[v] - 5 * (1 - std::exp(-row[t] / 1e-3))));
    std::printf("  RC step max error (%s): %.2e V\n", trap ? "trap" : "be", maxErr);
    CHECK(maxErr < (trap ? 5e-3 : 2e-2), "RC step response matches analytic");
}

static void testRlStep() {
    // tau = L/R = 1 ms; i(t) = (V/R)(1 - exp(-t/tau)), V=1, R=1
    Result r = run("V1 in 0 PULSE(0 1 0 0 0 1 2)\nR1 in mid 1\nL1 mid 0 1m\n.tran 5u 5m\n");
    int t = col(r, "time"), i = col(r, "i(l1)");
    double maxErr = 0;
    for (auto& row : r.rows) maxErr = std::max(maxErr, std::fabs(row[i] - (1 - std::exp(-row[t] / 1e-3))));
    std::printf("  RL step max error: %.2e A\n", maxErr);
    CHECK(maxErr < 5e-3, "RL step response matches analytic");
}

static void testDiode() {
    Result r = run("V1 in 0 DC 5\nR1 in d 1k\nD1 d 0 dm\n.model dm D(IS=1e-14 N=1)\n.op\n");
    double vd = r.rows[0][col(r, "v(d)")];
    double iR = (5 - vd) / 1e3, iD = 1e-14 * (std::exp(vd / 0.025852) - 1);
    std::printf("  diode vd = %.4f V\n", vd);
    CHECK(vd > 0.6 && vd < 0.8, "diode forward drop plausible");
    CHECK(std::fabs(iR - iD) / iR < 1e-3, "diode satisfies KCL");
}

static void testCmosInverter() {
    // Symmetric inverter: beta_n = beta_p, so Vm = VDD/2 (matches cmos-inverter-explorer's formula).
    Result r = run(
        "VDD vdd 0 DC 1.8\nVIN in 0 DC 0\n"
        "M1 out in 0 0 nm W=2u L=1u\nM2 out in vdd vdd pm W=4u L=1u\n"
        ".model nm NMOS(VTO=0.4 KP=100u LAMBDA=0.05)\n.model pm PMOS(VTO=-0.4 KP=50u LAMBDA=0.05)\n"
        ".dc VIN 0 1.8 0.01\n");
    int vin = col(r, "vin"), vo = col(r, "v(out)");
    (void)vin;
    CHECK(near(r.rows.front()[vo], 1.8, 0.01), "Vout high when Vin low");
    CHECK(near(r.rows.back()[vo], 0.0, 0.01), "Vout low when Vin high");
    bool mono = true;
    double vm = -1;
    for (size_t k = 1; k < r.rows.size(); ++k) {
        if (r.rows[k][vo] > r.rows[k - 1][vo] + 1e-6) mono = false;
        double a = r.rows[k - 1][vo] - r.rows[k - 1][0], b = r.rows[k][vo] - r.rows[k][0];
        if (vm < 0 && a > 0 && b <= 0) vm = r.rows[k][0];
    }
    std::printf("  inverter Vm = %.3f V (expected 0.900)\n", vm);
    CHECK(mono, "VTC monotonic decreasing");
    CHECK(near(vm, 0.9, 0.02), "switching threshold near VDD/2");
}

static void testRectifier() {
    // Half-wave rectifier with filter cap: output should sit near peak minus one diode drop.
    Result r = run("V1 in 0 SIN(0 5 50)\nD1 in out dm\nC1 out 0 100u\nR1 out 0 10k\n"
                   ".model dm D(IS=1e-14)\n.tran 100u 80m\n");
    int v = col(r, "v(out)");
    double vmax = 0;
    for (auto& row : r.rows) vmax = std::max(vmax, row[v]);
    std::printf("  rectifier peak = %.3f V\n", vmax);
    CHECK(vmax > 4.0 && vmax < 4.6, "rectifier peak about 5 V minus diode drop");
}


static int nearestRow(const Result& r, double f) {
    int best = 0;
    for (size_t k = 0; k < r.rows.size(); ++k)
        if (std::fabs(std::log(r.rows[k][0] / f)) < std::fabs(std::log(r.rows[best][0] / f))) best = static_cast<int>(k);
    return best;
}

static void testAcRc() {
    // R = 1k, C = 159.1549 nF -> fc = 1 kHz: -3.01 dB / -45 deg at fc, -20.04 dB / -84.3 deg at 10 fc
    Result r = run("V1 in 0 DC 0 AC 1\nR1 in out 1k\nC1 out 0 159.1549n\n.ac dec 50 10 100k\n");
    int db = col(r, "vdb(out)"), ph = col(r, "vp(out)");
    int k1 = nearestRow(r, 1000), k2 = nearestRow(r, 10000);
    std::printf("  AC RC: fc=%.1f Hz  %.3f dB %.2f deg | 10fc: %.3f dB %.2f deg\n", r.rows[k1][0], r.rows[k1][db],
                r.rows[k1][ph], r.rows[k2][db], r.rows[k2][ph]);
    CHECK(near(r.rows[k1][db], -3.0103, 0.01), "RC gain at fc is -3.01 dB");
    CHECK(near(r.rows[k1][ph], -45.0, 0.1), "RC phase at fc is -45 deg");
    CHECK(near(r.rows[k2][db], -20.043, 0.02), "RC gain at 10 fc");
    CHECK(near(r.rows[k2][ph], -84.29, 0.1), "RC phase at 10 fc");
}

static void testAcRlc() {
    // Series RLC, output across C. Q = sqrt(L/C)/R = 3.162; peak gain = Q/sqrt(1-1/(4Q^2)) = 3.202 = 10.11 dB at f0*sqrt(1-1/(2Q^2)) = 4906 Hz
    Result r = run("V1 in 0 DC 0 AC 1\nR1 in a 10\nL1 a out 1m\nC1 out 0 1u\n.ac dec 200 100 100k\n");
    int db = col(r, "vdb(out)");
    double peak = -1e9, fpk = 0;
    for (auto& row : r.rows) if (row[db] > peak) { peak = row[db]; fpk = row[0]; }
    std::printf("  AC RLC peak: %.2f dB at %.0f Hz (closed form: 10.11 dB at 4906 Hz)\n", peak, fpk);
    CHECK(near(peak, 10.105, 0.05), "RLC resonant peak matches closed form");
    CHECK(near(fpk, 4906, 60), "RLC resonant frequency matches closed form");
}

static void testAcMosAmp() {
    // Common-source amplifier: gm = beta*Vov = 1 mS, gain = -gm*RD = -8 (18.06 dB, 180 deg). Needs the bias point.
    Result r = run("VDD vdd 0 DC 3\nVG g 0 DC 1 AC 1\nRD vdd d 8k\nM1 d g 0 0 nm W=10u L=1u\n"
                   ".model nm NMOS(VTO=0.5 KP=200u LAMBDA=0)\n.ac lin 1 1k 1k\n");
    double db = r.rows[0][col(r, "vdb(d)")], ph = r.rows[0][col(r, "vp(d)")];
    std::printf("  AC CS amp: %.3f dB (expect 18.062), phase %.1f deg\n", db, ph);
    CHECK(near(db, 18.062, 0.01), "common-source gain = -gm*RD");
    CHECK(std::fabs(ph) > 179.9, "common-source output is inverted");
}

static void testAcNeedsStimulus() {
    bool threw = false;
    try { run("V1 in 0 DC 1\nR1 in 0 1k\n.ac dec 10 1 1k\n"); } catch (const SimError&) { threw = true; }
    CHECK(threw, ".ac without an AC source is rejected");
}

static void testSubckt() {
    Result r = run(".subckt divider in out\nR1 in out 1k\nR2 out 0 1k\n.ends\n"
                   "V1 a 0 DC 8\nX1 a b divider\nX2 b c divider\n.op\n");
    CHECK(near(r.rows[0][col(r, "v(b)")], 3.2, 1e-6), "loaded subckt divider v(b) = 3.2");
    CHECK(near(r.rows[0][col(r, "v(c)")], 1.6, 1e-6), "second instance v(c) = 1.6");

    // internal nodes get hierarchical names, nested subckts expand
    Result n = run(".subckt rc in out\nR1 in mid 1k\nR2 mid out 1k\n.ends\n.subckt two a b\nX1 a m rc\nX2 m b rc\n.ends\n"
                   "V1 s 0 DC 4\nXT s t two\nR9 t 0 4k\n.op\n");
    bool hasMid = false;
    for (auto& c : n.columns) if (c == "v(xt.x1.mid)") hasMid = true;
    CHECK(hasMid, "hierarchical internal node name xt.x1.mid");
    CHECK(near(n.rows[0][col(n, "v(t)")], 2.0, 1e-6), "nested subckt chain = 4k series + 4k load");

    bool t1 = false, t2 = false, t3 = false;
    try { parseNetlist("X1 a b nosuch\n"); } catch (const ParseError&) { t1 = true; }
    try { parseNetlist(".subckt s a b\n.ends\nX1 a s\n"); } catch (const ParseError&) { t2 = true; }
    try { parseNetlist(".subckt loop n\nX1 n loop\n.ends\nX0 1 loop\n"); } catch (const ParseError&) { t3 = true; }
    CHECK(t1, "unknown subcircuit rejected");
    CHECK(t2, "wrong port count rejected");
    CHECK(t3, "recursive subcircuit rejected");
}

static void testProbe() {
    Result r = run("V1 in 0 DC 6\nR1 in out 1k\nR2 out 0 2k\n.probe v(out)\n.op\n");
    CHECK(r.columns.size() == 1 && r.columns[0] == "v(out)", ".probe limits the output columns");
    CHECK(near(r.rows[0][0], 4.0, 1e-6), "probed value correct");
    bool threw = false;
    try { run("V1 in 0 DC 1\nR1 in 0 1k\n.probe v(nothere)\n.op\n"); } catch (const SimError&) { threw = true; }
    CHECK(threw, "unknown probe rejected");
}

static std::string gridNetlist(int n, bool dense) {
    std::string s = "VDD vdd 0 DC 1.0\nRPAD vdd n0_0 0.05\n";
    auto nm = [](int i, int j) { return "n" + std::to_string(i) + "_" + std::to_string(j); };
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (j + 1 < n) s += "RH" + std::to_string(i) + "_" + std::to_string(j) + " " + nm(i, j) + " " + nm(i, j + 1) + " 0.5\n";
            if (i + 1 < n) s += "RV" + std::to_string(i) + "_" + std::to_string(j) + " " + nm(i, j) + " " + nm(i + 1, j) + " 0.5\n";
            s += "I" + std::to_string(i) + "_" + std::to_string(j) + " " + nm(i, j) + " 0 DC 50u\n";
        }
    if (dense) s += ".options solver=dense\n";
    return s + ".op\n";
}

static void testSparseMatchesDense() {
    Result a = run(gridNetlist(20, false).c_str()), b = run(gridNetlist(20, true).c_str());
    double worst = 0;
    for (size_t i = 0; i < a.rows[0].size(); ++i) worst = std::max(worst, std::fabs(a.rows[0][i] - b.rows[0][i]));
    std::printf("  20x20 grid: max |sparse - dense| = %.2e V (%d unknowns)\n", worst, a.size);
    CHECK(worst < 1e-9, "sparse and dense solvers agree on a resistor grid");
}

static void testLargeGrid() {
    // 60x60 = 3,602 unknowns. Conservation: the supply must deliver exactly the total load current.
    Result r = run(gridNetlist(60, false).c_str());
    double isrc = r.rows[0][col(r, "i(vdd)")], expected = -3600 * 50e-6;
    double vmin = 1e9;
    for (size_t i = 0; i < r.rows[0].size() - 1; ++i) vmin = std::min(vmin, r.rows[0][i]);
    std::printf("  60x60 grid (%d unknowns): i(vdd)=%.6f A (expect %.6f), worst node %.4f V, fill %.1fx, factor %.3fs\n",
                r.size, isrc, expected, vmin, double(r.nnzLU) / double(r.nnzA), r.factorSeconds);
    CHECK(std::fabs(isrc - expected) / std::fabs(expected) < 1e-6, "supply current equals total load (KCL)");
    CHECK(vmin < 1.0 && vmin > 0.5, "IR drop is plausible");
}

static void testFactorizationReuse() {
    // Linear transient: only the t=0 DC solve, one backward-Euler step and the trapezoidal matrix are ever factored.
    Result r = run("V1 in 0 PULSE(0 5 0 0 0 1 2)\nR1 in out 1k\nC1 out 0 1u\n.tran 5u 5m\n");
    std::printf("  RC tran: %ld factorizations for %ld solves\n", r.factorizations, r.solves);
    CHECK(r.factorizations <= 3, "linear transient reuses factorisations");
    CHECK(r.solves >= 1000, "one solve per time step");
}

static void testErrors() {
    bool threw = false;
    try { run("V1 a 0 DC 1\nV2 a 0 DC 2\n.op\n"); } catch (const SimError&) { threw = true; }
    CHECK(threw, "parallel voltage sources rejected as singular");
    threw = false;
    try { parseNetlist("X1 a b c\n"); } catch (const ParseError&) { threw = true; }
    CHECK(threw, "unknown element rejected");
    threw = false;
    try { parseNetlist("M1 d g s b nomodel\n"); } catch (const ParseError&) { threw = true; }
    CHECK(threw, "unknown model rejected");
}

int main() {
    testParseValues();
    testDivider();
    testRcStep(true);
    testRcStep(false);
    testRlStep();
    testDiode();
    testCmosInverter();
    testRectifier();
    testAcRc();
    testAcRlc();
    testAcMosAmp();
    testAcNeedsStimulus();
    testSubckt();
    testProbe();
    testSparseMatchesDense();
    testLargeGrid();
    testFactorizationReuse();
    testErrors();
    if (failures) { std::printf("%d check(s) failed\n", failures); return 1; }
    std::printf("all tests passed\n");
    return 0;
}
