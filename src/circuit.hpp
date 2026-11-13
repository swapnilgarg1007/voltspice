#pragma once
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace vs {

enum class ElemType { R, C, L, V, I, D, M };

// Time-dependent value of an independent source.
struct Waveform {
    enum Kind { DC, PULSE, SIN } kind = DC;
    double dc = 0;
    double p[7] = {0, 0, 0, 0, 0, 0, 0};  // PULSE: v1 v2 td tr tf pw per | SIN: vo va freq td
    double acMag = 0, acPhase = 0;        // small-signal stimulus for .ac (phase in degrees)
    double dcValue() const;               // value used for operating-point analysis
    double eval(double t) const;          // value used during transient analysis
};

struct Model {
    std::string type;  // d, nmos, pmos
    std::map<std::string, double> params;
    double get(const std::string& k, double def) const {
        auto it = params.find(k);
        return it == params.end() ? def : it->second;
    }
};

struct Element {
    ElemType type;
    std::string name;
    std::vector<int> nodes;  // R,C,L,V,I,D: {n+, n-}; M: {drain, gate, source}
    double value = 0;        // R ohms, C farads, L henries
    Waveform wf;             // V and I sources
    std::string model;       // D and M
    double w = 1e-6, l = 1e-6;
};

struct Analysis {
    enum Kind { NONE, OP, DC, TRAN, AC } kind = NONE;
    std::string sweepSource;
    double start = 0, stop = 0, step = 0;  // .dc
    double tstep = 0, tstop = 0;           // .tran
    std::string acSweep = "dec";           // .ac: dec | oct | lin
    int acPoints = 10;                     // points per decade/octave, or total for lin
    double fstart = 0, fstop = 0;
};

class Circuit {
public:
    Circuit() { nodeNames_.push_back("0"); nodeIndex_["0"] = 0; }
    int node(const std::string& name);  // creates the node on first use; "0"/"gnd" -> ground
    int numNodes() const { return static_cast<int>(nodeNames_.size()); }
    const std::string& nodeName(int i) const { return nodeNames_[i]; }

    std::vector<Element> elements;
    std::map<std::string, Model> models;
    Analysis analysis;
    bool trapezoidal = true;  // .options method=trap|be
    bool denseSolver = false; // .options solver=dense (reference solver, for cross-checking)
    std::vector<std::string> probes;  // .probe v(node) i(vname): limits which signals are output

private:
    std::vector<std::string> nodeNames_;
    std::map<std::string, int> nodeIndex_;
};

struct ParseError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

}  // namespace vs
