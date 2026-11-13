#include "circuit.hpp"
#include <cmath>

namespace vs {

int Circuit::node(const std::string& name) {
    std::string n = (name == "gnd") ? "0" : name;
    auto it = nodeIndex_.find(n);
    if (it != nodeIndex_.end()) return it->second;
    int idx = static_cast<int>(nodeNames_.size());
    nodeNames_.push_back(n);
    nodeIndex_[n] = idx;
    return idx;
}

double Waveform::dcValue() const {
    switch (kind) {
        case DC: return dc;
        case PULSE: return p[0];  // v1
        case SIN: return p[0];    // vo
    }
    return 0;
}

double Waveform::eval(double t) const {
    const double kPi = 3.14159265358979323846;
    switch (kind) {
        case DC: return dc;
        case PULSE: {
            double v1 = p[0], v2 = p[1], td = p[2], tr = p[3], tf = p[4], pw = p[5], per = p[6];
            if (t < td) return v1;
            double tt = t - td;
            if (per > 0) tt = std::fmod(tt, per);
            if (tt < tr) return v1 + (v2 - v1) * tt / tr;
            if (tt < tr + pw) return v2;
            if (tt < tr + pw + tf) return v2 + (v1 - v2) * (tt - tr - pw) / tf;
            return v1;
        }
        case SIN: {
            double vo = p[0], va = p[1], f = p[2], td = p[3];
            if (t < td) return vo;
            return vo + va * std::sin(2 * kPi * f * (t - td));
        }
    }
    return 0;
}

}  // namespace vs
