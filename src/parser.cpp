#include "parser.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <sstream>

namespace vs {

double parseValue(const std::string& tok) {
    if (tok.empty()) throw ParseError("empty number");
    char* end = nullptr;
    double v = std::strtod(tok.c_str(), &end);
    if (end == tok.c_str()) throw ParseError("bad number: '" + tok + "'");
    std::string suf(end);
    if (suf.empty()) return v;
    double m = 1;
    if (suf.compare(0, 3, "meg") == 0) m = 1e6;
    else if (suf.compare(0, 3, "mil") == 0) m = 25.4e-6;
    else switch (suf[0]) {
        case 't': m = 1e12; break;
        case 'g': m = 1e9; break;
        case 'k': m = 1e3; break;
        case 'm': m = 1e-3; break;
        case 'u': m = 1e-6; break;
        case 'n': m = 1e-9; break;
        case 'p': m = 1e-12; break;
        case 'f': m = 1e-15; break;
        default: throw ParseError("bad number suffix: '" + tok + "'");
    }
    return v * m;  // trailing unit letters ("1kohm", "10uF") are ignored, as in SPICE
}

namespace {

using Lines = std::vector<std::pair<int, std::string>>;

struct Subckt {
    std::vector<std::string> ports;
    Lines body;
};

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> tokenize(std::string line) {
    line = lower(line);
    for (char& c : line)
        if (c == '(' || c == ')' || c == ',' || c == '=') c = ' ';
    std::istringstream ss(line);
    std::vector<std::string> t;
    for (std::string w; ss >> w;) t.push_back(w);
    return t;
}

std::map<std::string, double> keyValues(const std::vector<std::string>& t, size_t from) {
    std::map<std::string, double> kv;
    for (size_t i = from; i + 1 < t.size(); i += 2) kv[t[i]] = parseValue(t[i + 1]);
    return kv;
}

void parseSource(Element& e, const std::vector<std::string>& t) {
    if (t.size() < 4) throw ParseError("source needs: name n+ n- value");
    size_t i = 3;
    if (t[i] == "ac") {
        e.wf.kind = Waveform::DC;
        e.wf.dc = 0;
    } else if (t[i] == "dc") {
        if (t.size() < 5) throw ParseError("DC needs a value");
        e.wf.kind = Waveform::DC;
        e.wf.dc = parseValue(t[i + 1]);
    } else if (t[i] == "pulse") {
        e.wf.kind = Waveform::PULSE;
        if (t.size() < i + 1 + 7) throw ParseError("PULSE needs: v1 v2 td tr tf pw per");
        for (int k = 0; k < 7; ++k) e.wf.p[k] = parseValue(t[i + 1 + k]);
    } else if (t[i] == "sin") {
        e.wf.kind = Waveform::SIN;
        if (t.size() < i + 1 + 3) throw ParseError("SIN needs: vo va freq [td]");
        for (int k = 0; k < 4 && i + 1 + k < t.size(); ++k) {
            try { e.wf.p[k] = parseValue(t[i + 1 + k]); } catch (const ParseError&) { break; }
        }
    } else {
        e.wf.kind = Waveform::DC;
        e.wf.dc = parseValue(t[i]);
    }
    for (size_t k = i; k < t.size(); ++k) {  // optional: AC mag [phase_degrees]
        if (t[k] != "ac") continue;
        if (k + 1 >= t.size()) throw ParseError("AC needs a magnitude");
        e.wf.acMag = parseValue(t[k + 1]);
        if (k + 2 < t.size()) {
            try { e.wf.acPhase = parseValue(t[k + 2]); } catch (const ParseError&) {}
        }
        break;
    }
}

class Builder {
public:
    Builder(Circuit& c, const std::map<std::string, Subckt>& subs) : c_(c), subs_(subs) {}

    void process(const Lines& lines, const std::string& prefix, const std::map<std::string, std::string>& ports,
                 int depth) {
        if (depth > 20) throw ParseError("subcircuit nesting too deep (recursive definition?)");
        for (const auto& [no, line] : lines) {
            auto t = tokenize(line);
            if (t.empty()) continue;
            try {
                if (t[0][0] == '.') {
                    if (t[0] == ".end") { if (depth == 0) return; continue; }
                    if (t[0] == ".model" || depth == 0) directive(t, line);
                    continue;
                }
                element(t, prefix, ports, depth);
            } catch (const ParseError& err) {
                std::string msg = err.what();
                if (msg.rfind("line ", 0) != 0) msg = "line " + std::to_string(no) + ": " + msg;
                throw ParseError(msg);
            }
        }
    }

private:
    Circuit& c_;
    const std::map<std::string, Subckt>& subs_;

    void directive(const std::vector<std::string>& t, const std::string& rawLine) {
        const std::string& head = t[0];
        Analysis& a = c_.analysis;
        if (head == ".op") {
            a.kind = Analysis::OP;
        } else if (head == ".dc") {
            if (t.size() < 5) throw ParseError(".dc needs: source start stop step");
            a.kind = Analysis::DC;
            a.sweepSource = t[1];
            a.start = parseValue(t[2]); a.stop = parseValue(t[3]); a.step = parseValue(t[4]);
            if (a.step == 0) throw ParseError(".dc step must be non-zero");
        } else if (head == ".tran") {
            if (t.size() < 3) throw ParseError(".tran needs: tstep tstop");
            a.kind = Analysis::TRAN;
            a.tstep = parseValue(t[1]); a.tstop = parseValue(t[2]);
            if (a.tstep <= 0 || a.tstop <= 0) throw ParseError(".tran times must be positive");
        } else if (head == ".ac") {
            if (t.size() < 5) throw ParseError(".ac needs: dec|oct|lin points fstart fstop");
            if (t[1] != "dec" && t[1] != "oct" && t[1] != "lin") throw ParseError(".ac sweep must be dec, oct or lin");
            a.kind = Analysis::AC;
            a.acSweep = t[1];
            a.acPoints = static_cast<int>(parseValue(t[2]));
            a.fstart = parseValue(t[3]); a.fstop = parseValue(t[4]);
            if (a.acPoints < 1 || a.fstart <= 0 || a.fstop < a.fstart) throw ParseError("bad .ac range");
        } else if (head == ".model") {
            if (t.size() < 3) throw ParseError(".model needs: name type [params]");
            Model m;
            m.type = t[2];
            m.params = keyValues(t, 3);
            c_.models[t[1]] = m;
        } else if (head == ".probe") {
            std::istringstream ss(lower(rawLine));  // keep "v(out)" intact
            std::string w;
            ss >> w;
            while (ss >> w) c_.probes.push_back(w);
        } else if (head == ".options") {
            for (size_t i = 1; i + 1 < t.size(); i += 2) {
                if (t[i] == "method") c_.trapezoidal = (t[i + 1] != "be");
                if (t[i] == "solver") c_.denseSolver = (t[i + 1] == "dense");
            }
        }
        // other directives (.title, .print, ...) are ignored
    }

    void element(const std::vector<std::string>& t, const std::string& prefix,
                 const std::map<std::string, std::string>& ports, int depth) {
        const std::string& head = t[0];
        auto nodeOf = [&](const std::string& n) -> int {
            if (n == "0" || n == "gnd") return c_.node("0");
            auto it = ports.find(n);
            if (it != ports.end()) return c_.node(it->second);
            return c_.node(prefix + n);
        };
        auto need = [&](size_t n) { if (t.size() < n) throw ParseError("too few fields for " + head); };

        if (head[0] == 'x') {  // subcircuit instance: xname n1 ... nk subname
            need(3);
            const std::string& subName = t.back();
            auto it = subs_.find(subName);
            if (it == subs_.end()) throw ParseError("unknown subcircuit '" + subName + "'");
            const Subckt& sc = it->second;
            if (t.size() - 2 != sc.ports.size())
                throw ParseError(head + ": subcircuit '" + subName + "' expects " + std::to_string(sc.ports.size()) +
                                 " nodes, got " + std::to_string(t.size() - 2));
            std::map<std::string, std::string> map;
            for (size_t k = 0; k < sc.ports.size(); ++k) map[sc.ports[k]] = c_.nodeName(nodeOf(t[1 + k]));
            process(sc.body, prefix + head + ".", map, depth + 1);
            return;
        }

        Element e;
        e.name = prefix + head;
        switch (head[0]) {
            case 'r': case 'c': case 'l':
                need(4);
                e.type = head[0] == 'r' ? ElemType::R : head[0] == 'c' ? ElemType::C : ElemType::L;
                e.nodes = {nodeOf(t[1]), nodeOf(t[2])};
                e.value = parseValue(t[3]);
                if (e.type == ElemType::R && e.value == 0) throw ParseError("resistor cannot be 0 ohms");
                break;
            case 'v': case 'i':
                e.type = head[0] == 'v' ? ElemType::V : ElemType::I;
                e.nodes = {nodeOf(t[1]), nodeOf(t[2])};
                parseSource(e, t);
                break;
            case 'd':
                need(4);
                e.type = ElemType::D;
                e.nodes = {nodeOf(t[1]), nodeOf(t[2])};
                e.model = t[3];
                break;
            case 'm': {
                need(6);
                e.type = ElemType::M;
                e.nodes = {nodeOf(t[1]), nodeOf(t[2]), nodeOf(t[3])};  // bulk (t[4]) is ignored
                e.model = t[5];
                auto kv = keyValues(t, 6);
                if (kv.count("w")) e.w = kv["w"];
                if (kv.count("l")) e.l = kv["l"];
                break;
            }
            default:
                throw ParseError("unknown element '" + head + "'");
        }
        if (!names_.insert(e.name).second) throw ParseError("duplicate element name '" + e.name + "'");
        c_.elements.push_back(e);
    }

    std::set<std::string> names_;
};

}  // namespace

Circuit parseNetlist(const std::string& text) {
    // Pass 1: join continuation lines, strip comments.
    Lines lines;
    std::istringstream in(text);
    std::string raw;
    int lineNo = 0;
    while (std::getline(in, raw)) {
        ++lineNo;
        size_t semi = raw.find(';');
        if (semi != std::string::npos) raw.erase(semi);
        size_t first = raw.find_first_not_of(" \t\r");
        if (first == std::string::npos || raw[first] == '*') continue;
        if (raw[first] == '+' && !lines.empty()) {
            lines.back().second += " " + raw.substr(first + 1);
            continue;
        }
        lines.emplace_back(lineNo, raw);
    }

    // Pass 2: pull out .subckt ... .ends blocks.
    std::map<std::string, Subckt> subs;
    Lines top;
    Subckt* cur = nullptr;
    std::string curName;
    for (auto& ln : lines) {
        auto t = tokenize(ln.second);
        if (t.empty()) continue;
        try {
            if (t[0] == ".subckt") {
                if (cur) throw ParseError("nested .subckt definitions are not supported");
                if (t.size() < 2) throw ParseError(".subckt needs a name");
                curName = t[1];
                if (subs.count(curName)) throw ParseError("duplicate subcircuit '" + curName + "'");
                cur = &subs[curName];
                cur->ports.assign(t.begin() + 2, t.end());
            } else if (t[0] == ".ends") {
                if (!cur) throw ParseError(".ends without .subckt");
                cur = nullptr;
            } else if (cur) {
                cur->body.push_back(ln);
            } else {
                top.push_back(ln);
            }
        } catch (const ParseError& err) {
            throw ParseError("line " + std::to_string(ln.first) + ": " + err.what());
        }
    }
    if (cur) throw ParseError("missing .ends for subcircuit '" + curName + "'");

    Circuit c;
    Builder(c, subs).process(top, "", {}, 0);

    for (const auto& e : c.elements)
        if (e.type == ElemType::M && !c.models.count(e.model))
            throw ParseError("element " + e.name + ": unknown model '" + e.model + "'");
    return c;
}

}  // namespace vs
