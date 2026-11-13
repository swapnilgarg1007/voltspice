#include "parser.hpp"
#include "simulator.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

static std::string jsonEscape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if (c == '\n') o += "\\n";
        else if (static_cast<unsigned char>(c) >= 0x20) o += c;
    }
    return o;
}

static std::string num(double v) {
    if (!std::isfinite(v)) return "null";
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.9g", v);
    return buf;
}

int main(int argc, char** argv) {
    bool json = false, stats = false;
    std::string path;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--json") json = true;
        else if (a == "--stats") stats = true;
        else if (a == "-h" || a == "--help") {
            std::cout << "usage: voltspice [--json] [--stats] <netlist.cir | ->\n";
            return 0;
        } else path = a;
    }
    if (path.empty()) { std::cerr << "usage: voltspice [--json] [--stats] <netlist.cir | ->\n"; return 2; }

    std::stringstream buf;
    if (path == "-") buf << std::cin.rdbuf();
    else {
        std::ifstream f(path);
        if (!f) { std::cerr << "error: cannot open " << path << "\n"; return 2; }
        buf << f.rdbuf();
    }

    try {
        vs::Result r = vs::simulate(vs::parseNetlist(buf.str()));
        if (json) {
            std::cout << "{\"analysis\":\"" << r.analysis << "\",\"columns\":[";
            for (size_t i = 0; i < r.columns.size(); ++i)
                std::cout << (i ? "," : "") << '"' << jsonEscape(r.columns[i]) << '"';
            std::cout << "],\"rows\":[";
            for (size_t i = 0; i < r.rows.size(); ++i) {
                std::cout << (i ? ",[" : "[");
                for (size_t j = 0; j < r.rows[i].size(); ++j) std::cout << (j ? "," : "") << num(r.rows[i][j]);
                std::cout << "]";
            }
            std::cout << "],\"stats\":{\"unknowns\":" << r.size << ",\"nnzA\":" << r.nnzA << ",\"nnzLU\":" << r.nnzLU
                      << ",\"factorSeconds\":" << num(r.factorSeconds) << ",\"solveSeconds\":" << num(r.solveSeconds)
                      << "}}\n";
        } else {
            for (size_t i = 0; i < r.columns.size(); ++i) std::cout << (i ? "," : "") << r.columns[i];
            std::cout << "\n";
            for (auto& row : r.rows) {
                for (size_t j = 0; j < row.size(); ++j) std::cout << (j ? "," : "") << num(row[j]);
                std::cout << "\n";
            }
        }
        if (stats)
            std::fprintf(stderr,
                         "unknowns=%d nnz(A)=%zu nnz(LU)=%zu fill=%.2fx factorizations=%ld solves=%ld "
                         "factor=%.3fs solve=%.3fs\n",
                         r.size, r.nnzA, r.nnzLU, r.nnzA ? double(r.nnzLU) / double(r.nnzA) : 0.0, r.factorizations,
                         r.solves, r.factorSeconds, r.solveSeconds);
    } catch (const std::exception& e) {
        if (json) std::cout << "{\"error\":\"" << jsonEscape(e.what()) << "\"}\n";
        else std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
