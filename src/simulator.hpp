#pragma once
#include "circuit.hpp"
#include <stdexcept>
#include <string>
#include <vector>

namespace vs {

struct SimError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct Result {
    std::string analysis;                    // "op", "dc", "tran" or "ac"
    std::vector<std::string> columns;        // first column is the sweep variable (none for op)
    std::vector<std::vector<double>> rows;   // ac columns are vdb(x)/vp(x) = magnitude in dB, phase in degrees

    // Solver statistics
    int size = 0;                 // number of MNA unknowns
    size_t nnzA = 0, nnzLU = 0;   // matrix non-zeros and non-zeros after factorisation
    long factorizations = 0, solves = 0;
    double factorSeconds = 0, solveSeconds = 0;
};

// Runs the analysis requested in the circuit (defaults to .op). Throws SimError.
Result simulate(Circuit circuit);

}  // namespace vs
