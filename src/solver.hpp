#pragma once
#include <vector>

namespace vs {

// Solves A x = b in place by LU decomposition with partial pivoting.
// A is n*n, row-major. On success b holds x and true is returned.
// Returns false if the matrix is singular (e.g. a loop of voltage sources).
bool luSolve(std::vector<double>& A, std::vector<double>& b, int n);

}  // namespace vs
