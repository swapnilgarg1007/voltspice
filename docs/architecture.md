# Architecture

## Pipeline
`netlist text -> parser -> Circuit -> simulator (MNA + Newton + time stepping) -> Result -> CSV/JSON -> UI`

## Modified Nodal Analysis
Unknowns: one voltage per non-ground node, plus one current per voltage source and per inductor (branch currents). The system is `A x = b`.

Stamps (n+ / n- are the element's terminals, `a`/`b` below are matrix indices, ground rows are skipped):

| Element | Stamp |
|---|---|
| Resistor `R` | conductance `g = 1/R`: `+g` on (n+,n+) and (n-,n-), `-g` on (n+,n-) and (n-,n+) |
| Current source `I` | `b[n+] -= I`, `b[n-] += I` |
| Voltage source `V` | branch row: `v(n+) - v(n-) = V`; KCL columns `+1` at n+, `-1` at n- |
| Inductor `L` (DC) | like a 0 V source (short circuit) |
| Capacitor `C` (DC) | open circuit (no stamp) |

A tiny `gmin = 1e-12 S` from every node to ground keeps floating nodes (e.g. a node connected only through capacitors in DC) solvable.

## Transient companion models
With step `h`, previous voltage `v0` and current `i0`:

| | Backward Euler | Trapezoidal |
|---|---|---|
| Capacitor | `g = C/h`, `i = g(v - v0)` | `g = 2C/h`, `i = g(v - v0) - i0` |
| Inductor | `v = (L/h)(i - i0)` | `v = (2L/h)(i - i0) - v0` |

The first step always uses backward Euler: trapezoidal integration rings when a source jumps at t = 0. The initial state comes from a DC operating point with sources at their t = 0 value (`PULSE` -> v1, `SIN` -> offset).

## Nonlinear devices (Newton-Raphson)
Each iteration linearises a device at the current guess into a conductance plus an equivalent current source, stamps it, and re-solves. For a diode: `gd = dI/dV`, `ieq = I(v0) - gd*v0`.

- **Diode:** `I = Is(exp(v/nVt) - 1)`; above `v/nVt = 40` the exponential is replaced by its tangent line so it can't overflow.
- **MOSFET:** cutoff / triode / saturation with `(1 + lambda*Vds)`. The device is symmetric, so drain and source swap when `Vds` flips sign. PMOS reuses the NMOS equations by multiplying voltages by a polarity factor of -1.
- **Convergence:** `|dx| <= 1e-3*|x| + abstol` on every unknown, node-voltage steps limited to 1 V per iteration, 150 iterations max. If the operating point fails, sources are ramped from 0 to 100% in 25 steps.

## Solver
Dense LU with partial pivoting (`solver.cpp`). A pivot below `1e-14 * max|A|` is reported as a singular matrix, which catches voltage-source loops. Replacing this with sparse storage is the first roadmap item.

## Sparse LU (`src/sparse.hpp`)
MNA matrices are very sparse (a few entries per row) but have zeros on the diagonal where voltage sources and inductors sit, so a naive ordering fails. The solver:

1. Collects triplets, sorts them and sums duplicates into row-sparse storage.
2. At each step picks the active **column with the fewest entries** (tracked in an ordered set), then, among its rows that pass a threshold test (`|a| >= 0.1 * column max`), the row with the fewest entries. Few entries means little fill-in; the threshold keeps it numerically stable.
3. Eliminates that column from the other rows by merging sorted rows, creating fill-in only where needed.
4. Stores the multipliers (L) and pivot rows (U). `solve()` replays the multipliers on the right-hand side and back-substitutes, so one factorisation serves many right-hand sides.

The class is templated on `double` and `std::complex<double>`; AC analysis uses the complex version.

**When the factorisation is reused:** a linear circuit's matrix depends only on the analysis type, the integration method and the time step. So `.tran` factors the matrix once for the backward-Euler first step and once for the trapezoidal steps, and `.dc` sweeps of linear circuits factor once. Nonlinear circuits change the matrix every Newton iteration and refactor each time.

## AC analysis
1. Solve the DC operating point (this fixes the bias of diodes and MOSFETs).
2. For each frequency, stamp the circuit again with `omega = 2*pi*f`: resistors as before, capacitors as `j*omega*C` (imaginary conductance), inductors with `-j*omega*L` on the branch diagonal, diodes and MOSFETs as their linearised conductances (`gd`, `gm`, `gds`) with no equivalent current sources, and only the sources' `AC mag phase` as stimulus.
3. Solve the complex system with the complex sparse LU and report `20*log10|x|` and `arg(x)`.

## Subcircuits
`.subckt` blocks are collected in a first pass. Each `X` instance is expanded recursively while parsing: ports map to the instance's nodes, internal nodes and element names get the instance path as a prefix (`xt.x1.mid`, `xt.x1.r1`), and ground stays global. A depth limit of 20 stops recursive definitions.
