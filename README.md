# voltspice

A SPICE-style circuit simulator written from scratch in C++17, with a browser interface for plotting results. It reads a netlist, builds the circuit equations with **Modified Nodal Analysis (MNA)**, and solves them with its own **sparse LU solver** and Newton-Raphson iteration. It handles circuits with tens of thousands of nodes.

No simulation or linear-algebra libraries are used: the parser, matrix assembly, sparse solver, nonlinear iteration, time integration and AC analysis are all in `src/`.

## Features
- **Elements:** resistor, capacitor, inductor, independent V/I sources (DC, `PULSE`, `SIN`, plus `AC mag [phase]`), diode, MOSFET (square-law with channel-length modulation, NMOS and PMOS)
- **Analyses:** `.op`, `.dc` sweep, `.tran` (fixed step, trapezoidal or backward Euler), `.ac` (small-signal frequency sweep: `dec`, `oct`, `lin`)
- **Hierarchy:** `.subckt` / `.ends` with nested instances (`X` elements), flattened at parse time with names like `xt.x1.mid`
- **Sparse solver:** Markowitz-style pivoting with threshold pivoting, real and complex, with the factorisation reused across time steps and sweep points of linear circuits
- **Nonlinear solving:** Newton-Raphson with voltage-step limiting, an exponent-safe diode model, and a source-stepping fallback
- **Large outputs under control:** `.probe v(node) i(vname)` selects signals, and oversized outputs are refused with a clear message
- **Interfaces:** CLI (CSV or JSON, `--stats` for solver statistics) and a web UI with transient/DC plots, Bode plots, an example library and a hover readout

## Build and run
```bash
cmake -S . -B build && cmake --build build
cd build && ctest --output-on-failure && cd ..

./build/voltspice examples/rc_lowpass.cir               # CSV to stdout
./build/voltspice --json examples/ac_rc_bode.cir
./build/voltspice --stats examples/power_grid_ir_drop.cir >/dev/null

python3 server.py 8000                                  # web UI at http://localhost:8000
python3 tools/gen_grid.py 200 > grid200.cir             # 40,000-node benchmark circuit
```

## Netlist examples
```
* AC analysis: RC low-pass Bode plot
V1 in 0 DC 0 AC 1
R1 in out 1k
C1 out 0 159.1549n
.ac dec 40 10 1meg
```
```
* Hierarchy
.subckt divider in out
R1 in out 1k
R2 out 0 1k
.ends
V1 a 0 DC 8
X1 a b divider
X2 b c divider
.op
```
AC results are reported as `vdb(x)` (magnitude in dB) and `vp(x)` (phase in degrees). Source-current sign follows SPICE: `i(V1)` is positive when current flows into the `+` terminal, so a source delivering power shows a negative current.

## Validation
Measured by the test suite (`ctest`, 2 executables):

| Circuit | Check | Result |
|---|---|---|
| Resistive divider | exact `v(out)` | pass |
| RC step (tau = 1 ms) | vs `5(1 - e^(-t/tau))` | max error 6.2e-5 V (trapezoidal), 4.6e-3 V (backward Euler) |
| RL step | vs `(V/R)(1 - e^(-t/tau))` | max error 1.2e-5 A |
| Diode + resistor | KCL residual | below 1e-3 (Vd = 0.6925 V) |
| Symmetric CMOS inverter | switching threshold = VDD/2 | Vm = 0.900 V |
| RC low-pass, AC | -3.01 dB / -45 deg at fc, -20.04 dB / -84.3 deg at 10 fc | -3.010 dB / -45.00 deg, -20.043 dB / -84.29 deg |
| Series RLC, AC | resonant peak vs closed form `Q/sqrt(1-1/(4Q^2))` | 10.11 dB at 4,898 Hz (closed form: 10.11 dB at 4,906 Hz) |
| NMOS common-source amp, AC | gain = -gm*RD = -8 | 18.062 dB, 180 deg (needs the DC bias point) |
| Resistor grids | sparse vs dense solver, 402 unknowns | max difference 2.6e-13 V |
| 60x60 power grid (3,602 unknowns) | supply current = total load current | -0.180000 A vs -0.180000 A |
| Random MNA-style matrices | sparse LU vs dense LU, 200 systems | max difference 5.5e-14 |
| Subcircuits, probes, errors | hierarchy, nested names, recursion, bad input | pass |

**TODO:** compare waveforms against LTSpice for the same netlists and add overlay plots to `docs/`.

## Performance
Measured on a single core (`tools/bench.sh`, `tools/gen_grid.py`). Circuit: N x N resistive power grid with a load current at every node.

| Circuit | Unknowns | Time | Notes |
|---|---|---|---|
| DC power grid 100x100 | 10,002 | 0.30 s | fill-in 8.6x |
| DC power grid 200x200 | 40,002 | 1.8 s | 2.4M non-zeros in L+U |
| Transient RC mesh 100x100, 150 steps | 10,002 | 1.1 s | 3 factorisations, 151 solves |

Sparse vs the dense reference solver (`.options solver=dense`) on the same grids:

| Unknowns | Sparse | Dense |
|---|---|---|
| 2,502 | 0.05 s | 0.15 s |
| 4,902 | 0.12 s | 0.80 s |
| 8,102 | 0.23 s | 2.52 s |

The dense solver also needs n^2 memory (over 500 MB at 8,100 unknowns), so it cannot reach the sizes above. Fill-in of 8x to 12x is acceptable but not optimal; a proper minimum-degree or nested-dissection ordering would reduce it (see roadmap). Nonlinear circuits refactor every Newton iteration.

## Limitations
- Fixed time step only (no adaptive step control)
- Pivot ordering is a simple Markowitz heuristic, not AMD/COLAMD
- No `.param` / parameterised subcircuits; nested `.subckt` definitions are not supported
- No MOSFET body effect or capacitances, no diode capacitance or series resistance; the MOSFET bulk terminal is parsed but ignored
- No BJT; no noise or pole-zero analysis
- Newton failures on very stiff circuits may need manual help (smaller step, `.options method=be`)

## Roadmap
- [ ] Adaptive time stepping with local-truncation-error control
- [ ] Better fill-reducing ordering (AMD), and refactoring that reuses the symbolic ordering for nonlinear circuits
- [ ] WebAssembly build (Emscripten) so the UI runs without the Python server
- [ ] `.param` and parameterised subcircuits
- [ ] BJT model, MOSFET capacitances and a short-channel model
- [ ] LTSpice comparison plots and a regression suite
- [ ] Schematic editor in the browser that generates netlists

## Project layout
```
src/        circuit.*   data model + sources      parser.*    netlist parser (.subckt flattening)
            sparse.hpp  sparse LU (real/complex)  solver.*    dense LU (reference)
            simulator.* MNA, Newton, .op/.dc/.tran/.ac        main.cpp  CLI
web/        index.html, app.js (canvas plotter, no external dependencies)
server.py   dev server: static files + /api/simulate (wraps the CLI)
examples/   ready-to-run netlists (RC, RLC, rectifier, CMOS inverter, ring oscillator, amplifier, power grid)
tools/      gen_grid.py (benchmark netlist generator), bench.sh
tests/      validation suite
```

## License
MIT
