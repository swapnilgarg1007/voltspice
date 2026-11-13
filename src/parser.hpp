#pragma once
#include "circuit.hpp"
#include <string>

namespace vs {

// Parses an engineering-notation number such as "4.7k", "10meg", "2.2u". Throws ParseError.
double parseValue(const std::string& token);

// Parses a SPICE-style netlist (case-insensitive; '*' comment lines; '+' continuation lines).
// Elements: R C L V I D M X.  Sources: DC / PULSE / SIN plus an optional "AC mag [phase]".
// Directives: .model .subckt/.ends .op .dc .tran .ac .probe .options .end
// Subcircuit instances (X) are flattened at parse time; internal names become "x1.node".
Circuit parseNetlist(const std::string& text);

}  // namespace vs
