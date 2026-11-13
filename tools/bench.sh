#!/bin/bash
# Compare the sparse solver against the dense reference on N x N power grids.
#   cmake -S . -B build && cmake --build build && tools/bench.sh
cd "$(dirname "$0")/.."
for N in ${@:-50 70 90}; do
  python3 tools/gen_grid.py $N > /tmp/g$N.cir
  sed 's/^\.op/.options solver=dense\n.op/' /tmp/g$N.cir > /tmp/g${N}d.cir
  s=$(date +%s.%N); ./build/voltspice /tmp/g$N.cir > /dev/null;  sp=$(echo "$(date +%s.%N) - $s" | bc)
  s=$(date +%s.%N); ./build/voltspice /tmp/g${N}d.cir > /dev/null; dn=$(echo "$(date +%s.%N) - $s" | bc)
  printf "N=%-4d unknowns=%-6d sparse=%.2fs  dense=%.2fs\n" $N $((N*N+2)) $sp $dn
done
