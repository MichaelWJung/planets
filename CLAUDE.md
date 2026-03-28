# CLAUDE.md

## Project

Astronomical N-body simulation using Newtonian gravity (C++23). Primary scenario: our solar system. Potentially other scenarios later.

## Build

```bash
# Configure (once)
cd build && cmake -G Ninja -DCMAKE_TOOLCHAIN_FILE=clang21-toolchain.cmake ..

# Incremental build
ninja -C build

# Run tests
ninja -C build && ./build/tests/test_simulator
```

## Current Stage: Correctness & Clean Architecture

We are in stage 1 of a deliberate progression:
1. **Now:** Correctness, clean code, good architecture
2. Scale up body count until performance degrades
3. Profile and optimize (CPU)
4. GPU programming

**Do not skip stages.** Do not introduce optimizations (SIMD, parallelism, cache tricks, algorithmic shortcuts like Barnes-Hut) until the user has profiled, hit real limits, and decided to move to the next stage. When in doubt, prefer the clearer implementation over the faster one.

## Code Conventions

- C++23, no extensions
- Use mp-units for all physical quantities — never raw floats for physics values
- Use ISQ quantities and SI units (as already established in the codebase)
- Warnings are strict; code must compile cleanly
