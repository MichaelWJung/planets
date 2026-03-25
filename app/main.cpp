#include "simulator.h"

int main() {
  planets::Simulator sim{{}, 1.0};
  sim.step();
  return 0;
}
