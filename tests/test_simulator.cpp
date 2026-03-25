#include "simulator.h"

int main() {
  planets::Simulator sim{{}, 1.0 * mp_units::si::second};
  sim.step();
  return 0;
}
