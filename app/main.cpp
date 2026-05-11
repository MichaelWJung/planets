#include "plummer.h"
#include "simulator.h"
#include "visualization.h"

int main() {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  // Cluster: ~1e5 solar masses, scale radius ~10 pc
  constexpr auto M_total = 1.989e35 * isq::mass[kg];
  constexpr auto a       = 3.086e17 * isq::length[m];

  // 1 world unit = 1 scale radius; cluster spans ~±5 world units
  const double scale = a.numerical_value_in(m);

  // dt ≈ 0.0007 t_dyn; 5 steps/frame at 60 fps → ~300 kyr/s ≈ 0.2 t_dyn/s
  constexpr auto year            = 365 * d;
  constexpr auto dt              = 1000 * year;
  constexpr int  steps_per_frame = 5;

  planets::Simulator sim{planets::plummerBodies(1000, M_total, a), dt};
  planets::Visualization vis{1920, 1200, scale};

  while (!vis.shouldClose()) {
    for (int i = 0; i < steps_per_frame; ++i)
      sim.step();
    vis.render(sim.getState());
  }
}
