#include "plummer.h"
#include "simulator.h"
#include "trails.h"
#include "visualization.h"

int main() {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  // Cluster parameters: ~1e5 solar masses, scale radius ~10 pc
  constexpr auto M_total = 1.989e35 * isq::mass[kg];
  constexpr auto a       = 3.086e17 * isq::length[m];

  // Viewport: 1 scale radius ≈ 154 px; half-screen covers ~6 scale radii
  constexpr int    window_width     = 1920;
  constexpr int    window_height    = 1200;
  constexpr double meters_per_pixel = 2e15;

  // dt ≈ 0.0007 t_dyn; 5 steps/frame at 60 fps → ~300 kyr/s ≈ 0.2 t_dyn/s
  constexpr auto year           = 365 * d;
  constexpr auto dt             = 1000 * year;
  constexpr int  steps_per_frame = 5;

  // No warm-up needed: Plummer initial conditions are already in virial equilibrium
  planets::Simulator sim{planets::plummerBodies(1000, M_total, a), dt};
  planets::Trails    trails{100};  // short trail — 1000 bodies × 100 positions

  planets::Visualization vis{window_width, window_height, meters_per_pixel};
  while (!vis.shouldClose()) {
    for (int i = 0; i < steps_per_frame; ++i) {
      sim.step();
      trails.record(sim.getState());
    }
    vis.render(sim.getState(), trails);
  }
}
