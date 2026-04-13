#include "simulator.h"
#include "solar_system.h"
#include "trails.h"
#include "visualization.h"

int main() {
  using namespace mp_units::si::unit_symbols;

  // Viewport: all eight planets fit within the window.
  // Neptune (~30 AU) sits ~562 px from center; 1 AU ≈ 19 px.
  constexpr int window_width = 1920;
  constexpr int window_height = 1200;
  constexpr double meters_per_pixel = 8e9; // ~1 AU ≈ 19 px

  constexpr auto year = 365 * d;
  constexpr auto initial_delay = 10'000 * year;

  // 5 simulation days per rendered frame → ~300 sim-days/s at 60 fps.
  // Earth completes one orbit in ~1.2 s of real time.
  constexpr int steps_per_frame = 5;

  planets::Simulator sim{planets::solarSystemBodies(), 1.0 * d};
  planets::Trails trails{10000};

  for (auto i = 1 * d; i < initial_delay; i += 1 * d) {
    sim.step();
    trails.record(sim.getState());
  }

  planets::Visualization vis{window_width, window_height, meters_per_pixel};
  while (!vis.shouldClose()) {
    for (int i = 0; i < steps_per_frame; ++i) {
      sim.step();
      trails.record(sim.getState());
    }
    vis.render(sim.getState(), trails);
  }
}
