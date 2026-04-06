#include "simulator.h"
#include "solar_system.h"
#include "visualization.h"

int main() {
  using namespace mp_units::si::unit_symbols;

  // Viewport: zoom so Mercury (~0.39 AU) clears the Sun's drawn radius.
  // At this scale Jupiter (~5.2 AU) fits comfortably; Saturn (~9.5 AU) is
  // visible on the horizontal axis; Uranus and Neptune are off-screen.
  constexpr int window_width = 1920;
  constexpr int window_height = 1200;
  constexpr double meters_per_pixel = 1.6e9; // ~1 AU ≈ 94 px

  // 5 simulation days per rendered frame → ~300 sim-days/s at 60 fps.
  // Earth completes one orbit in ~1.2 s of real time.
  constexpr int steps_per_frame = 5;

  planets::Simulator sim{planets::solarSystemBodies(), 1.0 * d};
  planets::Visualization vis{window_width, window_height, meters_per_pixel};

  while (!vis.shouldClose()) {
    for (int i = 0; i < steps_per_frame; ++i) {
      sim.step();
    }
    vis.render(sim.getState());
  }
}
