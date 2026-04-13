#include "visualization.h"

#include <raylib.h>

#include <algorithm>
#include <cmath>

namespace planets {

using namespace mp_units::si::unit_symbols;

Visualization::Visualization(int width, int height, double meters_per_pixel)
    : width_(width), height_(height), meters_per_pixel_(meters_per_pixel) {
  InitWindow(width_, height_, "Planets");
  SetTargetFPS(60);
}

Visualization::~Visualization() { CloseWindow(); }

bool Visualization::shouldClose() const { return WindowShouldClose(); }

int Visualization::bodyRadius(double mass_kg) {
  // Log-linear mapping from mass to pixel radius.
  // Calibrated for solar system bodies: Mercury (~1e23 kg) → 1 px,
  // Sun (~2e30 kg) → 6 px.
  constexpr double log_min = 23.0;
  constexpr double log_max = 30.3;
  constexpr double min_px = 1.0;
  constexpr double max_px = 6.0;
  const double t =
      std::clamp((std::log10(mass_kg) - log_min) / (log_max - log_min), 0.0, 1.0);
  return static_cast<int>(std::max(1.0, min_px + t * (max_px - min_px)));
}

void Visualization::render(const std::vector<Body> &bodies,
                           const Trails &trails) {
  const auto toScreen = [&](const position_t &pos) -> std::pair<int, int> {
    const auto &vec =
        pos.quantity_ref_from(solar_system_center_of_mass).numerical_value_ref_in(m);
    return {static_cast<int>(width_ / 2 + vec[0] / meters_per_pixel_),
            static_cast<int>(height_ / 2 - vec[1] / meters_per_pixel_)};
  };

  BeginDrawing();
  ClearBackground(BLACK);

  for (const auto &trail : trails.data()) {
    const auto n = static_cast<double>(trail.size());
    for (std::size_t i = 0; i < trail.size(); ++i) {
      const double t = (n > 1) ? static_cast<double>(i) / (n - 1) : 1.0;
      const auto gray = static_cast<unsigned char>(t * 120);
      const auto [sx, sy] = toScreen(trail[i]);
      DrawPixel(sx, sy, Color{gray, gray, gray, 255});
    }
  }

  for (const Body &body : bodies) {
    const auto [sx, sy] = toScreen(body.position);
    const int r = bodyRadius(body.mass.numerical_value_in(kg));
    DrawCircle(sx, sy, static_cast<float>(r), WHITE);
  }

  EndDrawing();
}

} // namespace planets
