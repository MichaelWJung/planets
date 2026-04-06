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
  // Calibrated for solar system bodies: Mercury (~1e23 kg) → ~2 px,
  // Sun (~2e30 kg) → ~25 px.
  constexpr double log_min = 23.0;
  constexpr double log_max = 30.3;
  constexpr double min_px = 2.0;
  constexpr double max_px = 25.0;
  const double t =
      std::clamp((std::log10(mass_kg) - log_min) / (log_max - log_min), 0.0, 1.0);
  return static_cast<int>(std::max(1.0, min_px + t * (max_px - min_px)));
}

void Visualization::render(const std::vector<Body> &bodies) {
  BeginDrawing();
  ClearBackground(BLACK);

  for (const Body &body : bodies) {
    const auto &disp =
        body.position.quantity_ref_from(solar_system_center_of_mass);
    const auto &vec = disp.numerical_value_ref_in(m);
    const auto sx =
        static_cast<float>(width_ / 2 + vec[0] / meters_per_pixel_);
    const auto sy =
        static_cast<float>(height_ / 2 - vec[1] / meters_per_pixel_);
    const int r =
        bodyRadius(body.mass.numerical_value_in(kg));
    DrawCircle(static_cast<int>(sx), static_cast<int>(sy),
               static_cast<float>(r), WHITE);
  }

  EndDrawing();
}

} // namespace planets
