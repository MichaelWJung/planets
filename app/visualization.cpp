#include "visualization.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace planets {

using namespace mp_units::si::unit_symbols;

Visualization::Visualization(int width, int height, double scale)
    : width_(width), height_(height), scale_(scale) {
  InitWindow(width_, height_, "Planets");
  SetTargetFPS(60);
  camera_.target     = {0.0f, 0.0f, 0.0f};
  camera_.up         = {0.0f, 1.0f, 0.0f};
  camera_.fovy       = 45.0f;
  camera_.projection = CAMERA_PERSPECTIVE;
}

Visualization::~Visualization() { CloseWindow(); }

bool Visualization::shouldClose() const { return WindowShouldClose(); }

float Visualization::bodyRadius(double mass_kg, double log_min, double log_max) {
  constexpr float min_r = 0.02f;
  constexpr float max_r = 0.10f;
  const double t = (log_max > log_min)
      ? std::clamp((std::log10(mass_kg) - log_min) / (log_max - log_min), 0.0, 1.0)
      : 1.0;
  return min_r + static_cast<float>(t) * (max_r - min_r);
}

void Visualization::render(const std::vector<Body> &bodies) {
  // Orbital camera: left-drag rotates, scroll wheel zooms
  if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
    const Vector2 d = GetMouseDelta();
    azimuth_   -= d.x * 0.3f;
    elevation_ -= d.y * 0.3f;
    elevation_  = std::clamp(elevation_, -89.0f, 89.0f);
  }
  distance_ = std::max(1.0f, distance_ * (1.0f - 0.1f * GetMouseWheelMove()));

  const float az = azimuth_   * DEG2RAD;
  const float el = elevation_ * DEG2RAD;
  camera_.position = {distance_ * std::cos(el) * std::sin(az),
                      distance_ * std::sin(el),
                      distance_ * std::cos(el) * std::cos(az)};

  double log_mass_min = std::numeric_limits<double>::max();
  double log_mass_max = std::numeric_limits<double>::lowest();
  for (const auto &body : bodies) {
    const double log_m = std::log10(body.mass.numerical_value_in(kg));
    log_mass_min = std::min(log_mass_min, log_m);
    log_mass_max = std::max(log_mass_max, log_m);
  }

  BeginDrawing();
  ClearBackground(BLACK);
  BeginMode3D(camera_);

  for (const Body &body : bodies) {
    const auto &vec = body.position
        .quantity_ref_from(solar_system_center_of_mass)
        .numerical_value_ref_in(m);
    const Vector3 pos{static_cast<float>(vec[0] / scale_),
                      static_cast<float>(vec[1] / scale_),
                      static_cast<float>(vec[2] / scale_)};
    const float r = bodyRadius(body.mass.numerical_value_in(kg), log_mass_min, log_mass_max);
    DrawSphereEx(pos, r, 6, 6, WHITE);
  }

  EndMode3D();
  EndDrawing();
}

} // namespace planets
