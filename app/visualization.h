#pragma once

#include "simulator.h"
#include "trails.h"

namespace planets {

// Owns the raylib window and renders Body state each frame.
// metersPerPixel maps simulation distance (m) to screen pixels.
class Visualization {
public:
  Visualization(int width, int height, double meters_per_pixel);
  ~Visualization();

  [[nodiscard]] bool shouldClose() const;
  void render(const std::vector<Body> &bodies, const Trails &trails);

private:
  int width_;
  int height_;
  double meters_per_pixel_;

  // Maps body mass (kg) to a screen radius in pixels; minimum 1 px.
  // log_min/log_max are the log10 mass bounds of the current body set.
  [[nodiscard]] static int bodyRadius(double mass_kg, double log_min, double log_max);
};

} // namespace planets
