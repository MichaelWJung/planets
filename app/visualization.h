#pragma once

#include "simulator.h"

namespace planets {

// Owns the raylib window and renders Body state each frame.
// metersPerPixel maps simulation distance (m) to screen pixels.
class Visualization {
public:
  Visualization(int width, int height, double meters_per_pixel);
  ~Visualization();

  [[nodiscard]] bool shouldClose() const;
  void render(const std::vector<Body> &bodies);

private:
  int width_;
  int height_;
  double meters_per_pixel_;

  // Maps body mass (kg) to a screen radius in pixels; minimum 1 px.
  [[nodiscard]] static int bodyRadius(double mass_kg);
};

} // namespace planets
