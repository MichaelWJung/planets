#pragma once

#include "simulator.h"
#include <raylib.h>

namespace planets {

// Owns the raylib window and renders Body state each frame in 3D.
// scale converts simulation metres to Raylib world units.
// Left-drag to orbit, scroll to zoom.
class Visualization {
public:
  Visualization(int width, int height, double scale);
  ~Visualization();

  [[nodiscard]] bool shouldClose() const;
  void render(const std::vector<Body> &bodies);

private:
  int    width_;
  int    height_;
  double scale_;

  Camera3D camera_{};
  float    azimuth_   = 45.0f;
  float    elevation_ = 20.0f;
  float    distance_  = 20.0f;

  // Maps mass to sphere radius in world units; auto-scales to the body set's mass range.
  [[nodiscard]] static float bodyRadius(double mass_kg, double log_min, double log_max);
};

} // namespace planets
