#pragma once

#include <array>
#include <vector>

namespace planets {

class Vec3 {
public:
  Vec3(double x, double y, double z) : coords_{{x, y, z}} {}

  double &operator[](std::ptrdiff_t i) { return coords_[i]; }
  const double &operator[](std::ptrdiff_t i) const { return coords_[i]; }

private:
  std::array<double, 3> coords_;
};

Vec3 operator*(const Vec3 &v, double a);
Vec3 operator*(double a, const Vec3 &v);
Vec3 operator+(const Vec3 &u, const Vec3 &v);
Vec3 operator-(const Vec3 &u, const Vec3 &v);

double computeNorm(const Vec3 &v);

struct Body {
  double mass;
  Vec3 position;
  Vec3 velocity;
  Vec3 acceleration;
};

class Simulator {
public:
  Simulator(std::vector<Body> bodies, double dt);

  // get the current state of all bodies
  const std::vector<Body> &getState() const;

  // simulate one time step
  void step();

private:
  std::vector<Body> bodies_;
  double dt_;
};

} // namespace planets
