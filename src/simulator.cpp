#include "simulator.h"

#include <cmath>

namespace planets {

namespace {

constexpr double G = 1;

void updateVelocity(Body &b, double dt) {
  b.velocity = b.velocity + b.acceleration * dt;
}

void updatePosition(Body &b, double dt) {
  b.position = b.position + b.velocity * dt;
}

void recomputeAccelerations(std::vector<Body> &bodies) {
  for (Body &b : bodies) {
    b.acceleration = Vec3{0, 0, 0};
  }
  for (Body &b1 : bodies) {
    for (Body &b2 : bodies) {
      if (&b1 != &b2) {
        const auto r = b2.position - b1.position;
        const double norm = computeNorm(r);
        const double cubeNorm = norm * norm * norm;
        const auto force = -G * b1.mass * b2.mass / cubeNorm * r;
        b1.acceleration = b1.acceleration + 1.0 / b1.mass * force;
        b2.acceleration = b2.acceleration - 1.0 / b2.mass * force;
      }
    }
  }
}

} // namespace

Vec3 operator*(const Vec3 &v, double a) {
  return {v[0] * a, v[1] * a, v[2] * a};
}

Vec3 operator*(double a, const Vec3 &v) {
  return {a * v[0], a * v[1], a * v[2]};
}

Vec3 operator+(const Vec3 &u, const Vec3 &v) {
  return {u[0] + v[0], u[1] + v[1], u[2] + v[2]};
}

Vec3 operator-(const Vec3 &u, const Vec3 &v) {
  return {u[0] - v[0], u[1] - v[1], u[2] - v[2]};
}

double computeNorm(const Vec3 &v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

Simulator::Simulator(std::vector<Body> bodies, double dt)
    : bodies_{std::move(bodies)}, dt_{dt} {}

const std::vector<Body> &Simulator::getState() const { return bodies_; }

void Simulator::step() {
  const auto dt2 = dt_ / 2.0;
  for (Body &b : bodies_) {
    updateVelocity(b, dt2);
  }
  for (Body &b : bodies_) {
    updatePosition(b, dt_);
  }
  recomputeAccelerations(bodies_);
  for (Body &b : bodies_) {
    updateVelocity(b, dt2);
  }
}

} // namespace planets
