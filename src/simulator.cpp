#include "simulator.h"

namespace planets {

namespace {

using namespace mp_units;
using namespace mp_units::si;

// Gravitational constant in SI units: m^3 kg^-1 s^-2
constexpr auto G =
    6.674e-11 * (metre * metre * metre / (kilogram * second * second));

auto computeLength(auto r) {
  return r.numerical_value_in(metre).magnitude() * isq::length[metre];
}

void updateVelocity(Body &b, time_t dt) {
  b.velocity = quantity_cast<isq::velocity>(b.velocity + b.acceleration * dt);
}

void updatePosition(Body &b, time_t dt) {
  b.position =
      quantity_cast<isq::position_vector>(b.position + b.velocity * dt);
}

void recomputeAccelerations(std::vector<Body> &bodies) {
  for (Body &b : bodies) {
    b.acceleration = acceleration_t{};
  }
  for (Body &b1 : bodies) {
    for (Body &b2 : bodies) {
      if (&b1 != &b2) {
        const auto r = b2.position - b1.position;
        const auto norm = computeLength(r);
        const auto cubeNorm = norm * norm * norm;
        const auto force = -G * b1.mass * b2.mass / cubeNorm * r;
        b1.acceleration = b1.acceleration + quantity_cast<isq::acceleration>(1.0 / b1.mass * force);
        b2.acceleration = b2.acceleration - quantity_cast<isq::acceleration>(1.0 / b2.mass * force);
      }
    }
  }
}

} // namespace

Simulator::Simulator(std::vector<Body> bodies, time_t dt)
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
