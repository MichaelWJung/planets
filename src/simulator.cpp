#include "simulator.h"

#include <ranges>

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

void updateVelocity(Body &b, const acceleration_t &accel, time_t dt) {
  b.velocity = quantity_cast<isq::velocity>(b.velocity + accel * dt);
}

void updatePosition(Body &b, time_t dt) {
  b.position =
      quantity_cast<isq::position_vector>(b.position + b.velocity * dt);
}

void recomputeAccelerations(const std::vector<Body> &bodies,
                            std::vector<acceleration_t> &accelerations) {
  for (acceleration_t &a : accelerations) {
    a = acceleration_t{};
  }
  for (auto [body_i, accel_i] : std::views::zip(bodies, accelerations)) {
    for (auto [body_j, accel_j] : std::views::zip(bodies, accelerations)) {
      if (&body_i != &body_j) {
        const auto r = body_j.position - body_i.position;
        const auto norm = computeLength(r);
        const auto cubeNorm = norm * norm * norm;
        const auto force = G * body_i.mass * body_j.mass / cubeNorm * r;
        accel_i = accel_i + quantity_cast<isq::acceleration>(1.0 / body_i.mass * force);
        accel_j = accel_j - quantity_cast<isq::acceleration>(1.0 / body_j.mass * force);
      }
    }
  }
}

} // namespace

Simulator::Simulator(std::vector<Body> bodies, time_t dt)
    : bodies_{std::move(bodies)}, accelerations_(bodies_.size()), dt_{dt} {}

const std::vector<Body> &Simulator::getState() const { return bodies_; }

void Simulator::step() {
  const auto dt2 = dt_ / 2.0;
  for (auto [body, accel] : std::views::zip(bodies_, accelerations_)) {
    updateVelocity(body, accel, dt2);
  }
  for (Body &b : bodies_) {
    updatePosition(b, dt_);
  }
  recomputeAccelerations(bodies_, accelerations_);
  for (auto [body, accel] : std::views::zip(bodies_, accelerations_)) {
    updateVelocity(body, accel, dt2);
  }
}

} // namespace planets
