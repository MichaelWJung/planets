#include "simulator.h"

#include <ranges>

namespace planets {

namespace {

using namespace mp_units;
using namespace mp_units::si::unit_symbols;

auto computeLength(QuantityOf<isq::displacement> auto d) {
  const auto &r = d.numerical_value_ref_in(d.unit);
  return quantity_cast<isq::length>(r.norm() * d.unit);
}

void updateVelocity(Body &b, const acceleration_t &accel, time_t dt) {
  b.velocity += accel * dt;
}

void updatePosition(Body &b, time_t dt) { b.position += b.velocity * dt; }

void recomputeAccelerations(const std::vector<Body> &bodies,
                            std::vector<acceleration_t> &accelerations) {
  for (acceleration_t &a : accelerations) {
    a = {};
  }
  for (size_t i = 0; i < bodies.size(); ++i) {
    for (size_t j = i + 1; j < bodies.size(); ++j) {
      const auto r = bodies[j].position - bodies[i].position;
      const auto norm = computeLength(r);
      const auto cubeNorm = norm * norm * norm;
      const auto force = -G * bodies[i].mass * bodies[j].mass / cubeNorm * r;
      accelerations[i] -= 1.0 / bodies[i].mass * force;
      accelerations[j] += 1.0 / bodies[j].mass * force;
    }
  }
}

} // namespace

Simulator::Simulator(std::vector<Body> bodies, time_t dt)
    : bodies_{std::move(bodies)}, accelerations_(bodies_.size()), dt_{dt} {
  recomputeAccelerations(bodies_, accelerations_);
}

const std::vector<Body> &Simulator::getState() const { return bodies_; }

void Simulator::step() {
  const auto dt_half = dt_ / 2.0;
  for (auto [body, accel] : std::views::zip(bodies_, accelerations_)) {
    updateVelocity(body, accel, dt_half);
  }
  for (Body &b : bodies_) {
    updatePosition(b, dt_);
  }
  recomputeAccelerations(bodies_, accelerations_);
  for (auto [body, accel] : std::views::zip(bodies_, accelerations_)) {
    updateVelocity(body, accel, dt_half);
  }
}

} // namespace planets
