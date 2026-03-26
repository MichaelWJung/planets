#include "simulator.h"

#include <gtest/gtest.h>

using namespace planets;
using namespace mp_units;
using namespace mp_units::si;
using namespace mp_units::isq;

TEST(SimulatorTest, EmptySimulator) {
  Simulator sim{{}, 1.0 * second};
  sim.step();
  EXPECT_TRUE(sim.getState().empty());
}

TEST(SimulatorTest, SingleBodyAtRest) {
  Body body{
      .mass     = 1.0 * kilogram,
      .position = position_t{},
      .velocity = velocity_t{},
  };
  Simulator sim{{body}, 1.0 * second};
  sim.step();

  const auto &state = sim.getState();
  EXPECT_EQ(state[0].position, body.position);
  EXPECT_EQ(state[0].velocity, body.velocity);
}

TEST(SimulatorTest, SingleBodyWithVelocity) {
  // Body moving at 1 m/s along the x-axis, starting at the origin
  auto vel = velocity_t{cartesian_vector<double>{1.0, 0.0, 0.0}, isq::velocity[metre / second]};
  Body body{
      .mass     = 1.0 * kilogram,
      .position = position_t{},
      .velocity = vel,
  };
  Simulator sim{{body}, 1.0 * second};
  sim.step();

  // No forces on a single body, so velocity is unchanged.
  // Position = initial + velocity * dt = (0,0,0) + (1,0,0)*1s = (1,0,0) m
  auto expected_pos = position_t{cartesian_vector<double>{1.0, 0.0, 0.0}, isq::position_vector[metre]};
  const auto &state = sim.getState();
  EXPECT_EQ(state[0].position, expected_pos);
  EXPECT_EQ(state[0].velocity, vel);
}
