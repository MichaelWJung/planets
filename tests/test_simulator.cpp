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

TEST(SimulatorTest, TwoBodiesAttractEachOther) {
  // Body 0 (mass m) at (-d, 0, 0), body 1 (mass 2m) at (+d, 0, 0), both at rest.
  // After a step they should move toward each other.
  // Same gravitational force on both; body 0 has half the mass → exactly twice the acceleration.
  constexpr auto m  = 1.0e24 * isq::mass[kilogram];
  constexpr auto d  = 1.0e10 * isq::length[metre];
  constexpr auto dt = 1.0 * isq::time[day];

  Body body0{
      .mass     = m,
      .position = position_t{cartesian_vector<double>{-d.numerical_value_in(metre), 0.0, 0.0},
                              isq::position_vector[metre]},
      .velocity = velocity_t{},
  };
  Body body1{
      .mass     = 2.0 * m,
      .position = position_t{cartesian_vector<double>{d.numerical_value_in(metre), 0.0, 0.0},
                              isq::position_vector[metre]},
      .velocity = velocity_t{},
  };

  Simulator sim{{body0, body1}, dt};
  sim.step();

  const auto &state = sim.getState();

  // Same force on both; body 0 has half the mass → exactly twice the velocity change.
  // Checks all three components: y and z stay zero, x is in the correct -2:1 ratio.
  EXPECT_EQ(state[0].velocity, -2.0 * state[1].velocity);

  // Gravity is attractive: body 0 (at -x) gains a positive x velocity toward body 1 (at +x).
  // A sign check on one component is unavoidable here; mp-units has no ordering on vector quantities.
  EXPECT_GT(state[0].velocity.numerical_value_in(metre / second)[0], 0.0);

  // After another step, positions change (velocity is now nonzero).
  sim.step();
  const auto &state2 = sim.getState();

  // Displacement ratio is exactly -2:1, mirroring the velocity ratio (momentum conservation).
  // Also implicitly checks that y and z position components stay zero.
  EXPECT_EQ(state2[0].position - body0.position, -2.0 * (state2[1].position - body1.position));

  // Separation changed — it decreased because body 0 moved toward body 1.
  EXPECT_NE(state2[1].position - state2[0].position, body1.position - body0.position);
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
