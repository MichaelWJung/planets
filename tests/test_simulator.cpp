#include "simulator.h"

#include <gtest/gtest.h>

using namespace mp_units;
using namespace mp_units::si::unit_symbols;

namespace planets {
namespace {

TEST(SimulatorTest, EmptySimulator) {
  Simulator sim{{}, 1 * s};
  sim.step();
  EXPECT_TRUE(sim.getState().empty());
}

TEST(SimulatorTest, SingleBodyAtRest) {
  Body body{
      .mass = 1 * kg,
      .position = position_t{},
      .velocity = velocity_t{},
  };
  Simulator sim{{body}, 1 * s};
  sim.step();

  const auto &state = sim.getState();
  EXPECT_EQ(state[0].position, body.position);
  EXPECT_EQ(state[0].velocity, body.velocity);
}

TEST(SimulatorTest, SingleBodyWithVelocity) {
  // Body moving at 1 m/s along the x-axis, starting at the origin
  auto vel =
      velocity_t{cartesian_vector<double>{1.0, 0.0, 0.0}, isq::velocity[m / s]};
  Body body{
      .mass = 1 * kg,
      .position = position_t{},
      .velocity = vel,
  };
  Simulator sim{{body}, 1 * s};
  sim.step();

  // No forces on a single body, so velocity is unchanged.
  // Position = initial + velocity * dt = (0,0,0) + (1,0,0)*1s = (1,0,0) m
  auto expected_pos =
      position_t{cartesian_vector<double>{1.0, 0.0, 0.0} * isq::displacement[m],
                 solar_system_center_of_mass};
  const auto &state = sim.getState();
  EXPECT_EQ(state[0].position, expected_pos);
  EXPECT_EQ(state[0].velocity, vel);
}

TEST(SimulatorTest, TwoBodiesAttractEachOther) {
  // Body 0 (mass m) at (-x, 0, 0), body 1 (mass 2m) at (+x, 0, 0), both at
  // rest. After a step they should move toward each other. Same gravitational
  // force on both; body 0 has half the mass → exactly twice the acceleration.
  constexpr auto mass = 1.0e24 * isq::mass[kg];
  constexpr auto x = 1.0e10 * isq::length[m];
  constexpr auto dt = 1.0 * isq::time[d];

  Body body0{
      .mass = mass,
      .position = position_t{cartesian_vector<double>{-x.numerical_value_in(m),
                                                      0.0, 0.0} *
                                 isq::displacement[m],
                             solar_system_center_of_mass},
      .velocity = velocity_t{},
  };
  Body body1{
      .mass = 2.0 * mass,
      .position = position_t{cartesian_vector<double>{x.numerical_value_in(m),
                                                      0.0, 0.0} *
                                 isq::displacement[m],
                             solar_system_center_of_mass},
      .velocity = velocity_t{},
  };

  Simulator sim{{body0, body1}, dt};
  sim.step();

  const auto state = sim.getState();

  // Same force on both; body 0 has half the mass → exactly twice the velocity
  // change. Checks all three components: y and z stay zero, x is in the correct
  // -2:1 ratio.
  EXPECT_EQ(state[0].velocity, -2.0 * state[1].velocity);

  // Gravity is attractive: body 0 (at -x) gains a positive x velocity toward
  // body 1 (at +x).
  EXPECT_GT(state[0].velocity.numerical_value_in(m / s)[0], 0.0);

  // After another step, positions change (velocity is now nonzero).
  sim.step();
  const auto state2 = sim.getState();

  // Displacement ratio is exactly -2:1, mirroring the velocity ratio (momentum
  // conservation). Also implicitly checks that y and z position components stay
  // zero.
  EXPECT_EQ(state2[0].position - body0.position,
            -2.0 * (state2[1].position - body1.position));

  // TODO: check distance actually decreased
  // Separation changed — it decreased because body 0 moved toward body 1.
  EXPECT_NE(state2[1].position - state2[0].position,
            body1.position - body0.position);
}

} // namespace
} // namespace planets
