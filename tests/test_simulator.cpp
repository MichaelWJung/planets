#include "simulator.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mp-units/math.h>
#include <ranges>

using namespace mp_units;
using namespace mp_units::si::unit_symbols;

namespace planets {
namespace {

using testing::AllOf;
using testing::ElementsAre;
using testing::Field;
using testing::IsEmpty;

[[nodiscard]] displacement_t make_displacement(QuantityOf<isq::length> auto x,
                                               QuantityOf<isq::length> auto y,
                                               QuantityOf<isq::length> auto z) {
  return cartesian_vector<double>{x.numerical_value_in(m),
                                  y.numerical_value_in(m),
                                  z.numerical_value_in(m)} *
         isq::displacement[m];
}

[[nodiscard]] position_t make_position(QuantityOf<isq::length> auto x,
                                       QuantityOf<isq::length> auto y,
                                       QuantityOf<isq::length> auto z) {
  return position_t{make_displacement(x, y, z), solar_system_center_of_mass};
}

[[nodiscard]] velocity_t make_velocity(auto vx, auto vy, auto vz)
  requires requires {
    vx.numerical_value_in(m / s);
    vy.numerical_value_in(m / s);
    vz.numerical_value_in(m / s);
  }
{
  return cartesian_vector<double>{vx.numerical_value_in(m / s),
                                  vy.numerical_value_in(m / s),
                                  vz.numerical_value_in(m / s)} *
         isq::velocity[m / s];
}

TEST(SimulatorTest, EmptySimulator) {
  Simulator sim{{}, 1 * s};
  sim.step();
  EXPECT_THAT(sim.getState(), IsEmpty());
}

TEST(SimulatorTest, SingleBodyAtRest) {
  const Body body{
      .mass = 1 * kg,
      .position = position_t{},
      .velocity = velocity_t{},
  };
  Simulator sim{{body}, 1 * s};
  sim.step();

  EXPECT_THAT(sim.getState(),
              ElementsAre(AllOf(Field(&Body::position, body.position),
                                Field(&Body::velocity, body.velocity))));
}

TEST(SimulatorTest, SingleBodyWithVelocity) {
  // Body moving at 1 m/s along the x-axis, starting at the origin
  const auto vel = make_velocity(1 * (m / s), 0 * (m / s), 0 * (m / s));
  const Body body{
      .mass = 1 * kg,
      .position = position_t{},
      .velocity = vel,
  };
  Simulator sim{{body}, 1 * s};
  sim.step();

  // No forces on a single body, so velocity is unchanged.
  // Position = initial + velocity * dt = (0,0,0) + (1,0,0)*1s = (1,0,0) m
  const auto expected_pos = make_position(1 * m, 0 * m, 0 * m);
  EXPECT_THAT(sim.getState(),
              ElementsAre(AllOf(Field(&Body::position, expected_pos),
                                Field(&Body::velocity, vel))));
}

TEST(SimulatorTest, TwoBodiesAttractEachOther) {
  // Body 0 (mass m) at (-x, 0, 0), body 1 (mass 2m) at (+x, 0, 0), both at
  // rest. After a step they should move toward each other. Same gravitational
  // force on both; body 0 has half the mass → exactly twice the acceleration.
  constexpr auto mass = 1.0e24 * isq::mass[kg];
  constexpr auto x = 1.0e10 * isq::length[m];
  constexpr auto dt = 1.0 * isq::time[d];

  const Body body0{
      .mass = mass,
      .position = make_position(-x, 0 * m, 0 * m),
      .velocity = velocity_t{},
  };
  const Body body1{
      .mass = 2.0 * mass,
      .position = make_position(x, 0 * m, 0 * m),
      .velocity = velocity_t{},
  };

  Simulator sim{{body0, body1}, dt};
  sim.step();

  const auto state1 = sim.getState();

  // Same force on both; body 0 has half the mass → exactly twice the velocity
  // change. Checks all three components: y and z stay zero, x is in the correct
  // -2:1 ratio.
  EXPECT_EQ(state1[0].velocity, -2.0 * state1[1].velocity);

  // Gravity is attractive: body 0 (at -x) gains a positive x velocity toward
  // body 1 (at +x).
  EXPECT_GT(state1[0].velocity.numerical_value_in(m / s)[0], 0.0);

  // After another step, positions change (velocity is now nonzero).
  sim.step();
  const auto state2 = sim.getState();

  // Displacement ratio is -2:1, mirroring the velocity ratio (momentum
  // conservation). Also implicitly checks that y and z position components stay
  // near zero.
  {
    const auto disp0 = state2[0].position - body0.position;
    const auto disp1 = state2[1].position - body1.position;
    constexpr double tol = 1e-5;
    for (int k : std::views::iota(0, 3))
      EXPECT_NEAR(disp0.numerical_value_in(m)[k],
                  -2.0 * disp1.numerical_value_in(m)[k], tol);
  }

  const auto initial_sep = norm(body1.position - body0.position);
  const auto final_sep = norm(state2[1].position - state2[0].position);
  EXPECT_LE(final_sep, initial_sep);
}

TEST(SimulatorTest, TwoBodiesCircularOrbit) {
  // Planet and moon in a mutual orbit about their common center of mass.
  // Planet: (-r_p, 0, 0) with velocity (0, -v_p, 0)
  // Moon:   (+r_m, 0, 0) with velocity (0, +v_m, 0)
  constexpr auto mass_planet = 6.0e24 * isq::mass[kg];
  constexpr auto mass_moon = 7.0e22 * isq::mass[kg];
  constexpr auto distance = 1.0e10 * isq::length[m];
  constexpr auto mass_total = mass_planet + mass_moon;
  const auto r_planet = mass_moon / mass_total * distance;
  const auto r_moon = mass_planet / mass_total * distance;
  const auto v_planet =
      mass_moon / mass_total * sqrt(G * mass_total / distance);
  const auto v_moon =
      mass_planet / mass_total * sqrt(G * mass_total / distance);
  const auto orbital_period =
      2 * std::numbers::pi *
      sqrt(distance * distance * distance / G / mass_total);

  constexpr auto dt = 1.0 * isq::time[d];

  const Body planet{
      .mass = mass_planet,
      .position = make_position(-r_planet, 0 * m, 0 * m),
      .velocity = make_velocity(0 * (m / s), -v_planet, 0 * (m / s)),
  };
  const Body moon{
      .mass = mass_moon,
      .position = make_position(r_moon, 0 * m, 0 * m),
      .velocity = make_velocity(0 * (m / s), v_moon, 0 * (m / s)),
  };

  // --- Distance conservation over 5000 steps ---
  {
    const auto initial_dist = norm(moon.position - planet.position);

    Simulator sim{{planet, moon}, dt};
    // Leapfrog is symplectic: distance error stays bounded (no secular drift).
    // With dt ≈ T/3600, induced eccentricity ≈ (dt/T)² ≈ 1e-7; 1e-4 gives ample
    // margin.
    constexpr double dist_tolerance_rel = 1e-4;
    for (int step = 0; step < 5000; ++step) {
      sim.step();
      const auto &state = sim.getState();
      const auto dist = norm(state[1].position - state[0].position);
      EXPECT_LE(abs(dist - initial_dist), dist_tolerance_rel * initial_dist)
          << "separation changed at step " << step;
    }
  }

  // --- Sign-flip checks ---
  // x flips sign at T/4 (and again at 3T/4); y flips sign at T/2.
  // We verify each flip by checking the sign just before and just after.
  {
    const int quarter_period_steps =
        static_cast<int>(floor<one>(orbital_period / 4.0 / dt));
    const int half_period_steps =
        static_cast<int>(floor<one>(orbital_period / 2.0 / dt));

    Simulator sim{{planet, moon}, dt};

    // x-component flips sign at T/4: positive just before, negative just after.
    for (int step = 0; step < quarter_period_steps; ++step)
      sim.step();
    {
      const auto sep = sim.getState()[1].position - sim.getState()[0].position;
      EXPECT_GT(sep.numerical_value_in(m)[0], 0.0)
          << "x-component should still be positive just before T/4";
    }
    sim.step();
    {
      const auto sep = sim.getState()[1].position - sim.getState()[0].position;
      EXPECT_LT(sep.numerical_value_in(m)[0], 0.0)
          << "x-component should be negative just after T/4";
    }

    // y-component flips sign at T/2: positive just before, negative just after.
    for (int step = quarter_period_steps + 1; step < half_period_steps; ++step)
      sim.step();
    {
      const auto sep = sim.getState()[1].position - sim.getState()[0].position;
      EXPECT_GT(sep.numerical_value_in(m)[1], 0.0)
          << "y-component should still be positive just before T/2";
    }
    sim.step();
    {
      const auto sep = sim.getState()[1].position - sim.getState()[0].position;
      EXPECT_LT(sep.numerical_value_in(m)[1], 0.0)
          << "y-component should be negative just after T/2";
    }
  }
}

} // namespace
} // namespace planets
