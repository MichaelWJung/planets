#pragma once

#include <mp-units/cartesian_vector.h>
#include <mp-units/systems/isq/space_and_time.h>
#include <mp-units/systems/si.h>
#include <vector>

namespace planets {

// Gravitational constant
inline constexpr auto G = [] {
  using namespace mp_units::si::unit_symbols;
  return 6.6743e-11 * (m3 / (kg * s2));
}();

inline constexpr struct solar_system_center_of_mass final
    : mp_units::absolute_point_origin<mp_units::isq::displacement> {
} solar_system_center_of_mass;

using mass_t =
    mp_units::quantity<mp_units::isq::mass[mp_units::si::kilogram], double>;
using time_t =
    mp_units::quantity<mp_units::isq::time[mp_units::si::second], double>;
using position_t =
    mp_units::quantity_point<mp_units::isq::displacement[mp_units::si::metre],
                             solar_system_center_of_mass,
                             mp_units::cartesian_vector<double>>;
using displacement_t =
    mp_units::quantity<mp_units::isq::displacement[mp_units::si::metre],
                       mp_units::cartesian_vector<double>>;
using velocity_t = mp_units::quantity<
    mp_units::isq::velocity[mp_units::si::metre / mp_units::si::second],
    mp_units::cartesian_vector<double>>;
using acceleration_t = mp_units::quantity<
    mp_units::isq::acceleration[mp_units::si::metre / mp_units::si::second /
                                mp_units::si::second],
    mp_units::cartesian_vector<double>>;

struct Body {
  mass_t mass;
  position_t position;
  velocity_t velocity;
};

class Simulator {
public:
  Simulator(std::vector<Body> bodies, time_t dt);

  // get the current state of all bodies
  const std::vector<Body> &getState() const;

  // simulate one time step
  void step();

private:
  std::vector<Body> bodies_;
  std::vector<acceleration_t> accelerations_;
  time_t dt_;
};

} // namespace planets
