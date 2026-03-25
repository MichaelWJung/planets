#pragma once

#include <mp-units/cartesian_vector.h>
#include <mp-units/systems/si.h>
#include <vector>

namespace planets {

using mass_t =
    mp_units::quantity<mp_units::isq::mass[mp_units::si::kilogram], double>;
using time_t =
    mp_units::quantity<mp_units::isq::time[mp_units::si::second], double>;
using position_t =
    mp_units::quantity<mp_units::isq::position_vector[mp_units::si::metre],
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
  acceleration_t acceleration;
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
  time_t dt_;
};

} // namespace planets
