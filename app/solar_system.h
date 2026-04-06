#pragma once

#include "simulator.h"

namespace planets {

// Solar system bodies with approximate circular-orbit initial conditions.
// All planets placed on the +x axis; orbital velocity is in the +y direction.
// Masses and mean orbital radii from NASA planetary fact sheets.
inline std::vector<Body> solarSystemBodies() {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  const auto pos = [](double x, double y, double z) -> position_t {
    return position_t{cartesian_vector<double>{x, y, z} * isq::displacement[m],
                      solar_system_center_of_mass};
  };
  const auto vel = [](double vx, double vy, double vz) -> velocity_t {
    return cartesian_vector<double>{vx, vy, vz} * isq::velocity[m / s];
  };

  return {
      // Sun
      Body{.mass = 1.989e30 * kg, .position = pos(0, 0, 0), .velocity = vel(0, 0, 0)},
      // Mercury  r=0.387 AU  v_circ=47.87 km/s
      Body{.mass = 3.302e23 * kg, .position = pos(5.791e10, 0, 0), .velocity = vel(0, 47870, 0)},
      // Venus    r=0.723 AU  v_circ=35.02 km/s
      Body{.mass = 4.869e24 * kg, .position = pos(1.082e11, 0, 0), .velocity = vel(0, 35020, 0)},
      // Earth    r=1.000 AU  v_circ=29.78 km/s
      Body{.mass = 5.972e24 * kg, .position = pos(1.496e11, 0, 0), .velocity = vel(0, 29780, 0)},
      // Mars     r=1.524 AU  v_circ=24.13 km/s
      Body{.mass = 6.419e23 * kg, .position = pos(2.279e11, 0, 0), .velocity = vel(0, 24130, 0)},
      // Jupiter  r=5.203 AU  v_circ=13.06 km/s
      Body{.mass = 1.899e27 * kg, .position = pos(7.784e11, 0, 0), .velocity = vel(0, 13060, 0)},
      // Saturn   r=9.537 AU  v_circ= 9.68 km/s
      Body{.mass = 5.685e26 * kg, .position = pos(1.432e12, 0, 0), .velocity = vel(0, 9680, 0)},
      // Uranus   r=19.19 AU  v_circ= 6.80 km/s
      Body{.mass = 8.682e25 * kg, .position = pos(2.867e12, 0, 0), .velocity = vel(0, 6800, 0)},
      // Neptune  r=30.07 AU  v_circ= 5.43 km/s
      Body{.mass = 1.024e26 * kg, .position = pos(4.515e12, 0, 0), .velocity = vel(0, 5430, 0)},
  };
}

} // namespace planets
