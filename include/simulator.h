#pragma once

#include <array>
#include <vector>
#include <mp-units/math.h>
#include <mp-units/systems/si.h>

namespace planets {

template<typename T>
class Vec3 {
public:
  Vec3(T x, T y, T z) : coords_{{x, y, z}} {}

  T &operator[](std::ptrdiff_t i) { return coords_[i]; }
  const T &operator[](std::ptrdiff_t i) const { return coords_[i]; }

private:
  std::array<T, 3> coords_;
};

template<typename T>
Vec3<T> operator+(const Vec3<T> &u, const Vec3<T> &v) {
  return {u[0] + v[0], u[1] + v[1], u[2] + v[2]};
}

template<typename T>
Vec3<T> operator-(const Vec3<T> &u, const Vec3<T> &v) {
  return {u[0] - v[0], u[1] - v[1], u[2] - v[2]};
}

template<typename T, typename S>
auto operator*(const Vec3<T> &v, S a) {
  return Vec3<decltype(std::declval<T>() * std::declval<S>())>{v[0] * a, v[1] * a, v[2] * a};
}

template<typename S, typename T>
auto operator*(S a, const Vec3<T> &v) {
  return Vec3<decltype(std::declval<S>() * std::declval<T>())>{a * v[0], a * v[1], a * v[2]};
}

template<typename T>
auto computeNorm(const Vec3<T> &v) {
  return mp_units::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

using mass_t         = mp_units::quantity<mp_units::si::kilogram, double>;
using length_t       = mp_units::quantity<mp_units::si::metre, double>;
using time_t         = mp_units::quantity<mp_units::si::second, double>;
using velocity_t     = mp_units::quantity<mp_units::si::metre / mp_units::si::second, double>;
using acceleration_t = mp_units::quantity<mp_units::si::metre / mp_units::si::second / mp_units::si::second, double>;

struct Body {
  mass_t mass;
  Vec3<length_t> position;
  Vec3<velocity_t> velocity;
  Vec3<acceleration_t> acceleration;
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
