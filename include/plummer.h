#pragma once

#include "simulator.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <ranges>

namespace planets {

// Generate N bodies distributed according to the Plummer sphere model.
// Masses are drawn from the Kroupa (2001) IMF and rescaled so their sum equals
// M_total. Positions and velocities follow the exact Plummer distribution
// function, placing the system in virial equilibrium. The returned system is
// in the center-of-mass frame (zero total momentum, center of mass at origin).
inline std::vector<Body> plummerBodies(
    int N,
    mass_t M_total,
    mp_units::quantity<mp_units::isq::length[mp_units::si::metre], double> a,
    std::uint64_t seed = 42)
{
    using namespace mp_units;
    using namespace mp_units::si::unit_symbols;

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // --- Kroupa (2001) IMF: piecewise power law ξ(m) ∝ m^{-α} ---
    //   Segment 1: α=1.3 on [0.08, 0.5]  M_sun  (abundant low-mass stars)
    //   Segment 2: α=2.3 on [0.5,  150]  M_sun  (solar-type and massive)
    // Continuity enforced at the break: k·m_break^{-2.3} = m_break^{-1.3}
    //   → k = m_break^{α2-α1} = m_break^1
    constexpr double m_lo    = 0.08;
    constexpr double m_break = 0.5;
    constexpr double m_hi    = 150.0;
    constexpr double beta1   = 1.0 - 1.3;   // -0.3  (exponent in antiderivative)
    constexpr double beta2   = 1.0 - 2.3;   // -1.3
    constexpr double k       = m_break;      // continuity factor

    // Segment weights = ∫ξ(m)dm over each segment
    const double W1 = (std::pow(m_break, beta1) - std::pow(m_lo,    beta1)) / beta1;
    const double W2 = k * (std::pow(m_hi, beta2) - std::pow(m_break, beta2)) / beta2;
    const double p1 = W1 / (W1 + W2);

    // Exact inverse-CDF sampler for ∝ m^{-α} on [m_a, m_b]
    const auto sample_segment = [&](double m_a, double m_b, double beta) {
        const double u = uniform(rng);
        return std::pow(u * (std::pow(m_b, beta) - std::pow(m_a, beta))
                            + std::pow(m_a, beta),
                        1.0 / beta);
    };

    std::vector<double> masses(N);
    for (auto& mass : masses)
        mass = (uniform(rng) < p1)
                 ? sample_segment(m_lo, m_break, beta1)
                 : sample_segment(m_break, m_hi, beta2);

    // Rescale raw samples so their sum equals M_total
    const double M_kg       = M_total.numerical_value_in(kg);
    const double raw_sum    = std::ranges::fold_left(masses, 0.0, std::plus{});
    const double mass_scale = M_kg / raw_sum;
    for (auto& mass : masses) mass *= mass_scale;   // now in kg

    // --- Plummer positions: inverse CDF of cumulative mass M(r)/M ---
    // M(r)/M = r^3/(r^2+a^2)^{3/2}  →  r = a / sqrt(X^{-2/3} - 1)
    //
    // --- Plummer velocities: rejection-sample speed from the DF ---
    // f(E) ∝ (-E)^{7/2}  →  p(q|r) ∝ q^2(1-q^2)^{7/2},  q = v/v_esc
    // Maximum of q^2(1-q^2)^{7/2} is at q^2=2/9
    const double a_m   = a.numerical_value_in(m);
    const double GM    = G.numerical_value_in(m3 / (kg * s2)) * M_kg;
    const double g_max = (2.0 / 9.0) * std::pow(7.0 / 9.0, 3.5);

    std::vector<Body> bodies;
    bodies.reserve(N);

    for (int i = 0; i < N; ++i) {
        // Position radius
        const double X   = uniform(rng);
        const double r   = a_m / std::sqrt(std::pow(X, -2.0 / 3.0) - 1.0);

        // Uniform direction on sphere (position)
        const double cos_t = 2.0 * uniform(rng) - 1.0;
        const double sin_t = std::sqrt(1.0 - cos_t * cos_t);
        const double phi   = 2.0 * std::numbers::pi * uniform(rng);

        // Speed via rejection sampling
        const double v_esc = std::sqrt(2.0 * GM / std::sqrt(r * r + a_m * a_m));
        double q;
        do {
            q = uniform(rng);
        } while (uniform(rng) * g_max > q * q * std::pow(1.0 - q * q, 3.5));
        const double v = q * v_esc;

        // Uniform direction on sphere (velocity)
        const double cos_p = 2.0 * uniform(rng) - 1.0;
        const double sin_p = std::sqrt(1.0 - cos_p * cos_p);
        const double psi   = 2.0 * std::numbers::pi * uniform(rng);

        bodies.push_back(Body{
            .mass = masses[i] * kg,
            .position = position_t{
                cartesian_vector<double>{r * sin_t * std::cos(phi),
                                         r * sin_t * std::sin(phi),
                                         r * cos_t} * isq::displacement[m],
                solar_system_center_of_mass},
            .velocity = cartesian_vector<double>{v * sin_p * std::cos(psi),
                                                  v * sin_p * std::sin(psi),
                                                  v * cos_p} * isq::velocity[m / s],
        });
    }

    // Move to center-of-mass frame
    auto total_mass = 0.0 * kg;
    for (const auto& body : bodies) total_mass += body.mass;

    auto v_cm = cartesian_vector<double>{0, 0, 0} * isq::velocity[m / s];
    for (const auto& body : bodies)
        v_cm += (body.mass / total_mass) * body.velocity;
    for (auto& body : bodies)
        body.velocity -= v_cm;

    auto r_cm = cartesian_vector<double>{0, 0, 0} * isq::displacement[m];
    for (const auto& body : bodies)
        r_cm += (body.mass / total_mass) * (body.position - position_t{});
    for (auto& body : bodies)
        body.position -= r_cm;

    return bodies;
}

} // namespace planets
