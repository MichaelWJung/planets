#include "plummer.h"

#include <cmath>
#include <gtest/gtest.h>
#include <mp-units/math.h>
#include <ranges>

using namespace mp_units;
using namespace mp_units::si::unit_symbols;

namespace planets {
namespace {

// Shared parameters for all Plummer tests: a smallish cluster, 1000 bodies.
constexpr int    N       = 1000;
constexpr auto   M_total = 1.989e35 * isq::mass[kg];   // ~1e5 solar masses
constexpr auto   a       = 3.086e17 * isq::length[m];  // ~10 pc
constexpr std::uint64_t seed = 42;

TEST(PlummerTest, BodyCount) {
    const auto bodies = plummerBodies(N, M_total, a, seed);
    EXPECT_EQ(static_cast<int>(bodies.size()), N);
}

TEST(PlummerTest, TotalMass) {
    const auto bodies = plummerBodies(N, M_total, a, seed);

    auto sum = 0.0 * kg;
    for (const auto& b : bodies) sum += b.mass;

    EXPECT_NEAR(sum.numerical_value_in(kg),
                M_total.numerical_value_in(kg),
                1e-10 * M_total.numerical_value_in(kg));
}

TEST(PlummerTest, TotalMomentumIsZero) {
    const auto bodies = plummerBodies(N, M_total, a, seed);

    double px = 0, py = 0, pz = 0, p_mag_sum = 0;
    for (const auto& b : bodies) {
        const double mass = b.mass.numerical_value_in(kg);
        const auto   vel  = b.velocity.numerical_value_in(m / s);
        px += mass * vel[0];
        py += mass * vel[1];
        pz += mass * vel[2];
        p_mag_sum += mass * std::sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
    }

    // Each component should be negligible relative to the average individual momentum.
    const double tol = 1e-10 * p_mag_sum / N;
    EXPECT_NEAR(px, 0.0, tol);
    EXPECT_NEAR(py, 0.0, tol);
    EXPECT_NEAR(pz, 0.0, tol);
}

TEST(PlummerTest, CenterOfMassAtOrigin) {
    const auto bodies = plummerBodies(N, M_total, a, seed);

    double rx = 0, ry = 0, rz = 0;
    double total_mass = 0;
    for (const auto& b : bodies) {
        const double mass = b.mass.numerical_value_in(kg);
        const auto   pos  = (b.position - position_t{}).numerical_value_in(m);
        rx += mass * pos[0];
        ry += mass * pos[1];
        rz += mass * pos[2];
        total_mass += mass;
    }
    rx /= total_mass;
    ry /= total_mass;
    rz /= total_mass;

    // CoM should be at the origin to within a tiny fraction of the scale radius.
    const double tol = 1e-10 * a.numerical_value_in(m);
    EXPECT_NEAR(rx, 0.0, tol);
    EXPECT_NEAR(ry, 0.0, tol);
    EXPECT_NEAR(rz, 0.0, tol);
}

TEST(PlummerTest, VirialEquilibrium) {
    // For an exact Plummer distribution, 2T + W = 0 (virial theorem).
    // With N=1000 and a Kroupa mass spectrum there are finite-N fluctuations,
    // so we accept 2T/|W| within 30% of 1.
    const auto bodies = plummerBodies(N, M_total, a, seed);

    double T = 0.0;
    for (const auto& b : bodies) {
        const double mass = b.mass.numerical_value_in(kg);
        const auto   vel  = b.velocity.numerical_value_in(m / s);
        T += 0.5 * mass * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
    }

    const double G_val = G.numerical_value_in(m3 / (kg * s2));
    double W = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            const double mi = bodies[i].mass.numerical_value_in(kg);
            const double mj = bodies[j].mass.numerical_value_in(kg);
            const auto   rij = norm(bodies[j].position - bodies[i].position);
            W -= G_val * mi * mj / rij.numerical_value_in(m);
        }
    }

    EXPECT_NEAR(2.0 * T / std::abs(W), 1.0, 0.3);
}

} // namespace
} // namespace planets
