#include "cpems/cpeds-direction.h"
#include "gtest/gtest.h"

namespace cpems {
namespace {

TEST(astroCoordTest, canDoPrecession) {
    double alpha = 100; // deg
    double delta = 20;  // deg
    double jd1 = 2451545;
    double jd2 = 2460002.316149432677776;
    DirectionRaDec d = DirectionRaDec(alpha * PI180, delta * PI180, jd1);
    d.print_direction("ref dir jd1");
    d.toJD(jd2);
    d.print_direction("ref dir jd2");
    d.toDeg();
    double alpha_ref = d.ra();
    double delta_ref = d.dec();

    double test_alpha = RT4_PrecessionRA(jd1, jd2, alpha, delta);
    double test_delta = RT4_PrecessionDEC(jd1, jd2, alpha, delta);

    DirectionRaDec test_dir =
        DirectionRaDec(test_alpha * PI180, test_delta * PI180);
    test_dir.print_direction("test dir");

    EXPECT_LE(fabs(test_alpha - alpha_ref), 0.0001);
    EXPECT_LE(fabs(test_delta - delta_ref), 0.0001);
}

TEST(astroCoordTest, canDoPrecessionAndNutationAndAberration) {
    double alpha = 100; // deg
    double delta = 20;  // deg
    double jd1 = 2451545;
    double jd2 = 2460002.316149432677776;
    DirectionRaDec d = DirectionRaDec(alpha * PI180, delta * PI180, jd1);
    d.print_direction("ref dir jd1");
    d.toJD(jd2);
    d.nutate(0);
    d.aberrate(jd2);
    d.print_direction("ref dir jd2");
    d.toDeg();
    double alpha_ref = d.ra();
    double delta_ref = d.dec();

    double test_alpha = RT4_PrecessionRA(jd1, jd2, alpha, delta);
    double test_delta = RT4_PrecessionDEC(jd1, jd2, alpha, delta);

    // cpeds_rt4_nutation(jd2, &test_alpha, &test_delta);
    cpeds_rt4_nutation_and_aberration(jd2, &test_alpha, &test_delta);

    DirectionRaDec test_dir =
        DirectionRaDec(test_alpha * PI180, test_delta * PI180, jd2);
    test_dir.print_direction("test dir");

    EXPECT_LE(fabs(test_alpha - alpha_ref), 0.0001);
    EXPECT_LE(fabs(test_delta - delta_ref), 0.0001);
}

} // namespace
} // namespace cpems
