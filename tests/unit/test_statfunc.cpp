#include <gtest/gtest.h>

#include "StatFunc.h"

// Light StatFunc smoke tests. Avoid LOGGER.e death paths.

TEST(StatFunc, PchisqDf1KnownValues) {
    // chi^2(1) critical values: P(X > 3.841) ≈ 0.05, P(X > 6.635) ≈ 0.01
    EXPECT_NEAR(StatFunc::pchisq(3.8414588, 1.0), 0.05, 1e-4);
    EXPECT_NEAR(StatFunc::pchisq(6.6348966, 1.0), 0.01, 1e-4);
    EXPECT_NEAR(StatFunc::pchisq(0.0, 1.0), 1.0, 1e-10);
}

TEST(StatFunc, QchisqRoundTrip) {
    const double p = 0.05;
    const double q = StatFunc::qchisq(p, 1.0);
    EXPECT_NEAR(StatFunc::pchisq(q, 1.0), p, 1e-6);
}
