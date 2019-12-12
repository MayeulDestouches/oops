/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_EXPECT_H_
#define OOPS_UTIL_EXPECT_H_

#include <sstream>

// IMPORTANT: To use the macro below, it is also necessary to include "eckit/testing/Test.h",
// after defining ECKIT_TESTING_SELF_REGISTER_CASES if needed.

// Provides more informative output on failure than the raw EXPECT() macro.
#define EXPECT_EQUAL(expr, expected) \
    do { \
        if (!((expr) == (expected))) { \
            std::stringstream str; \
            str << ("EXPECT condition '" #expr " == " #expected "' failed. ") \
                << "(Received: " << expr << "; expected: " << expected << ")"; \
            throw eckit::testing::TestException(str.str(), Here()); \
        } \
    } while (false)

#endif  // OOPS_UTIL_EXPECT_H_