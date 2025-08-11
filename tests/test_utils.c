// File: tests/test_utils.c

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"

#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fprintf(stderr, "\n--- TEST FAILED: %s ---\n", message); \
        exit(1); \
    } else { \
        printf("PASSED: %s\n", message); \
    }

int main() {
    printf("\n--- Running Tests for Utils Module ---\n");

    // --- Test 1: Deterministic Seeding ---
    // Seeding the RNG with the same seed should produce the same sequence of numbers.
    rng_init(12345);
    int first_val = rng_uniform_int(1000);
    double second_val = rng_uniform_double();

    rng_init(12345); // Re-seed with the same value
    TEST_ASSERT(rng_uniform_int(1000) == first_val, "Seeding is deterministic (integer)");
    TEST_ASSERT(rng_uniform_double() == second_val, "Seeding is deterministic (double)");
    printf("Deterministic seeding tests seem OK.\n\n");

    // --- Test 2: Range Checks ---
    printf("--- Test 2: Range Checks ---\n");
    rng_init(time(NULL)); // Use a different seed for range checks

    // Check double range [0.0, 1.0)
    for (int i = 0; i < 50000; ++i) {
        double val = rng_uniform_double();
        if (val < 0.0 || val >= 1.0) {
            // Fail explicitly if out of range
            TEST_ASSERT(0, "Double value was out of the expected [0.0, 1.0) range");
        }
    }
    TEST_ASSERT(1, "All doubles were in the expected [0.0, 1.0) range");

    // Check integer range [0, max-1]
    int max_val = 100;
    for (int i = 0; i < 50000; ++i) {
        int val = rng_uniform_int(max_val);
        if (val < 0 || val >= max_val) {
            // Fail explicitly if out of range
            TEST_ASSERT(0, "Integer value was out of the expected [0, max) range");
        }
    }
    TEST_ASSERT(1, "All integers were in the expected [0, max) range");
    printf("Range check tests seem OK.\n\n");

    printf("--- All Utils tests passed. ---\n");
    return 0;
}