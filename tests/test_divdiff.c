// File: tests/test_divdiff.c
// Purpose: A standalone program to test the divdiff module.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "divdiff.h"

#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fprintf(stderr, "TEST FAILED: %s\n", message); \
        exit(1); \
    } else { \
        printf("TEST PASSED: %s\n", message); \
    }

// Helper to check for near-equality of doubles
void check_close(double a, double b, const char* message) {
    TEST_ASSERT(fabs(a - b) < 1e-9, message);
}

int main() {
    printf("\n--- Running Tests for OPTIMIZED DivDiff Module ---\n");
    divdiff_global_init();

    // --- Test 1: Incremental Add/Remove Cycle ---
    printf("--- Test 1: Incremental Add/Remove Cycle ---\n");
    DivDiff* dd = divdiff_create(10, 50); // max_q=10, s_max=50
    
    double E0 = -1.5, E1 = -0.5, E2 = -1.0;
    
    // State 0: Empty
    TEST_ASSERT(dd->current_len == 0, "Initial length is 0");

    // State 1: Add E0
    divdiff_add_element(dd, E0);
    double expected_d0 = exp(E0);
    check_close(exex_to_double(dd->results[0]), expected_d0, "d[E0] is correct");

    // State 2: Add E1
    divdiff_add_element(dd, E1);
    double expected_d1 = (exp(E0) - exp(E1)) / (E0 - E1);
    check_close(exex_to_double(dd->results[1]), expected_d1, "d[E0, E1] is correct after add");

    // State 3: Add E2
    divdiff_add_element(dd, E2);
    double d01 = (exp(E0)-exp(E1))/(E0-E1);
    double d12 = (exp(E1)-exp(E2))/(E1-E2);
    double expected_d2 = (d01 - d12) / (E0 - E2);
    check_close(exex_to_double(dd->results[2]), expected_d2, "d[E0, E1, E2] is correct after add");

    // State 4: Remove E2
    divdiff_remove_last_element(dd);
    TEST_ASSERT(dd->current_len == 2, "Length is 2 after remove");
    // The previous result should still be accessible and correct
    check_close(exex_to_double(dd->results[1]), expected_d1, "d[E0, E1] is correct after remove");

    divdiff_free(dd);
    printf("Incremental tests seem OK.\n\n");

    divdiff_global_cleanup();
    printf("--- All Optimized DivDiff tests passed. ---\n");
    return 0;
}