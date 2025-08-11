// File: tests/test_divdiff.c
// Purpose: A standalone program to test the divdiff module.
// VERSION: Corrected to test for factorial-scaled results and expanded with more points.

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
    double diff = fabs(a - b);
    double tolerance = 1e-9;
    if (diff > tolerance) {
        fprintf(stderr, "DEBUG: Expected %.8e, got %.8e, diff = %.8e\n", b, a, diff);
    }
    TEST_ASSERT(diff < tolerance, message);
}

// Recursive helper to calculate the expected UN-SCALED divided difference
// This provides the ground truth for our tests.
static double calculate_expected_naive(const double* energies, int n) {
    if (n == 0) {
        return exp(energies[0]);
    }
    // Note: energies+1 is a pointer to the start of the sub-array {E1, E2, ...}
    double term1 = calculate_expected_naive(energies, n - 1);
    double term2 = calculate_expected_naive(energies + 1, n - 1);
    return (term1 - term2) / (energies[0] - energies[n]);
}


int main() {
    printf("\n--- Running Tests for OPTIMIZED DivDiff Module ---\n");
    divdiff_global_init();

    // --- Precompute factorials for the test ---
    double factorials[11];
    factorials[0] = 1.0;
    for(int i = 1; i <= 10; ++i) {
        factorials[i] = factorials[i-1] * i;
    }

    // --- Test 1: Longer Incremental Add/Remove Cycle ---
    printf("\n--- Test 1: Longer Incremental Add/Remove Cycle ---\n");
    DivDiff* dd = divdiff_create(10, 50);
    
    double energies[] = {-1.5, -0.5, -1.0, 0.5, -2.0};
    
    // --- Add elements one by one and check at each step ---
    printf("--- Adding elements ---\n");
    for (int i = 0; i < 5; ++i) {
        char msg[100];
        sprintf(msg, "d[E0...E%d] is correct after add", i);
        
        divdiff_add_element(dd, energies[i]);
        
        double expected_val = calculate_expected_naive(energies, i);
        double actual_val = exex_to_double(dd->results[i]) / factorials[i];
        
        check_close(actual_val, expected_val, msg);
    }

    // --- Remove elements one by one and check at each step ---
    printf("\n--- Removing elements ---\n");
    for (int i = 4; i >= 1; --i) {
        char msg[100];
        sprintf(msg, "d[E0...E%d] is correct after remove", i - 1);

        divdiff_remove_last_element(dd);
        TEST_ASSERT(dd->current_len == i, "Length is correct after remove");

        double expected_val = calculate_expected_naive(energies, i - 1);
        double actual_val = exex_to_double(dd->results[i - 1]) / factorials[i - 1];

        check_close(actual_val, expected_val, msg);
    }

    divdiff_free(dd);
    printf("Longer incremental tests seem OK.\n\n");

    // --- Test 2: Repeated Energies ---
    printf("--- Test 2: Repeated Energies ---\n");
    DivDiff* dd_rep = divdiff_create(10, 50);
    double E0 = -1.5, E1 = -0.5;

    divdiff_add_element(dd_rep, E0); // d[E0]
    divdiff_add_element(dd_rep, E1); // d[E0, E1]
    divdiff_add_element(dd_rep, E1); // d[E0, E1, E1]

    // The formula for repeated points involves derivatives:
    // d[z0, z1, z1] = (d[z0, z1] - d[z1, z1]) / (z0 - z1)
    // where d[z1, z1] = f'(z1). For f(x)=exp(x), f'(x)=exp(x).
    double d_01 = (exp(E0) - exp(E1)) / (E0 - E1);
    double d_11 = exp(E1); // The derivative
    double expected_d2_rep = (d_01 - d_11) / (E0 - E1);
    
    double actual_d2_rep = exex_to_double(dd_rep->results[2]) / factorials[2];
    check_close(actual_d2_rep, expected_d2_rep, "d[E0, E1, E1] with repeated energy is correct");

    divdiff_free(dd_rep);
    printf("Repeated energy tests seem OK.\n\n");


    divdiff_global_cleanup();
    printf("----------------------------------------\n");
    printf("All Optimized DivDiff tests passed!\n");
    printf("----------------------------------------\n");
    return 0;
}