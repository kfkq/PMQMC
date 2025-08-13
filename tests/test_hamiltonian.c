// File: tests/test_hamiltonian.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datatypes.h"
#include "hamiltonian.h"

#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fprintf(stderr, "\n--- TEST FAILED: %s ---\n", message); \
        exit(1); \
    } else { \
        printf("PASSED: %s\n", message); \
    }

void test_hamiltonian_with_worm(const char* test_filename) {
    printf("\n--- Running Tests for Hamiltonian Module ---\n");
    
    SimParams params = {0};
    Hamiltonian* h = hamiltonian_create_and_load(test_filename, &params);
    TEST_ASSERT(h != NULL, "hamiltonian_create_and_load runs without error");

    // Test parameters from the new file
    TEST_ASSERT(params.N == 8, "Parameter N is correct (8)");
    TEST_ASSERT(params.NOP == 16, "Parameter NOP is correct (16)");
    TEST_ASSERT(params.NCYCLES == 8, "Parameter NCYCLES is correct (8)");
    TEST_ASSERT(fabs(params.BETA - 1.0) < 1e-9, "Parameter BETA is correct (1.0)");
    TEST_ASSERT(params.TSTEPS == 10000, "Parameter TSTEPS is correct (10000)");
    TEST_ASSERT(params.STEPS == 1000000, "Parameter STEPS is correct (1000000)");
    TEST_ASSERT(params.STEPS_PER_MEASUREMENT == 10, "Parameter STEPS_PER_MEASUREMENT is correct (10)");
    TEST_ASSERT(params.QMAX == 1000, "Parameter QMAX is correct (1000)");
    TEST_ASSERT(params.NBINS == 150, "Parameter NBINS is correct (150)");
    TEST_ASSERT(params.WORM == 0, "Parameter WORM is correct");

    // Test a few data points from P_MATRIX
    TEST_ASSERT(h->P_matrix != NULL, "P_matrix was allocated");
    TEST_ASSERT(bitset_get(h->P_matrix[2], 0) == 1, "P_matrix[2] bit 0 is correct ('...0011')");
    TEST_ASSERT(bitset_get(h->P_matrix[2], 1) == 1, "P_matrix[2] bit 1 is correct ('...0011')");
    TEST_ASSERT(bitset_get(h->P_matrix[2], 2) == 0, "P_matrix[2] bit 2 is correct ('...0011')");
    TEST_ASSERT(bitset_count(h->P_matrix[15]) == 2, "P_matrix[15] count is correct ('11000000')");

    // Test a few data points from DIAGONAL_TERM
    TEST_ASSERT(h->D0_size == 8, "D0_size is correct (8)");
    TEST_ASSERT(creal(h->D0_coeff[7]) == -1.0, "D0_coeff[7] is correct");
    TEST_ASSERT(bitset_count(h->D0_product[7]) == 2, "D0_product[7] count is correct ('10000001')");

    // Test a few data points from OFF_DIAGONAL_TERMS
    TEST_ASSERT(h->D_sizes != NULL, "D_sizes was allocated");
    TEST_ASSERT(h->D_sizes[2] == 2, "D_sizes[2] is correct (2)");
    TEST_ASSERT(h->D_sizes[15] == 2, "D_sizes[15] is correct (2)");
    
    TEST_ASSERT(h->D_coeffs != NULL, "D_coeffs was allocated");
    TEST_ASSERT(fabs(creal(h->D_coeffs[2][0]) - -1.0) < 1e-9, "D_coeffs[2][0] real part is correct");
    TEST_ASSERT(fabs(creal(h->D_coeffs[2][1]) - 1.0) < 1e-9, "D_coeffs[2][1] real part is correct");
    
    TEST_ASSERT(h->D_products != NULL, "D_products was allocated");
    TEST_ASSERT(bitset_count(h->D_products[2][1]) == 2, "D_products[2][1] count is correct ('00000011')");

    // Cleanup
    hamiltonian_free(h, &params);
}

int main() {
    // Test with WORM = True
    test_hamiltonian_with_worm("tests/test_pmqmc.in");
    printf("--- All Hamiltonian tests passed ---\n");
    
    return 0;
}