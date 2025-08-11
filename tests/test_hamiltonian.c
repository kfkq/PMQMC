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

const char* test_filename = "tests/test_pmqmc.in";

int main() {
    printf("\n--- Running Tests for Hamiltonian Module ---\n");
    
    SimParams params = {0};
    int status = hamiltonian_load(test_filename, &params);
    TEST_ASSERT(status == 0, "hamiltonian_load runs without error");

    // Test parameters from the new file
    TEST_ASSERT(params.N == 8, "Parameter N is correct (8)");
    TEST_ASSERT(params.NOP == 16, "Parameter NOP is correct (16)");
    TEST_ASSERT(params.NCYCLES == 8, "Parameter NCYCLES is correct (8)");
    TEST_ASSERT(fabs(params.BETA - 10.0) < 1e-9, "Parameter BETA is correct (10.0)");

    // Test a few data points from P_MATRIX
    TEST_ASSERT(P_matrix != NULL, "P_matrix was allocated");
    TEST_ASSERT(bitset_get(P_matrix[2], 0) == 1, "P_matrix[2] bit 0 is correct ('...0011')");
    TEST_ASSERT(bitset_get(P_matrix[2], 1) == 1, "P_matrix[2] bit 1 is correct ('...0011')");
    TEST_ASSERT(bitset_get(P_matrix[2], 2) == 0, "P_matrix[2] bit 2 is correct ('...0011')");
    TEST_ASSERT(bitset_count(P_matrix[15]) == 2, "P_matrix[15] count is correct ('11000000')");

    // Test a few data points from DIAGONAL_TERM
    TEST_ASSERT(D0_size == 8, "D0_size is correct (8)");
    TEST_ASSERT(creal(D0_coeff[7]) == -1.0, "D0_coeff[7] is correct");
    TEST_ASSERT(bitset_count(D0_product[7]) == 2, "D0_product[7] count is correct ('10000001')");

    // Test a few data points from OFF_DIAGONAL_TERMS
    TEST_ASSERT(D_sizes != NULL, "D_sizes was allocated");
    TEST_ASSERT(D_sizes[2] == 2, "D_sizes[2] is correct (2)");
    TEST_ASSERT(D_sizes[15] == 2, "D_sizes[15] is correct (2)");
    
    TEST_ASSERT(D_coeffs != NULL, "D_coeffs was allocated");
    TEST_ASSERT(fabs(creal(D_coeffs[2][0]) - -1.0) < 1e-9, "D_coeffs[2][0] real part is correct");
    TEST_ASSERT(fabs(creal(D_coeffs[2][1]) - 1.0) < 1e-9, "D_coeffs[2][1] real part is correct");
    
    TEST_ASSERT(D_products != NULL, "D_products was allocated");
    TEST_ASSERT(bitset_count(D_products[2][1]) == 2, "D_products[2][1] count is correct ('00000011')");

    // Cleanup
    hamiltonian_free(&params);

    printf("--- All Hamiltonian tests passed. ---\n");
    return 0;
}