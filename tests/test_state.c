// File: tests/test_state.c
// VERSION: Corrected function calls and includes.

#include <stdio.h>
#include <stdlib.h>
#include <math.h> // BUG FIX: Added for fabs()
#include "datatypes.h"
#include "hamiltonian.h"
#include "state.h"

#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fprintf(stderr, "\n--- TEST FAILED: %s ---\n", message); \
        exit(1); \
    } else { \
        printf("PASSED: %s\n", message); \
    }

const char* test_filename = "tests/test_pmqmc.in";

int main() {
    printf("\n--- Running Tests for State Module ---\n");

    // --- Setup: Load a known Hamiltonian ---
    SimParams params = {0};
    divdiff_global_init();
    Hamiltonian* h = hamiltonian_create_and_load(test_filename, &params);
    TEST_ASSERT(h != NULL, "Prerequisite: Hamiltonian loaded successfully");

    // --- Test 1: State Initialization ---
    QMCState* state = state_create(&params);
    TEST_ASSERT(state != NULL, "state_create runs without error");
    TEST_ASSERT(state->lattice != NULL, "Lattice was allocated");
    TEST_ASSERT(state->lattice->num_bits == params.N, "Lattice has correct size (N)");
    TEST_ASSERT(state->q == 0, "Initial operator sequence length q is 0");
    TEST_ASSERT(state->Energies != NULL, "Energies array was allocated");
    TEST_ASSERT(state->weight_calculator != NULL, "Weight calculator was allocated");
    
    // Perform initial calculation
    state_recalculate_props(state, h); // BUG FIX: Call the correctly named function
    double initial_energy = state_calculate_classical_energy(h, state->lattice);
    TEST_ASSERT(fabs(state->Energies[0] - initial_energy) < 1e-9, "Initial energy was calculated correctly");
    
    printf("State initialization tests seem OK.\n\n");

    // --- Test 2: Property Recalculation ---
    printf("--- Test 2: Property Recalculation ---\n");
    // Manually create a sequence to test the update function
    state->q = 2;
    state->Sq[0] = 2; // P_matrix["00000011"]
    state->Sq[1] = 3; // P_matrix["00000100"]

    state_recalculate_props(state, h); // BUG FIX: Call the correctly named function
    
    // Create the expected intermediate lattice state
    bitset_t* intermediate_lattice = bitset_create(params.N);
    bitset_copy(intermediate_lattice, state->lattice);
    bitset_xor(intermediate_lattice, h->P_matrix[2]);
    double expected_E1 = state_calculate_classical_energy(h, intermediate_lattice);

    TEST_ASSERT(fabs(state->Energies[1] - expected_E1) < 1e-9, "Energy[1] is correct after one operator");

    bitset_free(intermediate_lattice);
    printf("Property recalculation tests seem OK.\n\n");

    // --- Cleanup ---
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();

    printf("--- All State tests passed. ---\n");
    return 0;
}