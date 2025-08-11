// File: tests/test_state.c

#include <stdio.h>
#include <stdlib.h>
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
    int load_status = hamiltonian_load(test_filename, &params);
    TEST_ASSERT(load_status == 0, "Prerequisite: Hamiltonian loaded successfully");

    // --- Test 1: State Initialization ---
    int init_status = state_init(&params);
    TEST_ASSERT(init_status == 0, "state_init runs without error");
    TEST_ASSERT(lattice != NULL, "Lattice was allocated");
    TEST_ASSERT(lattice->num_bits == params.N, "Lattice has correct size (N)");
    TEST_ASSERT(q == 0, "Initial operator sequence length q is 0");
    TEST_ASSERT(Energies != NULL, "Energies array was allocated");
    TEST_ASSERT(weight_calculator != NULL, "Weight calculator was allocated");
    
    // Test initial energy calculation
    double initial_energy = state_calculate_classical_energy(lattice);
    TEST_ASSERT(Energies[0] == initial_energy, "Initial energy was calculated correctly");
    
    // Test initial weight calculation
    TEST_ASSERT(weight_calculator->current_len == 1, "Weight calculator has one initial energy");
    TEST_ASSERT(weight_calculator->energies[0] == initial_energy, "Weight calculator stores correct initial energy");

    printf("State initialization tests seem OK.\n\n");

    // --- Test 2: Energy List Update ---
    printf("--- Test 2: Energy List Update ---\n");
    // Manually create a sequence to test the update function
    q = 2;
    Sq[0] = 2; // P_matrix["00000011"]
    Sq[1] = 3; // P_matrix["00000100"]

    state_update_energy_list();
    
    // Create the expected intermediate lattice state
    bitset_t* intermediate_lattice = bitset_create(params.N);
    bitset_copy(intermediate_lattice, lattice);
    bitset_xor(intermediate_lattice, P_matrix[2]);
    double expected_E1 = state_calculate_classical_energy(intermediate_lattice);

    TEST_ASSERT(Energies[1] == expected_E1, "Energy[1] is correct after one operator");

    bitset_free(intermediate_lattice);
    printf("Energy list update tests seem OK.\n\n");

    // --- Cleanup ---
    state_free();
    hamiltonian_free(&params);
    divdiff_global_cleanup();

    printf("--- All State tests passed. ---\n");
    return 0;
}