// File: tests/test_updates.c
// Purpose: A standalone program to test the QMC updates module.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "datatypes.h"
#include "hamiltonian.h"
#include "state.h"
#include "updates.h"
#include "utils.h"

#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fprintf(stderr, "\n--- TEST FAILED: %s ---\n", message); \
        exit(1); \
    } else { \
        printf("PASSED: %s\n", message); \
    }

// Helper to print the current state for debugging
void print_state(const QMCState* state) {
    printf("  q = %d, Sq = [", state->q);
    for (int i = 0; i < state->q; ++i) printf("%d ", state->Sq[i]);
    printf("], lattice = ");
    bitset_print(state->lattice);
    printf("\n");
}

int main() {
    printf("\n--- Running Tests for All Update Types ---\n");

    // --- 1. Setup ---
    rng_init(42);
    divdiff_global_init();

    SimParams params = {0};
    Hamiltonian* h = hamiltonian_create_and_load("tests/test_pmqmc.in", &params);
    TEST_ASSERT(h != NULL, "Prerequisite: Hamiltonian loaded successfully");

    QMCState* state = state_create(&params);
    TEST_ASSERT(state != NULL, "Prerequisite: QMC state created successfully");

    state_recalculate_props(state, h);
    rebuild_divdiff_from_energies(state, &params);
    
    printf("Initial State:\n");
    print_state(state);
    TEST_ASSERT(state->q == 0, "Initial q is 0");

    // --- 2. Run a series of composite updates ---
    printf("\n--- Testing Composite Update Stability ---\n");
    int num_updates = 500;
    int q_sum = 0;
    int lattice_flips = 0;
    bitset_t* initial_lattice = bitset_create(params.N);
    bitset_copy(initial_lattice, state->lattice);

    for (int i = 0; i < num_updates; ++i) {
        // In a real simulation, one would choose between composite and worm updates.
        // Here, we test them separately.
        do_worm_update(state, h, &params);
        q_sum += state->q;
        // Check if the lattice has changed from the start
        int changed = 0;
        for(int b=0; b<params.N; ++b) {
            if(bitset_get(initial_lattice, b) != bitset_get(state->lattice, b)) {
                changed = 1;
                break;
            }
        }
        if(changed) lattice_flips++;
    }
    
    printf("Ran %d worm updates.\n", num_updates);
    printf("Average q = %.2f\n", (double)q_sum / num_updates);
    printf("Lattice state was different from initial state in %d of %d steps.\n", lattice_flips, num_updates);
    
    TEST_ASSERT(q_sum > 0, "q was able to increase from 0");

    printf("Final State after worm updates:\n");
    print_state(state);
    
    // --- 3. Test Composite Update (for comparison) ---
    printf("\n--- Testing Composite Update (Stub) ---\
");
    do_composite_update(state, h, &params);
    printf("Composite update placeholder executed without crashing.\n");
    print_state(state);

    // --- 4. Final Cleanup ---
    bitset_free(initial_lattice);
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();

    printf("\n--- All updates tests completed successfully! ---\n");
    return 0;
}