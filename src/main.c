// File: src/main.c
// Purpose: Main entry point for the PMR-QMC simulation. Dumps raw data.
// EDIT: Added a progress indicator for long simulations.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "datatypes.h"
#include "hamiltonian.h"
#include "state.h"
#include "updates.h"
#include "utils.h"
#include "measurements.h"

// --- Define the required input and output filenames ---
const char* INPUT_FILENAME = "hamiltonian.in";
const char* RAW_DATA_FILENAME = "raw.dat";

int main(void) {
    // --- Check for the existence of the input file first ---
    FILE* check_file = fopen(INPUT_FILENAME, "r");
    if (!check_file) {
        fprintf(stderr, "Error: Could not find or open the required input file '%s'.\n", INPUT_FILENAME);
        fprintf(stderr, "Please ensure the file exists in the current directory.\n");
        return 1;
    }
    fclose(check_file); // Close the file, we just needed to check it

    divdiff_global_init();
    rng_init(time(NULL));

    SimParams params = {0};
    Hamiltonian* h = hamiltonian_create_and_load(INPUT_FILENAME, &params);
    if (!h) return 1;

    QMCState* state = state_create(&params);
    if (!state) {
        hamiltonian_free(h, &params);
        return 1;
    }

    // Open the output file for writing raw data
    FILE* raw_out_file = fopen(RAW_DATA_FILENAME, "w");
    if (!raw_out_file) {
        perror("Error opening raw data file for writing");
        state_free(state);
        hamiltonian_free(h, &params);
        return 1;
    }
    // Write a header for clarity
    fprintf(raw_out_file, "# sgn H H2 q\n");

    // Initial calculation
    state_recalculate_props(state, h);
    rebuild_divdiff_from_energies(state, &params);

    printf("Starting simulation for %lld total steps...\n", params.STEPS);

    // --- NEW: Progress indicator setup ---
    const long long updates_to_print = 100; // How many times to print progress
    const long long print_interval = (params.STEPS > updates_to_print) ? (params.STEPS / updates_to_print) : 1;

    for (long long i = 0; i < params.STEPS; ++i) {
        if (params.WORM) {
            do_worm_update(state, h, &params);
        } else {
            do_composite_update(state, h, &params);
        }

        if (i % params.STEPS_PER_MEASUREMENT == 0) {
            ExExFloat weight = get_full_weight(state);
            double sgn = (exex_to_double(weight) > 0) ? 1.0 : -1.0;
            if (exex_to_double(weight) == 0.0) sgn = 0.0;

            double H_val = measure_H(state, &params);
            double H2_val = measure_H2(state, &params);

            // Write raw data directly to the file
            fprintf(raw_out_file, "%.8f %.8f %.8f %d\n", sgn, H_val, H2_val, state->q);
        }

        // --- NEW: Print progress to the console ---
        if (i > 0 && i % print_interval == 0) {
            double percent_done = (double)i / params.STEPS * 100.0;
            printf("\rProgress: [%.2f%%]", percent_done);
            fflush(stdout); // Force the output to be written to the terminal
        }
    }

    // --- NEW: Finalize progress indicator ---
    printf("\rProgress: [100.00%%]\n");

    printf("Simulation finished. Raw data written to %s\n", RAW_DATA_FILENAME);

    // Cleanup
    fclose(raw_out_file);
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();

    return 0;
}