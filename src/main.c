// File: src/main.c
// Purpose: Main entry point for the PMR-QMC simulation.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "datatypes.h"
#include "hamiltonian.h"
#include "state.h"
#include "updates.h"
#include "utils.h"
#include "measurements.h" // <-- ADD
#include "statistics.h"   // <-- ADD

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    divdiff_global_init();
    rng_init(time(NULL));

    SimParams params = {0};
    Hamiltonian* h = hamiltonian_create_and_load(argv[1], &params);
    if (!h) return 1;

    QMCState* state = state_create(&params);
    if (!state) {
        hamiltonian_free(h, &params);
        return 1;
    }

    Stats* stats = stats_create(&params); // <-- ADD: Create stats object
    if (!stats) {
        state_free(state);
        hamiltonian_free(h, &params);
        return 1;
    }

    // Initial calculation
    state_recalculate_props(state, h);
    rebuild_divdiff_from_energies(state, &params);

    printf("Starting thermalization...\n");
    for (long long i = 0; i < params.TSTEPS; ++i) {
        if (params.WORM) {
            do_worm_update(state, h, &params);
        } else {
            do_composite_update(state, h, &params);
        }
    }

    printf("Starting measurements...\n");
    long long total_measurements = params.STEPS / params.STEPS_PER_MEASUREMENT;
    for (long long i = 0; i < total_measurements; ++i) {
        for(int j = 0; j < params.STEPS_PER_MEASUREMENT; ++j) {
            if (params.WORM) {
                do_worm_update(state, h, &params);
            } else {
                do_composite_update(state, h, &params);
            }
        }
        
        // --- NEW MEASUREMENT AND ACCUMULATION LOGIC ---
        ExExFloat weight = get_full_weight(state);
        double sgn = (exex_to_double(weight) > 0) ? 1.0 : -1.0;
        if (exex_to_double(weight) == 0.0) sgn = 0.0;

        double H_val = measure_H(state, &params);
        double H2_val = measure_H2(state, &params);

        stats_accumulate(stats, sgn, H_val, H2_val);
        // --- END NEW LOGIC ---
    }

    // --- Final analysis and print results ---
    stats_finalize_and_print(stats, &params);

    // Cleanup
    stats_free(stats); // <-- ADD
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();

    return 0;
}