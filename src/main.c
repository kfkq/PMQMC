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

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    divdiff_global_init();
    rng_init(time(NULL));  // Or fixed seed for repro

    SimParams params = {0};
    Hamiltonian* h = hamiltonian_create_and_load(argv[1], &params);
    if (!h) return 1;

    QMCState* state = state_create(&params);
    if (!state) {
        hamiltonian_free(h, &params);
        return 1;
    }

    state_recalculate_props(state, h);
    rebuild_divdiff_from_energies(state, &params);

    // Thermalization
    for (long long i = 0; i < params.TSTEPS; ++i) {
        if (params.WORM) {
            do_worm_update(state, h, &params);
        } else {
            do_composite_update(state, h, &params);
        }
    }

    // Measurements (stub: add bins/observables)
    double energy_sum = 0.0;
    for (long long i = 0; i < params.STEPS; ++i) {
        if (params.WORM) {
            do_worm_update(state, h, &params);
        } else {
            do_composite_update(state, h, &params);
        }
        if (i % params.STEPS_PER_MEASUREMENT == 0) {
            // Measure: e.g., <H> = -d(ln Z)/d beta, but stub with classical E
            energy_sum += state->Energies[0];  // Replace with full <H>
        }
    }
    printf("Average classical energy: %.4f\n", energy_sum / (params.STEPS / params.STEPS_PER_MEASUREMENT));

    // Cleanup
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();

    return 0;
}