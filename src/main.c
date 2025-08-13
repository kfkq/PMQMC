// File: src/main.c
// Purpose: Main entry point for the PMR-QMC simulation.
// FINAL CORRECTED VERSION

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "datatypes.h"
#include "hamiltonian.h"
#include "state.h"
#include "updates.h"
#include "utils.h"
#include "measurements.h"
#include "hdf5.h"
#include "io.h"

// --- Define filenames ---
const char* INPUT_FILENAME = "hamiltonian.in";
const char* HDF5_DATA_FILENAME = "raw_data.h5";

// --- HDF5 Data Structure and Buffer ---
static Measurement hdf5_buffer[HDF5_BUFFER_SIZE];
static int buffer_count = 0;
static hsize_t total_measurements_written = 0;

int main(void) {
    // --- 1. Initialization and File Checks ---
    FILE* check_file = fopen(INPUT_FILENAME, "r");
    if (!check_file) {
        fprintf(stderr, "Error: Could not find or open the required input file '%s'.\n", INPUT_FILENAME);
        return 1;
    }
    fclose(check_file);

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

    // --- 2. HDF5 File and Dataset Setup ---
    hid_t file_id, memtype_id;
    hid_t dset_id = io_hdf5_setup(HDF5_DATA_FILENAME, &params, &file_id, &memtype_id);
    if (dset_id < 0) {
        fprintf(stderr, "Error: Could not setup HDF5 file.\n");
        hamiltonian_free(h, &params);
        state_free(state);
        divdiff_global_cleanup();
        return 1;
    }

    // --- 3. Simulation Loop ---
    state_recalculate_props(state, h);
    rebuild_divdiff_from_energies(state, &params);
    printf("Starting simulation for %lld total steps...\n", params.STEPS);
    const long long print_interval = (params.STEPS > 100) ? (params.STEPS / 100) : 1;

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
            double H = 0.0, H2 = 0.0, Z_mag = 0.0;
            if (params.MEASURE_H) H = measure_H(state, &params);
            if (params.MEASURE_H2) H2 = measure_H2(state, &params);
            if (params.MEASURE_Z_MAGNETIZATION) Z_mag = measure_Z_magnetization(state, &params);
            int q = state->q;
            
            io_hdf5_write_measurement(hdf5_buffer, &buffer_count, sgn, H, H2, Z_mag, q,
                                     params.MEASURE_H, params.MEASURE_H2, params.MEASURE_Z_MAGNETIZATION);
            
            if (buffer_count == HDF5_BUFFER_SIZE) {
                io_hdf5_flush_buffer(dset_id, memtype_id, hdf5_buffer, &buffer_count, &total_measurements_written);
            }
        }
        if (i > 0 && i % print_interval == 0) {
            printf("\rProgress: [%.2f%%]", (double)i / params.STEPS * 100.0);
            fflush(stdout);
        }
    }

    io_hdf5_flush_buffer(dset_id, memtype_id, hdf5_buffer, &buffer_count, &total_measurements_written);
    printf("\rProgress: [100.00%%]\n");
    printf("Simulation finished. Raw data written to %s\n", HDF5_DATA_FILENAME);

    // --- 4. Cleanup ---
    io_hdf5_cleanup(file_id, dset_id, memtype_id);
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();
    return 0;
}