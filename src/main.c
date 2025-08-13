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
#include "processing.h"
#include <mpi.h>

// --- Define filenames ---
const char* INPUT_FILENAME = "hamiltonian.in";
const char* HDF5_DATA_FILENAME_TEMPLATE = "raw_mpi_data/data%02d.h5";
char HDF5_DATA_FILENAME[256];

// --- MPI Variables ---
static int mpi_rank = 0;
static int mpi_size = 1;

// --- HDF5 Data Structure and Buffer ---
static Measurement hdf5_buffer[HDF5_BUFFER_SIZE];
static int buffer_count = 0;
static hsize_t total_measurements_written = 0;

int main(void) {
    // Initialize MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    // Set up rank-specific filename
    snprintf(HDF5_DATA_FILENAME, sizeof(HDF5_DATA_FILENAME), HDF5_DATA_FILENAME_TEMPLATE, mpi_rank);
    
    // Only rank 0 creates the directory
    if (mpi_rank == 0) {
        if (processing_init_mpi_data_directory() != 0) {
            fprintf(stderr, "Error: Could not initialize MPI data directory.\n");
            MPI_Finalize();
            return 1;
        }
    }
    
    // Synchronize all processes
    MPI_Barrier(MPI_COMM_WORLD);

    // --- 1. Initialization and File Checks ---
    FILE* check_file = fopen(INPUT_FILENAME, "r");
    if (!check_file) {
        fprintf(stderr, "Error: Could not find or open the required input file '%s'.\n", INPUT_FILENAME);
        MPI_Finalize();
        return 1;
    }
    fclose(check_file);

    divdiff_global_init();
    // Use a different seed for each MPI rank to ensure different random sequences
    rng_init(time(NULL) + mpi_rank);

    SimParams params = {0};
    Hamiltonian* h = hamiltonian_create_and_load(INPUT_FILENAME, &params);
    if (!h) {
        MPI_Finalize();
        return 1;
    }

    QMCState* state = state_create(&params);
    if (!state) {
        hamiltonian_free(h, &params);
        MPI_Finalize();
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
        MPI_Finalize();
        return 1;
    }

    // --- 3. Simulation Loop ---
    state_recalculate_props(state, h);
    rebuild_divdiff_from_energies(state, &params);
    // Initialize display for this rank
    printf("\033[%d;1HRank %d Progress: [0.00%%]", mpi_rank + 1, mpi_rank);
    fflush(stdout);
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
            // Move cursor to the line for this rank and update progress
            printf("\033[%d;1HRank %d Progress: [%.2f%%]", mpi_rank + 1, mpi_rank, (double)i / params.STEPS * 100.0);
            fflush(stdout);
        }
    }

    io_hdf5_flush_buffer(dset_id, memtype_id, hdf5_buffer, &buffer_count, &total_measurements_written);
    // Move to a new line after completion and show completion message
    printf("\033[%d;1H\033[KRank %d Progress: [100.00%%] - Completed\n", mpi_rank + 1, mpi_rank);
    printf("Rank %d simulation finished. Raw data written to %s\n", mpi_rank, HDF5_DATA_FILENAME);

    // --- 4. Cleanup ---
    io_hdf5_cleanup(file_id, dset_id, memtype_id);
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();
    
    // Synchronize all processes before post-processing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Only rank 0 performs post-processing
    if (mpi_rank == 0) {
        if (processing_combine_mpi_data(mpi_size, "raw_data.h5") == 0) {
            printf("Successfully combined data from all %d ranks.\n", mpi_size);
            processing_cleanup_mpi_data(mpi_size);
            printf("Cleaned up temporary MPI data files.\n");
        } else {
            fprintf(stderr, "Error: Failed to combine MPI data.\n");
        }
    }
    
    MPI_Finalize();
    return 0;
}