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

// --- Define filenames ---
const char* INPUT_FILENAME = "hamiltonian.in";
const char* HDF5_DATA_FILENAME = "raw_data.h5";

// --- HDF5 Data Structure and Buffer ---
typedef struct {
    double sgn;
    double H;
    double H2;
    double Z_mag;
    int q;
} Measurement;

#define HDF5_BUFFER_SIZE 1024
static Measurement hdf5_buffer[HDF5_BUFFER_SIZE];
static int buffer_count = 0;
static hsize_t total_measurements_written = 0;

// --- HDF5 Helper Functions ---
static void flush_buffer_to_hdf5(hid_t dset_id, hid_t memtype_id) {
    if (buffer_count == 0) return;
    hsize_t buffer_dims[1] = { (hsize_t)buffer_count };
    hid_t memspace_id = H5Screate_simple(1, buffer_dims, NULL);
    hsize_t new_total_dims[1] = { total_measurements_written + buffer_count };
    H5Dset_extent(dset_id, new_total_dims);
    hid_t filespace_id = H5Dget_space(dset_id);
    hsize_t offset[1] = { total_measurements_written };
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, buffer_dims, NULL);
    H5Dwrite(dset_id, memtype_id, memspace_id, filespace_id, H5P_DEFAULT, hdf5_buffer);
    H5Sclose(filespace_id);
    H5Sclose(memspace_id);
    total_measurements_written += buffer_count;
    buffer_count = 0;
}

static void write_int_attribute(hid_t obj_id, const char* name, int value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(obj_id, name, H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    H5Sclose(attr_space);
}

static void write_long_long_attribute(hid_t obj_id, const char* name, long long value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(obj_id, name, H5T_NATIVE_LLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_LLONG, &value);
    H5Aclose(attr_id);
    H5Sclose(attr_space);
}

static void write_double_attribute(hid_t obj_id, const char* name, double value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(obj_id, name, H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr_id);
    H5Sclose(attr_space);
}

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
    hid_t file_id = H5Fcreate(HDF5_DATA_FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(Measurement));
    H5Tinsert(memtype_id, "sgn", HOFFSET(Measurement, sgn), H5T_NATIVE_DOUBLE);
    if (params.MEASURE_H) H5Tinsert(memtype_id, "H", HOFFSET(Measurement, H), H5T_NATIVE_DOUBLE);
    if (params.MEASURE_H2) H5Tinsert(memtype_id, "H2", HOFFSET(Measurement, H2), H5T_NATIVE_DOUBLE);
    if (params.MEASURE_Z_MAGNETIZATION) H5Tinsert(memtype_id, "Z_mag", HOFFSET(Measurement, Z_mag), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "q", HOFFSET(Measurement, q), H5T_NATIVE_INT);

    hsize_t dims[1] = {0};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hid_t filespace_id = H5Screate_simple(1, dims, maxdims);
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[1] = {HDF5_BUFFER_SIZE};
    H5Pset_chunk(dcpl_id, 1, chunk_dims);
    H5Pset_deflate(dcpl_id, 6);
    hid_t dset_id = H5Dcreate2(file_id, "/measurements", memtype_id, filespace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    // --- THIS IS THE CRITICAL BLOCK THAT WAS MISSING FROM YOUR EXECUTABLE ---
    write_int_attribute(dset_id, "N", params.N);
    write_double_attribute(dset_id, "BETA", params.BETA);
    write_long_long_attribute(dset_id, "STEPS", params.STEPS);
    write_int_attribute(dset_id, "STEPS_PER_MEASUREMENT", params.STEPS_PER_MEASUREMENT);
    write_long_long_attribute(dset_id, "SKIP_MEASUREMENTS", params.SKIP_MEASUREMENTS);
    write_int_attribute(dset_id, "NBINS", params.NBINS);
    write_int_attribute(dset_id, "WORM", params.WORM);
    write_int_attribute(dset_id, "QMAX", params.QMAX);
    // --------------------------------------------------------------------

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
            Measurement* current = &hdf5_buffer[buffer_count];
            ExExFloat weight = get_full_weight(state);
            current->sgn = (exex_to_double(weight) > 0) ? 1.0 : -1.0;
            if (exex_to_double(weight) == 0.0) current->sgn = 0.0;
            if (params.MEASURE_H) current->H = measure_H(state, &params);
            if (params.MEASURE_H2) current->H2 = measure_H2(state, &params);
            if (params.MEASURE_Z_MAGNETIZATION) current->Z_mag = measure_Z_magnetization(state, &params);
            current->q = state->q;
            buffer_count++;
            if (buffer_count == HDF5_BUFFER_SIZE) {
                flush_buffer_to_hdf5(dset_id, memtype_id);
            }
        }
        if (i > 0 && i % print_interval == 0) {
            printf("\rProgress: [%.2f%%]", (double)i / params.STEPS * 100.0);
            fflush(stdout);
        }
    }

    flush_buffer_to_hdf5(dset_id, memtype_id);
    printf("\rProgress: [100.00%%]\n");
    printf("Simulation finished. Raw data written to %s\n", HDF5_DATA_FILENAME);

    // --- 4. Cleanup ---
    H5Pclose(dcpl_id);
    H5Dclose(dset_id);
    H5Sclose(filespace_id);
    H5Tclose(memtype_id);
    H5Fclose(file_id);
    state_free(state);
    hamiltonian_free(h, &params);
    divdiff_global_cleanup();
    return 0;
}