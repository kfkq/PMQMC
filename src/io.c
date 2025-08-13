// File: src/io.c
// Purpose: Implementation of HDF5 I/O functions for the PM-QMC simulator

#include "io.h"
#include "datatypes.h"
#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>

// Create a new HDF5 file
hid_t io_hdf5_create_file(const char* filename) {
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error: Could not create HDF5 file '%s'.", filename);
        return -1;
    }
    return file_id;
}

// Create the measurement dataset in the HDF5 file
hid_t io_hdf5_create_measurement_dataset(hid_t file_id, const SimParams* params) {
    // Create the compound datatype for measurements
    hid_t memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(Measurement));
    H5Tinsert(memtype_id, "sgn", HOFFSET(Measurement, sgn), H5T_NATIVE_DOUBLE);
    if (params->MEASURE_H) H5Tinsert(memtype_id, "H", HOFFSET(Measurement, H), H5T_NATIVE_DOUBLE);
    if (params->MEASURE_H2) H5Tinsert(memtype_id, "H2", HOFFSET(Measurement, H2), H5T_NATIVE_DOUBLE);
    if (params->MEASURE_Z_MAGNETIZATION) H5Tinsert(memtype_id, "Z_mag", HOFFSET(Measurement, Z_mag), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "q", HOFFSET(Measurement, q), H5T_NATIVE_INT);

    // Create dataspace with unlimited dimensions
    hsize_t dims[1] = {0};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hid_t filespace_id = H5Screate_simple(1, dims, maxdims);
    
    // Create dataset creation property list for chunking and compression
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[1] = {HDF5_BUFFER_SIZE};
    H5Pset_chunk(dcpl_id, 1, chunk_dims);
    H5Pset_deflate(dcpl_id, 6);
    
    // Create the dataset
    hid_t dset_id = H5Dcreate2(file_id, "/measurements", memtype_id, filespace_id, 
                               H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    
    // Write simulation parameters as attributes
    io_hdf5_write_int_attribute(dset_id, "N", params->N);
    io_hdf5_write_double_attribute(dset_id, "BETA", params->BETA);
    io_hdf5_write_long_long_attribute(dset_id, "STEPS", params->STEPS);
    io_hdf5_write_int_attribute(dset_id, "STEPS_PER_MEASUREMENT", params->STEPS_PER_MEASUREMENT);
    io_hdf5_write_long_long_attribute(dset_id, "SKIP_MEASUREMENTS", params->SKIP_MEASUREMENTS);
    io_hdf5_write_int_attribute(dset_id, "NBINS", params->NBINS);
    io_hdf5_write_int_attribute(dset_id, "WORM", params->WORM);
    io_hdf5_write_int_attribute(dset_id, "QMAX", params->QMAX);
    
    // Close resources
    H5Pclose(dcpl_id);
    H5Sclose(filespace_id);
    H5Tclose(memtype_id);
    
    return dset_id;
}

// Setup HDF5 file and dataset
hid_t io_hdf5_setup(const char* filename, const SimParams* params, hid_t* file_id, hid_t* memtype_id) {
    *file_id = io_hdf5_create_file(filename);
    if (*file_id < 0) {
        return -1;
    }
    
    hid_t dset_id = io_hdf5_create_measurement_dataset(*file_id, params);
    if (dset_id < 0) {
        H5Fclose(*file_id);
        return -1;
    }
    
    // Create memory type for writing data
    *memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(Measurement));
    H5Tinsert(*memtype_id, "sgn", HOFFSET(Measurement, sgn), H5T_NATIVE_DOUBLE);
    if (params->MEASURE_H) H5Tinsert(*memtype_id, "H", HOFFSET(Measurement, H), H5T_NATIVE_DOUBLE);
    if (params->MEASURE_H2) H5Tinsert(*memtype_id, "H2", HOFFSET(Measurement, H2), H5T_NATIVE_DOUBLE);
    if (params->MEASURE_Z_MAGNETIZATION) H5Tinsert(*memtype_id, "Z_mag", HOFFSET(Measurement, Z_mag), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype_id, "q", HOFFSET(Measurement, q), H5T_NATIVE_INT);
    
    return dset_id;
}

// Cleanup HDF5 resources
void io_hdf5_cleanup(hid_t file_id, hid_t dset_id, hid_t memtype_id) {
    if (memtype_id >= 0) H5Tclose(memtype_id);
    if (dset_id >= 0) H5Dclose(dset_id);
    if (file_id >= 0) H5Fclose(file_id);
}

// Flush the buffer to the HDF5 file
void io_hdf5_flush_buffer(hid_t dset_id, hid_t memtype_id, Measurement* buffer, int* buffer_count, hsize_t* total_written) {
    if (*buffer_count == 0) return;
    
    hsize_t buffer_dims[1] = { (hsize_t)(*buffer_count) };
    hid_t memspace_id = H5Screate_simple(1, buffer_dims, NULL);
    
    hsize_t new_total_dims[1] = { *total_written + *buffer_count };
    H5Dset_extent(dset_id, new_total_dims);
    
    hid_t filespace_id = H5Dget_space(dset_id);
    hsize_t offset[1] = { *total_written };
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, buffer_dims, NULL);
    
    H5Dwrite(dset_id, memtype_id, memspace_id, filespace_id, H5P_DEFAULT, buffer);
    
    H5Sclose(filespace_id);
    H5Sclose(memspace_id);
    
    *total_written += *buffer_count;
    *buffer_count = 0;
}

// Add a measurement to the buffer
void io_hdf5_write_measurement(Measurement* buffer, int* buffer_count,
                              double sgn, double H, double H2, double Z_mag, int q,
                              int measure_H, int measure_H2, int measure_Z_mag) {
    Measurement* current = &buffer[*buffer_count];
    current->sgn = sgn;
    if (measure_H) current->H = H;
    if (measure_H2) current->H2 = H2;
    if (measure_Z_mag) current->Z_mag = Z_mag;
    current->q = q;
    (*buffer_count)++;
}

// Write an integer attribute to an HDF5 object
void io_hdf5_write_int_attribute(hid_t obj_id, const char* name, int value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(obj_id, name, H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    H5Sclose(attr_space);
}

// Write a long long attribute to an HDF5 object
void io_hdf5_write_long_long_attribute(hid_t obj_id, const char* name, long long value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(obj_id, name, H5T_NATIVE_LLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_LLONG, &value);
    H5Aclose(attr_id);
    H5Sclose(attr_space);
}

// Write a double attribute to an HDF5 object
void io_hdf5_write_double_attribute(hid_t obj_id, const char* name, double value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(obj_id, name, H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr_id);
    H5Sclose(attr_space);
}