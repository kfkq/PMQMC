// File: src/io.h
// Purpose: HDF5 I/O functions for the PM-QMC simulator

#ifndef IO_H
#define IO_H

#include <hdf5.h>
#include "datatypes.h"

// Measurement structure for HDF5 output
typedef struct {
    double sgn;
    double H;
    double H2;
    double Z_mag;
    int q;
} Measurement;

// HDF5 buffer constants
#define HDF5_BUFFER_SIZE 1024

// HDF5 file handling functions
hid_t io_hdf5_create_file(const char* filename);
hid_t io_hdf5_create_measurement_dataset(hid_t file_id, const SimParams* params);
hid_t io_hdf5_setup(const char* filename, const SimParams* params, hid_t* file_id, hid_t* memtype_id);
void io_hdf5_cleanup(hid_t file_id, hid_t dset_id, hid_t memtype_id);

// Buffer management functions
void io_hdf5_flush_buffer(hid_t dset_id, hid_t memtype_id, Measurement* buffer, int* buffer_count, hsize_t* total_written);
void io_hdf5_write_measurement(Measurement* buffer, int* buffer_count, 
                              double sgn, double H, double H2, double Z_mag, int q,
                              int measure_H, int measure_H2, int measure_Z_mag);

// Attribute writing functions
void io_hdf5_write_int_attribute(hid_t obj_id, const char* name, int value);
void io_hdf5_write_long_long_attribute(hid_t obj_id, const char* name, long long value);
void io_hdf5_write_double_attribute(hid_t obj_id, const char* name, double value);

#endif // IO_H