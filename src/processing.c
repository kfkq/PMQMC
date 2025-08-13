// File: src/processing.c
// Purpose: Implementation of post-processing functions for MPI-QMC simulation data

#include "processing.h"
#include "datatypes.h"
#include "io.h"
#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

// Create directory for MPI raw data files
int processing_init_mpi_data_directory(void) {
    const char* dirname = "raw_mpi_data";
    
    // Try to create directory
    if (mkdir(dirname, 0755) == 0) {
        return 0; // Success
    }
    
    // If directory already exists, that's fine
    if (errno == EEXIST) {
        return 0;
    }
    
    // Other error
    fprintf(stderr, "Error: Could not create directory '%s'\n", dirname);
    return -1;
}

// Generate rank-specific filename
void processing_get_rank_filename(char* filename, int rank, size_t filename_size) {
    snprintf(filename, filename_size, "raw_mpi_data/data%02d.h5", rank);
}

// Read data from a single rank's HDF5 file
static int read_rank_data(const char* filename, Measurement** measurements, hsize_t* num_measurements) {
    // Open file
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error: Could not open HDF5 file '%s'\n", filename);
        return -1;
    }
    
    // Open dataset
    hid_t dset_id = H5Dopen2(file_id, "/measurements", H5P_DEFAULT);
    if (dset_id < 0) {
        fprintf(stderr, "Error: Could not open dataset in file '%s'\n", filename);
        H5Fclose(file_id);
        return -1;
    }
    
    // Get dataspace and dimensions
    hid_t dataspace_id = H5Dget_space(dset_id);
    if (dataspace_id < 0) {
        fprintf(stderr, "Error: Could not get dataspace from file '%s'\n", filename);
        H5Dclose(dset_id);
        H5Fclose(file_id);
        return -1;
    }
    
    hsize_t dims[1];
    int ndims = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    if (ndims != 1) {
        fprintf(stderr, "Error: Unexpected dataspace dimensions in file '%s'\n", filename);
        H5Sclose(dataspace_id);
        H5Dclose(dset_id);
        H5Fclose(file_id);
        return -1;
    }
    
    *num_measurements = dims[0];
    
    // Allocate memory for measurements
    *measurements = malloc(*num_measurements * sizeof(Measurement));
    if (!*measurements) {
        fprintf(stderr, "Error: Could not allocate memory for measurements\n");
        H5Sclose(dataspace_id);
        H5Dclose(dset_id);
        H5Fclose(file_id);
        return -1;
    }
    
    // Create memory datatype
    hid_t memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(Measurement));
    H5Tinsert(memtype_id, "sgn", HOFFSET(Measurement, sgn), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "H", HOFFSET(Measurement, H), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "H2", HOFFSET(Measurement, H2), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "Z_mag", HOFFSET(Measurement, Z_mag), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "q", HOFFSET(Measurement, q), H5T_NATIVE_INT);
    
    // Read data
    if (H5Dread(dset_id, memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, *measurements) < 0) {
        fprintf(stderr, "Error: Could not read data from file '%s'\n", filename);
        free(*measurements);
        *measurements = NULL;
        H5Tclose(memtype_id);
        H5Sclose(dataspace_id);
        H5Dclose(dset_id);
        H5Fclose(file_id);
        return -1;
    }
    
    // Cleanup
    H5Tclose(memtype_id);
    H5Sclose(dataspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    
    return 0;
}

// Combine MPI data from all ranks
int processing_combine_mpi_data(int num_ranks, const char* output_filename) {
    Measurement** all_measurements = malloc(num_ranks * sizeof(Measurement*));
    hsize_t* num_measurements = calloc(num_ranks, sizeof(hsize_t));  // Use calloc to initialize to 0
    
    if (!all_measurements || !num_measurements) {
        fprintf(stderr, "Error: Could not allocate memory for data arrays\n");
        free(all_measurements);
        free(num_measurements);
        return -1;
    }
    
    // Read data from all ranks
    for (int rank = 0; rank < num_ranks; rank++) {
        char filename[256];
        processing_get_rank_filename(filename, rank, sizeof(filename));
        
        if (read_rank_data(filename, &all_measurements[rank], &num_measurements[rank]) < 0) {
            fprintf(stderr, "Error: Could not read data from rank %d\n", rank);
            // Clean up already allocated memory
            for (int j = 0; j < rank; j++) {
                free(all_measurements[j]);
            }
            free(all_measurements);
            free(num_measurements);
            return -1;
        }
    }
    
    // Verify all ranks have the same number of measurements
    hsize_t common_num_measurements = num_measurements[0];
    for (int rank = 1; rank < num_ranks; rank++) {
        if (num_measurements[rank] != common_num_measurements) {
            fprintf(stderr, "Error: Rank %d has %llu measurements, but rank 0 has %llu\n", 
                    rank, (unsigned long long)num_measurements[rank], 
                    (unsigned long long)common_num_measurements);
            // Clean up
            for (int j = 0; j < num_ranks; j++) {
                free(all_measurements[j]);
            }
            free(all_measurements);
            free(num_measurements);
            return -1;
        }
    }
    
    // Allocate memory for combined measurements
    Measurement* combined_measurements = malloc(common_num_measurements * sizeof(Measurement));
    if (!combined_measurements) {
        fprintf(stderr, "Error: Could not allocate memory for combined measurements\n");
        // Clean up
        for (int j = 0; j < num_ranks; j++) {
            free(all_measurements[j]);
        }
        free(all_measurements);
        free(num_measurements);
        return -1;
    }
    
    // Combine data by averaging measurements at each step
    for (hsize_t i = 0; i < common_num_measurements; i++) {
        // Initialize combined measurement
        combined_measurements[i].sgn = 0.0;
        combined_measurements[i].H = 0.0;
        combined_measurements[i].H2 = 0.0;
        combined_measurements[i].Z_mag = 0.0;
        combined_measurements[i].q = 0;
        
        // Sum measurements from all ranks
        for (int rank = 0; rank < num_ranks; rank++) {
            combined_measurements[i].sgn += all_measurements[rank][i].sgn;
            combined_measurements[i].H += all_measurements[rank][i].H;
            combined_measurements[i].H2 += all_measurements[rank][i].H2;
            combined_measurements[i].Z_mag += all_measurements[rank][i].Z_mag;
            combined_measurements[i].q += all_measurements[rank][i].q;
        }
        
        // Average the measurements
        combined_measurements[i].sgn /= num_ranks;
        combined_measurements[i].H /= num_ranks;
        combined_measurements[i].H2 /= num_ranks;
        combined_measurements[i].Z_mag /= num_ranks;
        combined_measurements[i].q /= num_ranks; // Integer division
    }
    
    // Create output HDF5 file
    hid_t file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error: Could not create output HDF5 file '%s'\n", output_filename);
        free(combined_measurements);
        for (int j = 0; j < num_ranks; j++) {
            free(all_measurements[j]);
        }
        free(all_measurements);
        free(num_measurements);
        return -1;
    }
    
    // For now, we'll create a simple dataset without attributes
    // In a full implementation, you'd want to copy attributes from one of the input files
    
    // Create dataspace
    hsize_t dims[1] = {common_num_measurements};
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    // Create datatype
    hid_t memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(Measurement));
    H5Tinsert(memtype_id, "sgn", HOFFSET(Measurement, sgn), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "H", HOFFSET(Measurement, H), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "H2", HOFFSET(Measurement, H2), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "Z_mag", HOFFSET(Measurement, Z_mag), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "q", HOFFSET(Measurement, q), H5T_NATIVE_INT);
    
    // Create dataset
    hid_t dset_id = H5Dcreate2(file_id, "/measurements", memtype_id, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write data
    if (H5Dwrite(dset_id, memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, combined_measurements) < 0) {
        fprintf(stderr, "Error: Could not write combined data to output file\n");
        H5Dclose(dset_id);
        H5Tclose(memtype_id);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        free(combined_measurements);
        for (int j = 0; j < num_ranks; j++) {
            free(all_measurements[j]);
        }
        free(all_measurements);
        free(num_measurements);
        return -1;
    }
    
    // Cleanup
    H5Dclose(dset_id);
    H5Tclose(memtype_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);
    
    free(combined_measurements);
    for (int j = 0; j < num_ranks; j++) {
        free(all_measurements[j]);
    }
    free(all_measurements);
    free(num_measurements);
    
    printf("Successfully combined data from %d ranks into '%s'\n", num_ranks, output_filename);
    return 0;
}

// Clean up temporary MPI data files
int processing_cleanup_mpi_data(int num_ranks) {
    // Remove individual rank files
    for (int rank = 0; rank < num_ranks; rank++) {
        char filename[256];
        processing_get_rank_filename(filename, rank, sizeof(filename));
        remove(filename);
    }
    
    // Remove directory
    rmdir("raw_mpi_data");
    
    return 0;
}