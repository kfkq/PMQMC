// File: src/processing.h
// Purpose: Post-processing functions for MPI-QMC simulation data

#ifndef PROCESSING_H
#define PROCESSING_H

#include <hdf5.h>
#include "datatypes.h"
#include "io.h"

// Function to initialize MPI data collection
int processing_init_mpi_data_directory(void);

// Function to generate rank-specific filename
void processing_get_rank_filename(char* filename, int rank, size_t filename_size);

// Function to read and combine MPI data from all ranks
int processing_combine_mpi_data(int num_ranks, const char* output_filename);

// Function to clean up temporary MPI data files
int processing_cleanup_mpi_data(int num_ranks);

#endif // PROCESSING_H