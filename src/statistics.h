// File: src/statistics.h
// Purpose: Declares functions for statistical analysis of raw QMC data.

#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdio.h> // Needed for the FILE type
#include "datatypes.h"

/**
 * @brief Performs binning analysis on raw data arrays and writes the
 *        formatted results to the provided file stream.
 * @param output_stream An open, writable file stream (e.g., a file opened with "w").
 * @param sgn_data Array of measured signs.
 * @param H_data Array of measured <H> values.
 * @param H2_data Array of measured <H^2> values.
 * @param q_data Array of measured q values.
 * @param num_points The number of data points in the arrays.
 * @param params The simulation parameters (for NBINS, BETA, etc.).
 */
void perform_analysis_and_write_results(
    FILE* output_stream,
    const double* sgn_data,
    const double* H_data,
    const double* H2_data,
    const int* q_data,
    long long num_points,
    const SimParams* params
);

#endif // STATISTICS_H