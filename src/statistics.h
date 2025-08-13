// File: src/statistics.h
// Purpose: Declares structures and functions for statistical analysis (binning).

#ifndef STATISTICS_H
#define STATISTICS_H

#include "datatypes.h"

// Holds all data needed for binning analysis.
typedef struct {
    int nbins;
    long long bin_length;
    long long measurement_step;

    // Accumulators for the current bin
    double in_bin_sum_sgn;
    double in_bin_sum_H;
    double in_bin_sum_H2;
    double in_bin_sum_q;

    // Arrays of averages for each completed bin
    double* bin_mean_sgn;
    double* bin_mean_H;
    double* bin_mean_H2;
    double* bin_mean_q;

    // Other statistics
    int max_q;

} Stats;

/**
 * @brief Creates and initializes the statistics tracking object.
 * @param params The simulation parameters.
 * @return A pointer to the newly created Stats object, or NULL on failure.
 */
Stats* stats_create(const SimParams* params);

/**
 * @brief Frees all memory associated with the Stats object.
 * @param stats The Stats object to free.
 */
void stats_free(Stats* stats);

/**
 * @brief Accumulates new measurements into the current bin.
 *        Handles bin completion automatically.
 * @param stats The statistics object.
 * @param sgn The sign of the current configuration's weight.
 * @param H_val The measured value of <H>.
 * @param H2_val The measured value of <H^2>.
 */
void stats_accumulate(Stats* stats, double sgn, double H_val, double H2_val, int q_val);

/**
 * @brief Performs the final statistical analysis and prints the results.
 * @param stats The statistics object.
 * @param params The simulation parameters.
 */
void stats_finalize_and_print(const Stats* stats, const SimParams* params);

#endif // STATISTICS_H