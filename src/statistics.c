// File: src/statistics.c
// Purpose: Implements binning analysis for QMC results.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "statistics.h"

Stats* stats_create(const SimParams* params) {
    if (params->NBINS <= 1) {
        fprintf(stderr, "Error: NBINS must be greater than 1 for statistical analysis.\n");
        return NULL;
    }
    Stats* stats = (Stats*)calloc(1, sizeof(Stats));
    if (!stats) return NULL;

    stats->nbins = params->NBINS;
    stats->bin_length = params->STEPS / params->STEPS_PER_MEASUREMENT / params->NBINS;
    if (stats->bin_length == 0) {
        fprintf(stderr, "Warning: Total measurements are less than nbins. Statistics may be unreliable.\n");
        stats->bin_length = 1;
    }

    stats->bin_mean_sgn = (double*)calloc(stats->nbins, sizeof(double));
    stats->bin_mean_H = (double*)calloc(stats->nbins, sizeof(double));
    stats->bin_mean_H2 = (double*)calloc(stats->nbins, sizeof(double));

    if (!stats->bin_mean_sgn || !stats->bin_mean_H || !stats->bin_mean_H2) {
        stats_free(stats);
        return NULL;
    }
    return stats;
}

void stats_free(Stats* stats) {
    if (stats) {
        free(stats->bin_mean_sgn);
        free(stats->bin_mean_H);
        free(stats->bin_mean_H2);
        free(stats);
    }
}

void stats_accumulate(Stats* stats, double sgn, double H_val, double H2_val) {
    stats->in_bin_sum_sgn += sgn;
    stats->in_bin_sum_H += H_val * sgn;
    stats->in_bin_sum_H2 += H2_val * sgn;
    stats->measurement_step++;

    if (stats->measurement_step > 0 && stats->measurement_step % stats->bin_length == 0) {
        int bin_idx = (stats->measurement_step / stats->bin_length) - 1;
        if (bin_idx < stats->nbins) {
            stats->bin_mean_sgn[bin_idx] = stats->in_bin_sum_sgn / stats->bin_length;
            stats->bin_mean_H[bin_idx] = stats->in_bin_sum_H / stats->bin_length;
            stats->bin_mean_H2[bin_idx] = stats->in_bin_sum_H2 / stats->bin_length;
        }
        // Reset for the next bin
        stats->in_bin_sum_sgn = 0.0;
        stats->in_bin_sum_H = 0.0;
        stats->in_bin_sum_H2 = 0.0;
    }
}

// --- Final Analysis Helper Functions (static) ---

static double array_mean(const double* arr, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += arr[i];
    return sum / n;
}

static double array_variance(const double* arr, int n, double mean) {
    double sum_sq_diff = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_sq_diff += (arr[i] - mean) * (arr[i] - mean);
    }
    return sum_sq_diff / (n * (n - 1.0)); // Variance of the mean
}

static double array_covariance(const double* arr1, const double* arr2, int n, double mean1, double mean2) {
    double sum_prod_diff = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_prod_diff += (arr1[i] - mean1) * (arr2[i] - mean2);
    }
    return sum_prod_diff / (n * (n - 1.0)); // Covariance of the means
}

void stats_finalize_and_print(const Stats* stats, const SimParams* params) {
    printf("\n--- Final Results ---\n");

    // --- Basic Statistics from Bins ---
    double mean_sgn = array_mean(stats->bin_mean_sgn, stats->nbins);
    double var_sgn = array_variance(stats->bin_mean_sgn, stats->nbins, mean_sgn);
    double stdev_sgn = sqrt(var_sgn);

    double mean_H_sgn = array_mean(stats->bin_mean_H, stats->nbins);
    double var_H_sgn = array_variance(stats->bin_mean_H, stats->nbins, mean_H_sgn);
    double cov_H = array_covariance(stats->bin_mean_sgn, stats->bin_mean_H, stats->nbins, mean_sgn, mean_H_sgn);

    double mean_H2_sgn = array_mean(stats->bin_mean_H2, stats->nbins);
    double var_H2_sgn = array_variance(stats->bin_mean_H2, stats->nbins, mean_H2_sgn);
    double cov_H2 = array_covariance(stats->bin_mean_sgn, stats->bin_mean_H2, stats->nbins, mean_sgn, mean_H2_sgn);

    printf("Average Sign: %.8f +/- %.8f\n\n", mean_sgn, stdev_sgn);

    // --- Final Observable <H> ---
    // Formula from Appendix B of the paper (and datasummary.cpp)
    double final_H = (mean_H_sgn / mean_sgn) * (1.0 + var_sgn / (mean_sgn * mean_sgn)) - cov_H / (mean_sgn * mean_sgn);
    
    double stdev_H_term1 = var_H_sgn / (mean_H_sgn * mean_H_sgn);
    double stdev_H_term2 = var_sgn / (mean_sgn * mean_sgn);
    double stdev_H_term3 = 2.0 * cov_H / (mean_H_sgn * mean_sgn);
    double stdev_H = fabs(final_H) * sqrt(fmax(0.0, stdev_H_term1 + stdev_H_term2 - stdev_H_term3));

    printf("Observable <H>:\n");
    printf("  Mean      = %.8f\n", final_H);
    printf("  Std. Dev. = %.8f\n\n", stdev_H);

    // --- Final Observable <H^2> ---
    double final_H2 = (mean_H2_sgn / mean_sgn) * (1.0 + var_sgn / (mean_sgn * mean_sgn)) - cov_H2 / (mean_sgn * mean_sgn);
    
    double stdev_H2_term1 = var_H2_sgn / (mean_H2_sgn * mean_H2_sgn);
    double stdev_H2_term2 = var_sgn / (mean_sgn * mean_sgn);
    double stdev_H2_term3 = 2.0 * cov_H2 / (mean_H2_sgn * mean_sgn);
    double stdev_H2 = fabs(final_H2) * sqrt(fmax(0.0, stdev_H2_term1 + stdev_H2_term2 - stdev_H2_term3));

    printf("Observable <H^2>:\n");
    printf("  Mean      = %.8f\n", final_H2);
    printf("  Std. Dev. = %.8f\n\n", stdev_H2);

    // --- Derived Observable: Specific Heat ---
    if (final_H != 0 && final_H2 != 0) {
        double cv = (final_H2 - final_H * final_H) * params->BETA * params->BETA;
        // Note: Calculating the error bar for a derived quantity like this rigorously
        // requires the jackknife method (as in the C++ code) or error propagation formulas.
        // For now, we just print the value.
        printf("Derived Observable (Specific Heat C_v):\n");
        printf("  Value     = %.8f\n", cv);
    }
    printf("---------------------\n");
}