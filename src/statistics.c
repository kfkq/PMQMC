// File: src/statistics.c
// Purpose: Implements binning analysis on raw QMC data arrays and writes
//          the results to a file stream.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "statistics.h"

// --- Static Helper Functions for Statistical Calculations ---

static double array_mean(const double* arr, int n) {
    if (n == 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += arr[i];
    }
    return sum / n;
}

static double array_variance(const double* arr, int n, double mean) {
    if (n <= 1) return 0.0;
    double sum_sq_diff = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_sq_diff += (arr[i] - mean) * (arr[i] - mean);
    }
    return sum_sq_diff / (n * (n - 1.0));
}

static double array_covariance(const double* arr1, const double* arr2, int n, double mean1, double mean2) {
    if (n <= 1) return 0.0;
    double sum_prod_diff = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_prod_diff += (arr1[i] - mean1) * (arr2[i] - mean2);
    }
    return sum_prod_diff / (n * (n - 1.0));
}


// --- Main Public Analysis Function ---

void perform_analysis_and_write_results(
    FILE* output_stream,
    const double* sgn_data,
    const double* H_data,
    const double* H2_data,
    const int* q_data,
    long long num_points,
    const SimParams* params
) {
    if (params->NBINS <= 1) {
        fprintf(stderr, "Error: NBINS must be greater than 1 for statistical analysis.\n");
        return;
    }
    long long bin_length = num_points / params->NBINS;
    if (bin_length == 0) {
        fprintf(stderr, "Warning: Not enough data points (%lld) for %d bins. Analysis cannot proceed.\n", num_points, params->NBINS);
        return;
    }

    // --- 1. Create Bin Averages from Raw Data ---
    double* bin_mean_sgn = calloc(params->NBINS, sizeof(double));
    double* bin_mean_H_sgn = calloc(params->NBINS, sizeof(double));
    double* bin_mean_H2_sgn = calloc(params->NBINS, sizeof(double));
    double* bin_mean_q = calloc(params->NBINS, sizeof(double));
    int max_q = 0;

    if (!bin_mean_sgn || !bin_mean_H_sgn || !bin_mean_H2_sgn || !bin_mean_q) {
        fprintf(stderr, "Error: Failed to allocate memory for bin average arrays.\n");
        free(bin_mean_sgn); free(bin_mean_H_sgn); free(bin_mean_H2_sgn); free(bin_mean_q);
        return;
    }

    for (int i = 0; i < params->NBINS; ++i) {
        double sgn_sum = 0, H_sum = 0, H2_sum = 0, q_sum = 0;
        for (long long j = 0; j < bin_length; ++j) {
            long long idx = i * bin_length + j;
            sgn_sum += sgn_data[idx];
            H_sum += H_data[idx] * sgn_data[idx];
            H2_sum += H2_data[idx] * sgn_data[idx];
            q_sum += (double)q_data[idx];
            if (q_data[idx] > max_q) max_q = q_data[idx];
        }
        bin_mean_sgn[i] = sgn_sum / bin_length;
        bin_mean_H_sgn[i] = H_sum / bin_length;
        bin_mean_H2_sgn[i] = H2_sum / bin_length;
        bin_mean_q[i] = q_sum / bin_length;
    }

    // --- 2. Calculate Final Statistics and Write to File ---
    fprintf(output_stream, "--- Final Results ---\n");
    fprintf(output_stream, "Analyzed %lld data points, skipping %lld.\n", num_points, params->SKIP_MEASUREMENTS);
    fprintf(output_stream, "Number of Bins: %d\n", params->NBINS);
    fprintf(output_stream, "---------------------\n\n");

    double mean_sgn = array_mean(bin_mean_sgn, params->NBINS);
    double var_sgn = array_variance(bin_mean_sgn, params->NBINS, mean_sgn);
    double stdev_sgn = sqrt(var_sgn);

    double mean_H_sgn = array_mean(bin_mean_H_sgn, params->NBINS);
    double var_H_sgn = array_variance(bin_mean_H_sgn, params->NBINS, mean_H_sgn);
    double cov_H = array_covariance(bin_mean_sgn, bin_mean_H_sgn, params->NBINS, mean_sgn, mean_H_sgn);

    double mean_H2_sgn = array_mean(bin_mean_H2_sgn, params->NBINS);
    double var_H2_sgn = array_variance(bin_mean_H2_sgn, params->NBINS, mean_H2_sgn);
    double cov_H2 = array_covariance(bin_mean_sgn, bin_mean_H2_sgn, params->NBINS, mean_sgn, mean_H2_sgn);

    fprintf(output_stream, "Average Sign: %.8f +/- %.8f\n\n", mean_sgn, stdev_sgn);

    // Final Observable <H>
    double final_H = (mean_H_sgn / mean_sgn) * (1.0 + var_sgn / (mean_sgn * mean_sgn)) - cov_H / (mean_sgn * mean_sgn);
    double stdev_H_term1 = (mean_H_sgn != 0) ? var_H_sgn / (mean_H_sgn * mean_H_sgn) : 0;
    double stdev_H_term2 = (mean_sgn != 0) ? var_sgn / (mean_sgn * mean_sgn) : 0;
    double stdev_H_term3 = (mean_H_sgn != 0 && mean_sgn != 0) ? 2.0 * cov_H / (mean_H_sgn * mean_sgn) : 0;
    double stdev_H = fabs(final_H) * sqrt(fmax(0.0, stdev_H_term1 + stdev_H_term2 - stdev_H_term3));

    fprintf(output_stream, "Observable <H>:\n");
    fprintf(output_stream, "  Mean      = %.8f\n", final_H);
    fprintf(output_stream, "  Std. Dev. = %.8f\n\n", stdev_H);

    // Final Observable <H^2>
    double final_H2 = (mean_H2_sgn / mean_sgn) * (1.0 + var_sgn / (mean_sgn * mean_sgn)) - cov_H2 / (mean_sgn * mean_sgn);
    double stdev_H2_term1 = (mean_H2_sgn != 0) ? var_H2_sgn / (mean_H2_sgn * mean_H2_sgn) : 0;
    double stdev_H2_term2 = (mean_sgn != 0) ? var_sgn / (mean_sgn * mean_sgn) : 0;
    double stdev_H2_term3 = (mean_H2_sgn != 0 && mean_sgn != 0) ? 2.0 * cov_H2 / (mean_H2_sgn * mean_sgn) : 0;
    double stdev_H2 = fabs(final_H2) * sqrt(fmax(0.0, stdev_H2_term1 + stdev_H2_term2 - stdev_H2_term3));

    fprintf(output_stream, "Observable <H^2>:\n");
    fprintf(output_stream, "  Mean      = %.8f\n", final_H2);
    fprintf(output_stream, "  Std. Dev. = %.8f\n\n", stdev_H2);

    // Other Observables
    double final_q = array_mean(bin_mean_q, params->NBINS);
    double var_q = array_variance(bin_mean_q, params->NBINS, final_q);
    double stdev_q = sqrt(var_q);

    fprintf(output_stream, "Observable <Mean q>:\n");
    fprintf(output_stream, "  Mean      = %.8f\n", final_q);
    fprintf(output_stream, "  Std. Dev. = %.8f\n\n", stdev_q);

    fprintf(output_stream, "Observable <Max q>:\n");
    fprintf(output_stream, "  Value     = %d\n\n", max_q);

    // Derived Observable: Specific Heat
    if (fabs(mean_sgn) > 1e-9) {
        double cv = (final_H2 - final_H * final_H) * params->BETA * params->BETA;
        fprintf(output_stream, "Derived Observable (Specific Heat C_v):\n");
        fprintf(output_stream, "  Value     = %.8f\n", cv);
    }
    fprintf(output_stream, "---------------------\n");

    // --- 3. Cleanup ---
    free(bin_mean_sgn);
    free(bin_mean_H_sgn);
    free(bin_mean_H2_sgn);
    free(bin_mean_q);
}