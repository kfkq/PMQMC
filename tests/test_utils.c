// File: tests/test_utils.c

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"

#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fprintf(stderr, "\n--- TEST FAILED: %s ---\n", message); \
        exit(1); \
    } else { \
        printf("PASSED: %s\n", message); \
    }

#define MAX_BINS 100

// Helper function to calculate chi-square statistic for uniform distribution
double chi_square_test(int* observed, int expected, int bins) {
    double chi_square = 0.0;
    for (int i = 0; i < bins; i++) {
        double diff = observed[i] - expected;
        chi_square += (diff * diff) / expected;
    }
    return chi_square;
}

int main() {
    printf("\n--- Running Tests for Utils Module ---\n");

    // --- Test 1: Deterministic Seeding ---
    // Seeding the RNG with the same seed should produce the same sequence of numbers.
    rng_init(12345);
    int first_val = rng_uniform_int(1000);
    double second_val = rng_uniform_double();

    rng_init(12345); // Re-seed with the same value
    TEST_ASSERT(rng_uniform_int(1000) == first_val, "Seeding is deterministic (integer)");
    TEST_ASSERT(rng_uniform_double() == second_val, "Seeding is deterministic (double)");
    printf("Deterministic seeding tests seem OK.\n\n");

    // --- Test 2: Range Checks ---
    printf("--- Test 2: Range Checks ---\n");
    rng_init(time(NULL)); // Use a different seed for range checks

    // Check double range [0.0, 1.0)
    for (int i = 0; i < 50000; ++i) {
        double val = rng_uniform_double();
        if (val < 0.0 || val >= 1.0) {
            // Fail explicitly if out of range
            TEST_ASSERT(0, "Double value was out of the expected [0.0, 1.0) range");
        }
    }
    TEST_ASSERT(1, "All doubles were in the expected [0.0, 1.0) range");

    // Check integer range [0, max-1]
    int max_val = 100;
    for (int i = 0; i < 50000; ++i) {
        int val = rng_uniform_int(max_val);
        if (val < 0 || val >= max_val) {
            // Fail explicitly if out of range
            TEST_ASSERT(0, "Integer value was out of the expected [0, max) range");
        }
    }
    TEST_ASSERT(1, "All integers were in the expected [0, max) range");
    printf("Range check tests seem OK.\n\n");

    // Re-seed with a fixed value for deterministic statistical tests
    rng_init(42);

    // --- Test 3: Uniform Distribution of Integers ---
    printf("--- Test 3: Uniform Distribution of Integers ---\n");
    const int bins = 10;
    const int samples = 100000;
    int observed[MAX_BINS] = {0};
    
    // Generate random integers and count occurrences in each bin
    for (int i = 0; i < samples; ++i) {
        int val = rng_uniform_int(bins);
        observed[val]++;
    }
    
    // Calculate expected frequency
    int expected = samples / bins;
    
    // Calculate chi-square statistic
    double chi_square = chi_square_test(observed, expected, bins);
    
    // For 9 degrees of freedom (bins-1), critical value at 0.05 significance level is ~16.92
    double critical_value = 16.92;
    
    printf("Chi-square statistic: %.2f (critical value: %.2f)\n", chi_square, critical_value);
    TEST_ASSERT(chi_square < critical_value, "Integer distribution appears uniform (chi-square test)");
    printf("Integer distribution test seems OK.\n\n");

    // --- Test 4: Uniform Distribution of Doubles ---
    printf("--- Test 4: Uniform Distribution of Doubles ---\n");
    const int double_bins = 10;
    int double_observed[MAX_BINS] = {0};
    
    // Generate random doubles and count occurrences in each bin
    for (int i = 0; i < samples; ++i) {
        double val = rng_uniform_double();
        int bin = (int)(val * double_bins);
        // Handle edge case where val is exactly 1.0 (shouldn't happen but just in case)
        if (bin >= double_bins) bin = double_bins - 1;
        double_observed[bin]++;
    }
    
    // Calculate chi-square statistic
    double double_chi_square = chi_square_test(double_observed, expected, double_bins);
    
    printf("Chi-square statistic: %.2f (critical value: %.2f)\n", double_chi_square, critical_value);
    TEST_ASSERT(double_chi_square < critical_value, "Double distribution appears uniform (chi-square test)");
    printf("Double distribution test seems OK.\n\n");

    // --- Test 5: Mean and Variance ---
    printf("--- Test 5: Mean and Variance of Doubles ---\n");
    const int stat_samples = 100000;
    double sum = 0.0;
    double sum_sq = 0.0;
    
    // Generate samples and calculate sum and sum of squares
    for (int i = 0; i < stat_samples; ++i) {
        double val = rng_uniform_double();
        sum += val;
        sum_sq += val * val;
    }
    
    // Calculate mean and variance
    double mean = sum / stat_samples;
    double variance = (sum_sq / stat_samples) - (mean * mean);
    
    // For uniform distribution [0,1): expected mean = 0.5, expected variance = 1/12 ≈ 0.0833
    double expected_mean = 0.5;
    double expected_variance = 1.0/12.0;
    
    // Allow some tolerance for statistical variation
    double mean_tolerance = 0.01;
    double variance_tolerance = 0.001;
    
    printf("Sample mean: %.4f (expected: %.4f)\n", mean, expected_mean);
    printf("Sample variance: %.4f (expected: %.4f)\n", variance, expected_variance);
    
    TEST_ASSERT(fabs(mean - expected_mean) < mean_tolerance, "Mean is close to expected value");
    TEST_ASSERT(fabs(variance - expected_variance) < variance_tolerance, "Variance is close to expected value");
    printf("Statistical properties test seems OK.\n\n");

    // --- Test 6: Independence of Successive Values ---
    printf("--- Test 6: Independence of Successive Values ---\n");
    const int corr_samples = 50000;
    double x_sum = 0.0, y_sum = 0.0, xy_sum = 0.0;
    double x_sq_sum = 0.0, y_sq_sum = 0.0;
    
    // Generate pairs of successive values
    for (int i = 0; i < corr_samples; ++i) {
        double x = rng_uniform_double();
        double y = rng_uniform_double();
        
        x_sum += x;
        y_sum += y;
        xy_sum += x * y;
        x_sq_sum += x * x;
        y_sq_sum += y * y;
    }
    
    // Calculate correlation coefficient
    double x_mean = x_sum / corr_samples;
    double y_mean = y_sum / corr_samples;
    double numerator = xy_sum - corr_samples * x_mean * y_mean;
    double denominator = sqrt((x_sq_sum - corr_samples * x_mean * x_mean) * (y_sq_sum - corr_samples * y_mean * y_mean));
    double correlation = numerator / denominator;
    
    // For independent variables, correlation should be close to 0
    double correlation_tolerance = 0.02;
    
    printf("Correlation coefficient: %.4f (should be close to 0)\n", correlation);
    TEST_ASSERT(fabs(correlation) < correlation_tolerance, "Successive values appear independent");
    printf("Independence test seems OK.\n\n");

    printf("--- All Utils tests passed. ---\n");
    return 0;
}