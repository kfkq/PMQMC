// File: src/analyzer.c
// Purpose: Reads raw QMC data and performs statistical analysis.
// EDIT: Now automatically finds and reads "hamiltonian.in" and "raw.dat"
//       from the current directory without command-line arguments.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "datatypes.h"
#include "hamiltonian.h"
#include "statistics.h"

// --- Define the required input and output filenames ---
const char* INPUT_FILENAME = "hamiltonian.in";
const char* RAW_DATA_FILENAME = "raw.dat";
const char* RESULTS_FILENAME = "results.dat";

int main(void) { // No longer needs argc, argv
    // 1. Check for and load simulation parameters
    FILE* check_input_file = fopen(INPUT_FILENAME, "r");
    if (!check_input_file) {
        fprintf(stderr, "Error: Could not find or open the required parameter file '%s'.\n", INPUT_FILENAME);
        return 1;
    }
    fclose(check_input_file);

    SimParams params = {0};
    Hamiltonian* h = hamiltonian_create_and_load(INPUT_FILENAME, &params);
    if (!h) return 1;
    hamiltonian_free(h, &params); // We only needed the params

    // 2. Check for and open the raw data file
    FILE* raw_in_file = fopen(RAW_DATA_FILENAME, "r");
    if (!raw_in_file) {
        fprintf(stderr, "Error: Could not find or open the required data file '%s'.\n", RAW_DATA_FILENAME);
        return 1;
    }

    // 3. Count the number of data points
    printf("Reading raw data from '%s'...\n", RAW_DATA_FILENAME);
    long long total_lines = 0;
    char buffer[256];
    fgets(buffer, sizeof(buffer), raw_in_file); // Skip header
    while (fgets(buffer, sizeof(buffer), raw_in_file)) {
        total_lines++;
    }

    if (total_lines <= params.SKIP_MEASUREMENTS) {
        fprintf(stderr, "Error: SKIP_MEASUREMENTS (%lld) is >= total measurements (%lld).\n",
                params.SKIP_MEASUREMENTS, total_lines);
        fclose(raw_in_file);
        return 1;
    }
    long long points_to_analyze = total_lines - params.SKIP_MEASUREMENTS;

    // 4. Allocate memory and read the data
    double* sgn_data = malloc(points_to_analyze * sizeof(double));
    double* H_data = malloc(points_to_analyze * sizeof(double));
    double* H2_data = malloc(points_to_analyze * sizeof(double));
    int* q_data = malloc(points_to_analyze * sizeof(int));

    if (!sgn_data || !H_data || !H2_data || !q_data) {
        fprintf(stderr, "Error: Failed to allocate memory for data arrays.\n");
        fclose(raw_in_file);
        return 1;
    }

    rewind(raw_in_file);
    fgets(buffer, sizeof(buffer), raw_in_file); // Skip header
    for (long long i = 0; i < params.SKIP_MEASUREMENTS; ++i) {
        fgets(buffer, sizeof(buffer), raw_in_file);
    }

    // --- Read data with a progress indicator ---
    const long long updates_to_print = 100;
    const long long print_interval = (points_to_analyze > updates_to_print) ? (points_to_analyze / updates_to_print) : 1;

    for (long long i = 0; i < points_to_analyze; ++i) {
        if (fscanf(raw_in_file, "%lf %lf %lf %d", &sgn_data[i], &H_data[i], &H2_data[i], &q_data[i]) != 4) {
            fprintf(stderr, "\nError reading data at line %lld.\n", i + params.SKIP_MEASUREMENTS + 1);
            points_to_analyze = i;
            break;
        }
        if (i > 0 && i % print_interval == 0) {
            double percent_done = (double)i / points_to_analyze * 100.0;
            printf("\rReading data: [%.2f%%]", percent_done);
            fflush(stdout);
        }
    }
    printf("\rReading data: [100.00%%]\n");
    fclose(raw_in_file);

    // 5. Open the output file and perform analysis
    FILE* results_file = fopen(RESULTS_FILENAME, "w");
    if (!results_file) {
        perror("Error opening results.dat for writing");
        free(sgn_data); free(H_data); free(H2_data); free(q_data);
        return 1;
    }

    printf("Analyzing %lld data points (skipped %lld initial measurements).\n",
           points_to_analyze, params.SKIP_MEASUREMENTS);
    
    perform_analysis_and_write_results(results_file, sgn_data, H_data, H2_data, q_data, points_to_analyze, &params);

    fclose(results_file);
    printf("Analysis complete. Results written to %s\n", RESULTS_FILENAME);

    // 6. Cleanup
    free(sgn_data);
    free(H_data);
    free(H2_data);
    free(q_data);

    return 0;
}