// File: src/divdiff.h
// Purpose: Defines the API for the *optimized* O(q) divided differences
//          calculation engine.

#ifndef DIVDIFF_H
#define DIVDIFF_H

#include "datatypes.h"

// --- Custom High-Precision Floating Point Type ---
typedef struct {
    double mantissa;
    int exponent;
} ExExFloat;


// --- Optimized Divided Differences Calculator State ---
typedef struct {
    int current_len;
    int max_len;
    int s_max;

    // Internal state for O(q) updates
    double mu;
    int s;
    ExExFloat exp_mu;
    double* energies;
    ExExFloat* h;
    ExExFloat* ddd;
    
    // Final results are stored here
    ExExFloat* results;
} DivDiff;


// --- Global Initialization/Cleanup ---
void divdiff_global_init();
void divdiff_global_cleanup();


// --- DivDiff Object Management ---
DivDiff* divdiff_create(int max_q, int s_max);
void divdiff_free(DivDiff* dd);


// --- Core O(q) Calculation Functions ---
void divdiff_add_element(DivDiff* dd, double new_energy);
void divdiff_remove_last_element(DivDiff* dd);


// --- ExExFloat Helper Functions ---
ExExFloat exex_from_double(double d);
double exex_to_double(ExExFloat f);
ExExFloat exex_add(ExExFloat a, ExExFloat b);
ExExFloat exex_subtract(ExExFloat a, ExExFloat b);
ExExFloat exex_multiply(ExExFloat a, ExExFloat b);
ExExFloat exex_divide(ExExFloat a, ExExFloat b);
ExExFloat exex_multiply_double(ExExFloat a, double b);

#endif // DIVDIFF_H