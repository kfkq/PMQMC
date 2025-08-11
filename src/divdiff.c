// File: src/divdiff.c
// Purpose: Implements the *optimized* O(q) divided differences calculation engine.
// VERSION: Corrected via direct translation from the C++ reference.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "divdiff.h"

// --- Static Globals & Helpers ---
static double* invPowersOf2 = NULL;
static const int MAX_EXP_DIFF = 1024;
static const int EXTRA_LEN = 10; // For buffer in h array

// Forward declaration for internal functions
static void divdiff_recalculate_all(DivDiff* dd, int force_s, double force_mu);
static void divdiff_add_element_internal(DivDiff* dd, double new_energy, int force_s, double force_mu);


// Helper to calculate the mean of an array
static double calculate_mean(const double* arr, int len) {
    if (len <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < len; ++i) sum += arr[i];
    return sum / len;
}

// Helper to find the max absolute difference in an array
static double max_abs_diff(const double* arr, int len) {
    if (len <= 1) return 0.0;
    double min_val = arr[0], max_val = arr[0];
    for (int i = 1; i < len; ++i) {
        if (arr[i] < min_val) min_val = arr[i];
        if (arr[i] > max_val) max_val = arr[i];
    }
    return max_val - min_val;
}

// Helper to check if the scaling factor 's' needs to be updated
static int s_changed(const DivDiff* dd) {
    if (dd->current_len <= 1) return 0;
    return fabs(dd->energies[dd->current_len - 1] - dd->mu) / 3.5 > dd->s;
}

static void exex_normalize(ExExFloat* f) {
    if (f->mantissa == 0.0) { f->exponent = 0; return; }
    int exp_change;
    f->mantissa = frexp(f->mantissa, &exp_change);
    f->exponent += exp_change;
}

// Helper to initialize exp_mu from mu
static void exex_init_expmu(ExExFloat* f, double mu) {
    double e = mu * 1.4426950408889634; // mu / log(2)
    f->exponent = (int)ceil(e);
    f->mantissa = pow(2.0, e - ceil(e));
}

// --- ExExFloat Implementation (Correct) ---
ExExFloat exex_from_double(double d) { ExExFloat f; f.mantissa = d; f.exponent = 0; exex_normalize(&f); return f; }
double exex_to_double(ExExFloat f) { return ldexp(f.mantissa, f.exponent); }
ExExFloat exex_multiply_double(ExExFloat a, double b) { ExExFloat res = a; res.mantissa *= b; exex_normalize(&res); return res; }
ExExFloat exex_add(ExExFloat a, ExExFloat b) {
    ExExFloat res;
    if (a.exponent >= b.exponent) {
        int diff = a.exponent - b.exponent;
        if (diff >= MAX_EXP_DIFF) return a;
        res.mantissa = a.mantissa + b.mantissa * invPowersOf2[diff];
        res.exponent = a.exponent;
    } else {
        int diff = b.exponent - a.exponent;
        if (diff >= MAX_EXP_DIFF) return b;
        res.mantissa = b.mantissa + a.mantissa * invPowersOf2[diff];
        res.exponent = b.exponent;
    }
    exex_normalize(&res);
    return res;
}
ExExFloat exex_subtract(ExExFloat a, ExExFloat b) { b.mantissa *= -1.0; return exex_add(a, b); }
ExExFloat exex_multiply(ExExFloat a, ExExFloat b) {
    ExExFloat res;
    res.mantissa = a.mantissa * b.mantissa;
    res.exponent = a.exponent + b.exponent;
    exex_normalize(&res);
    return res;
}
ExExFloat exex_divide(ExExFloat a, ExExFloat b) {
    ExExFloat res;
    if (b.mantissa == 0.0) {
        res.mantissa = (a.mantissa > 0 ? 1.0 : -1.0) * INFINITY;
        res.exponent = 0;
        return res;
    }
    res.mantissa = a.mantissa / b.mantissa;
    res.exponent = a.exponent - b.exponent;
    exex_normalize(&res);
    return res;
}

// --- Global Init/Cleanup (Correct) ---
void divdiff_global_init() {
    if (invPowersOf2) return;
    invPowersOf2 = (double*)malloc(MAX_EXP_DIFF * sizeof(double));
    invPowersOf2[0] = 1.0;
    for (int i = 1; i < MAX_EXP_DIFF; ++i) invPowersOf2[i] = invPowersOf2[i - 1] * 0.5;
}
void divdiff_global_cleanup() {
    free(invPowersOf2);
    invPowersOf2 = NULL;
}

// --- DivDiff Object Management (Correct) ---
DivDiff* divdiff_create(int max_q, int s_max) {
    DivDiff* dd = (DivDiff*)malloc(sizeof(DivDiff));
    dd->max_len = max_q + 1;
    dd->s_max = s_max;
    dd->current_len = 0;
    dd->s = 1;
    dd->mu = 0.0;
    dd->energies = (double*)malloc(dd->max_len * sizeof(double));
    dd->h = (ExExFloat*)calloc(dd->max_len + EXTRA_LEN, sizeof(ExExFloat));
    dd->ddd = (ExExFloat*)malloc(dd->max_len * s_max * sizeof(ExExFloat));
    dd->results = (ExExFloat*)malloc(dd->max_len * sizeof(ExExFloat));
    dd->exp_mu = exex_from_double(1.0);
    return dd;
}
void divdiff_free(DivDiff* dd) {
    if (dd) {
        free(dd->energies);
        free(dd->h);
        free(dd->ddd);
        free(dd->results);
        free(dd);
    }
}

// --- Core O(q) Calculation Functions ---

// This is the internal worker function, a direct translation of the C++ AddElement
static void divdiff_add_element_internal(DivDiff* dd, double new_energy, int force_s, double force_mu) {
    int n = dd->current_len;
    int N_h_buffer = dd->max_len + EXTRA_LEN;
    ExExFloat curr;

    dd->energies[n] = new_energy;
    dd->current_len++;

    if (dd->current_len == 1) {
        // Base Case: Initialize state for the first element
        dd->s = (force_s == 0) ? 1 : force_s;
        dd->mu = (force_mu == 0) ? dd->energies[0] : force_mu;
        exex_init_expmu(&dd->exp_mu, dd->mu);

        dd->h[0] = exex_from_double(1.0);
        for (int k = 1; k < N_h_buffer; k++) {
            dd->h[k] = exex_divide(dd->h[k - 1], exex_from_double(dd->s));
        }

        if (dd->mu != dd->energies[0]) {
            for (int k = N_h_buffer - 1; k > 0; --k) {
                dd->h[k-1] = exex_add(dd->h[k-1], exex_multiply_double(dd->h[k], (dd->energies[0] - dd->mu) / k));
            }
        }
        curr = exex_multiply(dd->exp_mu, dd->h[0]);
        for (int k = 0; k < dd->s - 1; k++) {
            dd->ddd[k * dd->max_len] = curr;
            curr = exex_multiply(curr, dd->h[0]);
        }
        // The result for a single point is 0! * exp(E0)
        dd->results[0] = exex_from_double(exp(dd->energies[0]));

    } else if (s_changed(dd) || dd->current_len >= dd->max_len) {
        // Fallback Case: State is unstable, recalculate everything
        divdiff_recalculate_all(dd, force_s, force_mu);
    } else {
        // Fast Path: O(q) incremental update
        for (int k = N_h_buffer - 1; k > n; --k) {
            dd->h[k - 1] = exex_add(dd->h[k - 1], exex_multiply_double(dd->h[k], (dd->energies[n] - dd->mu) / k));
        }
        curr = exex_multiply(dd->exp_mu, dd->h[n]);

        // This in-place update is a direct translation of the C++ logic and is correct.
        for (int k = n; k >= 1; --k) {
            ExExFloat term1 = exex_multiply_double(dd->h[k - 1], n);
            ExExFloat term2 = exex_multiply_double(dd->h[k], dd->energies[n] - dd->energies[n - k]);
            dd->h[k - 1] = exex_divide(exex_add(term1, term2), exex_from_double(n - k + 1));
        }

        for (int k = 0; k < dd->s - 1; ++k) {
            dd->ddd[k * dd->max_len + n] = curr;
            curr = exex_multiply(dd->ddd[k * dd->max_len], dd->h[n]);
            for (int j = 1; j <= n; ++j) {
                curr = exex_add(curr, exex_multiply(dd->ddd[k * dd->max_len + j], dd->h[n - j]));
            }
        }
        dd->results[n] = curr;
    }
}

// Public-facing add function
void divdiff_add_element(DivDiff* dd, double new_energy) {
    divdiff_add_element_internal(dd, new_energy, 0, 0.0);
}

void divdiff_remove_last_element(DivDiff* dd) {
    if (dd->current_len > 0) {
        int n = dd->current_len - 1;
        int N_h_buffer = dd->max_len + EXTRA_LEN;

        for (int k = 1; k <= n; ++k) {
            ExExFloat term1 = exex_multiply_double(dd->h[k - 1], n - k + 1);
            ExExFloat term2 = exex_multiply_double(dd->h[k], dd->energies[n] - dd->energies[n - k]);
            dd->h[k - 1] = exex_divide(exex_subtract(term1, term2), exex_from_double(n));
        }
        for (int k = n + 1; k < N_h_buffer; ++k) {
            dd->h[k - 1] = exex_subtract(dd->h[k - 1], exex_multiply_double(dd->h[k], (dd->energies[n] - dd->mu) / k));
        }
        dd->current_len--;
    }
}

static void divdiff_recalculate_all(DivDiff* dd, int force_s, double force_mu) {
    int len = dd->current_len;
    if (len == 0) return;

    double* temp_energies = (double*)malloc(len * sizeof(double));
    memcpy(temp_energies, dd->energies, len * sizeof(double));

    dd->current_len = 0;
    
    int s_new = force_s;
    if (s_new == 0) {
        s_new = (int)ceil(max_abs_diff(temp_energies, len) / 3.5);
        if (s_new == 0) s_new = 1;
    }
    
    double mu_new = (force_mu == 0) ? calculate_mean(temp_energies, len) : force_mu;

    for (int i = 0; i < len; ++i) {
        divdiff_add_element_internal(dd, temp_energies[i], s_new, mu_new);
    }

    free(temp_energies);
}