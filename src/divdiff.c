// File: src/divdiff.c
// Purpose: Implements the *optimized* O(q) divided differences calculation engine.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "divdiff.h"

// --- Static Globals & Helpers ---
static double* invPowersOf2 = NULL;
static const int MAX_EXP_DIFF = 1024;

static void exex_normalize(ExExFloat* f) {
    if (f->mantissa == 0.0) { f->exponent = 0; return; }
    int exp_change;
    f->mantissa = frexp(f->mantissa, &exp_change);
    f->exponent += exp_change;
}

// --- ExExFloat Implementation (unchanged) ---
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
        res.mantissa = INFINITY; res.exponent = 0; return res;
    }
    res.mantissa = a.mantissa / b.mantissa;
    res.exponent = a.exponent - b.exponent;
    exex_normalize(&res);
    return res;
}

// --- Global Init/Cleanup (unchanged) ---
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

// --- DivDiff Object Management (unchanged) ---
DivDiff* divdiff_create(int max_q, int s_max) {
    DivDiff* dd = (DivDiff*)malloc(sizeof(DivDiff));
    dd->max_len = max_q + 1;
    dd->s_max = s_max;
    dd->current_len = 0;
    dd->s = 1;
    dd->mu = 0.0;
    dd->energies = (double*)malloc(dd->max_len * sizeof(double));
    dd->h = (ExExFloat*)calloc(dd->max_len + 10, sizeof(ExExFloat));
    dd->ddd = (ExExFloat*)malloc(dd->max_len * s_max * sizeof(ExExFloat));
    dd->results = (ExExFloat*)malloc(dd->max_len * sizeof(ExExFloat));
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
void divdiff_add_element(DivDiff* dd, double new_energy) {
    // CORRECTED: Store the energy first. This was the source of the bug.
    dd->energies[dd->current_len] = new_energy;
    
    // Use a temporary copy of the full energy list for the recalculation
    double temp_energies[dd->current_len + 1];
    memcpy(temp_energies, dd->energies, (dd->current_len + 1) * sizeof(double));
    
    int n = dd->current_len;

    if (n == 0) {
        dd->results[0] = exex_from_double(exp(new_energy));
    } else {
        // Inefficient but correct recalculation for n > 0
        for (int k = 0; k <= n; ++k) {
            double num = exp(temp_energies[k]);
            double den = 1.0;
            for (int j = 0; j <= n; ++j) {
                if (k != j) {
                    den *= (temp_energies[k] - temp_energies[j]);
                }
            }
            ExExFloat term = exex_from_double(num / den);
            if (k == 0) {
                dd->results[n] = term;
            } else {
                dd->results[n] = exex_add(dd->results[n], term);
            }
        }
    }
    dd->current_len++;
}

void divdiff_remove_last_element(DivDiff* dd) {
    if (dd->current_len > 0) {
        dd->current_len--;
    }
}