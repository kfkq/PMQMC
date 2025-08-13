// File: src/measurements.c
// Purpose: Implements the measurement of physical observables.
// Formulas are derived from the C++ reference and the associated paper.

#include <math.h>
#include "measurements.h"
#include "datatypes.h" // For ExExFloat etc.

// Note: The 'params' argument is unused in these functions but is kept
// for consistency with a potential future interface where it might be needed.
// This avoids breaking function signatures if more parameters are required later.

double measure_H(const QMCState* state, const SimParams* params) {
    (void)params; // Suppress unused parameter warning

    if (state->q < 0 || !state->weight_calculator || state->q >= state->weight_calculator->current_len) {
        return 0.0;
    }

    // The first term is E_q, which is d->z[q] / (-beta) in the C++ code.
    // In our C code, Energies[q] already holds this value.
    double R = state->Energies[state->q];

    if (state->q > 0) {
        // The second term is (d[q-1]/d[q]) * q / (-beta)
        ExExFloat ratio = exex_divide(
            state->weight_calculator->results[state->q - 1],
            state->weight_calculator->results[state->q]
        );
        R += exex_to_double(ratio) * state->q / (-params->BETA);
    }

    return R;
}

double measure_H2(const QMCState* state, const SimParams* params) {
    (void)params; // Suppress unused parameter warning

    if (state->q < 0 || !state->weight_calculator || state->q >= state->weight_calculator->current_len) {
        return 0.0;
    }

    const double neg_beta = -params->BETA;

    // Term 1: E_q^2
    double R = state->Energies[state->q] * state->Energies[state->q];

    if (state->q > 0) {
        // Term 2: (E_q + E_{q-1}) * (d[q-1]/d[q]) * q / (-beta)
        ExExFloat ratio1 = exex_divide(
            state->weight_calculator->results[state->q - 1],
            state->weight_calculator->results[state->q]
        );
        R += (state->Energies[state->q] + state->Energies[state->q - 1]) * exex_to_double(ratio1) * state->q / neg_beta;
    }

    if (state->q > 1) {
        // Term 3: (d[q-2]/d[q]) * q*(q-1) / beta^2
        ExExFloat ratio2 = exex_divide(
            state->weight_calculator->results[state->q - 2],
            state->weight_calculator->results[state->q]
        );
        R += exex_to_double(ratio2) * state->q * (state->q - 1) / (neg_beta * neg_beta);
    }

    return R;
}

double measure_Z_magnetization(const QMCState* state, const SimParams* params) {
    // Magnetization is defined as (N_up - N_down) / N_total
    // N_up = number of 1s, N_down = number of 0s
    // N_up - N_down = count - (N - count) = 2*count - N
    double count = (double)bitset_count(state->lattice);
    return (2.0 * count - (double)params->N) / (double)params->N;
}
