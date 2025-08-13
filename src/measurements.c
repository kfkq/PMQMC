// File: src/measurements.c
// Purpose: Implements the measurement of physical observables.
// Formulas are derived from the C++ reference and the associated paper.

#include <math.h>
#include "measurements.h"

double measure_H(const QMCState* state, const SimParams* params) {
    if (state->q < 0 || state->q >= state->weight_calculator->current_len) {
        return 0.0; // Should not happen in a valid state
    }

    // The first term is E_q, which is d->z[q] / (-beta) in the C++ code
    double R = state->Energies[state->q];

    if (state->q > 0) {
        // The second term is (d[q-1]/d[q]) * q / (-beta)
        ExExFloat ratio = exex_divide(
            state->weight_calculator->results[state->q - 1],
            state->weight_calculator->results[state->q]
        );
        R += exex_to_double(ratio) * state->q;
    }

    return R;
}

double measure_H2(const QMCState* state, const SimParams* params) {
    if (state->q < 0 || state->q >= state->weight_calculator->current_len) {
        return 0.0;
    }

    // Term 1: E_q^2
    double R = state->Energies[state->q] * state->Energies[state->q];

    if (state->q > 0) {
        // Term 2: (E_q + E_{q-1}) * (d[q-1]/d[q]) * q / (-beta)
        ExExFloat ratio1 = exex_divide(
            state->weight_calculator->results[state->q - 1],
            state->weight_calculator->results[state->q]
        );
        R += (state->Energies[state->q] + state->Energies[state->q - 1]) * exex_to_double(ratio1) * state->q;
    }

    if (state->q > 1) {
        // Term 3: (d[q-2]/d[q]) * q*(q-1) / beta^2
        ExExFloat ratio2 = exex_divide(
            state->weight_calculator->results[state->q - 2],
            state->weight_calculator->results[state->q]
        );
        R += exex_to_double(ratio2) * state->q * (state->q - 1);
    }

    return R;
}