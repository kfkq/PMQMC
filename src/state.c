// File: src/state.c
// Purpose: Implements the logic for managing the QMC simulation state.
// VERSION: Corrected sign logic in energy/d_k calculations to match C++ reference.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "state.h"
#include "utils.h"

// --- Private Helper Functions (static) ---

// Calculates d_k = <z | D_k | z> for a given configuration.
static complex_t state_calculate_d_k(const Hamiltonian* h, const bitset_t* config, int k) {
    complex_t sum = 0.0;
    for (int i = 0; i < h->D_sizes[k]; ++i) {
        // The C++ reference calculates parity based on sites where the Z-string is 1
        // and the lattice is 0. Let's replicate that logic exactly.
        int parity = 0;
        for(int bit = 0; bit < config->num_bits; ++bit) {
            if (bitset_get(h->D_products[k][i], bit) == 1 && bitset_get(config, bit) == 0) {
                parity++;
            }
        }

        // C++ logic: sum -= (2 * (parity % 2) - 1) * coeff
        if (parity % 2 == 0) { // Even parity
            sum += h->D_coeffs[k][i];
        } else { // Odd parity
            sum -= h->D_coeffs[k][i];
        }
    }
    return sum;
}

// --- Public API Functions ---

QMCState* state_create(const SimParams* params) {
    QMCState* state = (QMCState*)calloc(1, sizeof(QMCState));
    if (!state) return NULL;

    state->lattice = bitset_create(params->N);
    state->Sq = malloc(params->QMAX * sizeof(int));
    state->Energies = malloc((params->QMAX + 1) * sizeof(double));
    state->weight_calculator = divdiff_create(params->QMAX, 500);
    state->beta_pow_factorial = malloc(params->QMAX * sizeof(ExExFloat));
    state->worm = calloc(1, sizeof(Worm));

    if (!state->lattice || !state->Sq || !state->Energies || !state->weight_calculator || !state->beta_pow_factorial || !state->worm) {
        fprintf(stderr, "Error: Failed to allocate memory for QMC state internals.\n");
        state_free(state);
        return NULL;
    }

    state->worm->z_k = bitset_create(params->N);
    state->worm->z_l = bitset_create(params->N);
    if (!state->worm->z_k || !state->worm->z_l) {
        fprintf(stderr, "Error: Failed to allocate memory for worm bitsets.\n");
        state_free(state);
        return NULL;
    }

    for (int i = 0; i < params->N; ++i) {
        if (rng_uniform_int(2)) {
            bitset_set(state->lattice, i);
        }
    }

    state->q = 0;
    state->currD = 1.0 + 0.0 * I;

    state->beta_pow_factorial[0] = exex_from_double(1.0);
    for (int k = 1; k < params->QMAX; ++k) {
        state->beta_pow_factorial[k] = exex_multiply_double(state->beta_pow_factorial[k-1], -params->BETA / k);
    }

    return state;
}

void state_free(QMCState* state) {
    if (state) {
        bitset_free(state->lattice);
        free(state->Sq);
        free(state->Energies);
        divdiff_free(state->weight_calculator);
        free(state->beta_pow_factorial);
        if (state->worm) {
            bitset_free(state->worm->z_k);
            bitset_free(state->worm->z_l);
            free(state->worm);
        }
        free(state);
    }
}

double state_calculate_classical_energy(const Hamiltonian* h, const bitset_t* config) {
    complex_t total_energy = 0.0;
    for (int i = 0; i < h->D0_size; ++i) {
        // Replicating C++ logic: count sites where Z-string is 1 and lattice is 0.
        int parity = 0;
        for(int bit = 0; bit < config->num_bits; ++bit) {
            if (bitset_get(h->D0_product[i], bit) == 1 && bitset_get(config, bit) == 0) {
                parity++;
            }
        }
        
        // C++ logic: sum -= (2 * (parity % 2) - 1) * coeff
        if (parity % 2 == 0) { // Even parity
            total_energy += h->D0_coeff[i];
        } else { // Odd parity
            total_energy -= h->D0_coeff[i];
        }
    }
    return creal(total_energy);
}

void state_recalculate_props(QMCState* state, const Hamiltonian* h) {
    bitset_t* temp_lattice = bitset_create(h->P_matrix[0]->num_bits);
    bitset_copy(temp_lattice, state->lattice);

    state->Energies[0] = state_calculate_classical_energy(h, temp_lattice);
    state->currD = 1.0 + 0.0 * I;

    for (int i = 0; i < state->q; ++i) {
        int op_index = state->Sq[i];
        state->currD *= state_calculate_d_k(h, temp_lattice, op_index);
        
        bitset_xor(temp_lattice, h->P_matrix[op_index]);
        state->Energies[i + 1] = state_calculate_classical_energy(h, temp_lattice);
    }
    bitset_free(temp_lattice);
}

void rebuild_divdiff_from_energies(QMCState* state, const SimParams* params) {
    divdiff_free(state->weight_calculator);
    state->weight_calculator = divdiff_create(params->QMAX, 500);
    for (int i = 0; i <= state->q; ++i) {
        divdiff_add_element(state->weight_calculator, -params->BETA * state->Energies[i]);
    }
}

ExExFloat get_full_weight(const QMCState* state) {
    if (state->q < 0 || state->q >= state->weight_calculator->current_len) {
        return exex_from_double(0.0);
    }
    ExExFloat divdiff_part = state->weight_calculator->results[state->q];
    ExExFloat beta_part = state->beta_pow_factorial[state->q];
    ExExFloat combined = exex_multiply(divdiff_part, beta_part);
    return exex_multiply_double(combined, creal(state->currD));
}