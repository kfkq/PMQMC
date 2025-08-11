// File: src/state.c
// Purpose: Implements the logic for managing the QMC simulation state.
// VERSION: Pre-computes beta_pow_factorial array.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "state.h"

// (state_calculate_d_k function is unchanged)
static complex_t state_calculate_d_k(const Hamiltonian* h, const bitset_t* config, int k) {
    complex_t sum = 0.0;
    for (int i = 0; i < h->D_sizes[k]; ++i) {
        int overlap = 0;
        for(int bit = 0; bit < config->num_bits; ++bit) {
            if (bitset_get(config, bit) == 1 && bitset_get(h->D_products[k][i], bit) == 1) {
                overlap++;
            }
        }
        if (overlap % 2 != 0) {
            sum -= h->D_coeffs[k][i];
        } else {
            sum += h->D_coeffs[k][i];
        }
    }
    return sum;
}

QMCState* state_create(const SimParams* params) {
    QMCState* state = (QMCState*)calloc(1, sizeof(QMCState));
    if (!state) return NULL;

    state->lattice = bitset_create(params->N);
    state->Sq = malloc(params->QMAX * sizeof(int));
    state->Energies = malloc((params->QMAX + 1) * sizeof(double));
    state->weight_calculator = divdiff_create(params->QMAX, 500);
    state->beta_pow_factorial = malloc(params->QMAX * sizeof(ExExFloat));

    if (!state->lattice || !state->Sq || !state->Energies || !state->weight_calculator || !state->beta_pow_factorial) {
        fprintf(stderr, "Error: Failed to allocate memory for QMC state internals.\n");
        state_free(state);
        return NULL;
    }

    for (int i = 0; i < params->N; ++i) {
        if (rand() % 2) {
            bitset_set(state->lattice, i);
        }
    }

    state->q = 0;
    state->currD = 1.0 + 0.0 * I;

    // Pre-compute the (-beta)^k / k! terms
    state->beta_pow_factorial[0] = exex_from_double(1.0);
    for (int k = 1; k < params->QMAX; ++k) {
        ExExFloat temp = exex_multiply_double(state->beta_pow_factorial[k-1], -params->BETA / k);
        state->beta_pow_factorial[k] = temp;
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
        free(state);
    }
}

// (state_calculate_classical_energy is unchanged)
double state_calculate_classical_energy(const Hamiltonian* h, const bitset_t* config) {
    complex_t total_energy = 0.0;
    for (int i = 0; i < h->D0_size; ++i) {
        int overlap = 0;
        for(int bit = 0; bit < config->num_bits; ++bit) {
            if (bitset_get(config, bit) == 1 && bitset_get(h->D0_product[i], bit) == 1) {
                overlap++;
            }
        }
        if (overlap % 2 != 0) {
            total_energy -= h->D0_coeff[i];
        } else {
            total_energy += h->D0_coeff[i];
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