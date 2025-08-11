// File: src/state.c
// Purpose: Implements the logic for managing the QMC simulation state.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "state.h"
#include "hamiltonian.h"

QMCState* state_create(const SimParams* params) {
    QMCState* state = (QMCState*)calloc(1, sizeof(QMCState));
    if (!state) {
        fprintf(stderr, "Error: Failed to allocate memory for QMC state.\n");
        return NULL;
    }

    state->lattice = bitset_create(params->N);
    state->Sq = malloc(params->QMAX * sizeof(int));
    state->Energies = malloc((params->QMAX + 1) * sizeof(double));
    state->weight_calculator = divdiff_create(params->QMAX, 500);
    state->factorials = malloc((params->QMAX + 1) * sizeof(double));

    if (!state->lattice || !state->Sq || !state->Energies || !state->weight_calculator || !state->factorials) {
        fprintf(stderr, "Error: Failed to allocate memory for QMC state internals.\n");
        state_free(state); // Use the new free function for cleanup
        return NULL;
    }

    for (int i = 0; i < params->N; ++i) {
        if (rand() % 2) {
            bitset_set(state->lattice, i);
        }
    }

    state->q = 0;
    // Note: Can't call state_update_energy_list here as Hamiltonian is not passed.
    // The caller of state_create will need to do the first energy calculation.
    // For simplicity, we'll assume the caller handles the first weight calculation.

    state->factorials[0] = 1.0;
    for (int i = 1; i <= params->QMAX; ++i) {
        state->factorials[i] = state->factorials[i - 1] * i;
    }

    return state;
}

void state_free(QMCState* state) {
    if (state) {
        bitset_free(state->lattice);
        free(state->Sq);
        free(state->Energies);
        divdiff_free(state->weight_calculator);
        free(state->factorials);
        free(state);
    }
}

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

void state_update_energy_list(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    bitset_t* temp_lattice = bitset_create(params->N);
    bitset_copy(temp_lattice, state->lattice);
    state->Energies[0] = state_calculate_classical_energy(h, temp_lattice);
    for (int i = 0; i < state->q; ++i) {
        int op_index = state->Sq[i];
        bitset_xor(temp_lattice, h->P_matrix[op_index]);
        state->Energies[i + 1] = state_calculate_classical_energy(h, temp_lattice);
    }
    bitset_free(temp_lattice);
}