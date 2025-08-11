// File: src/state.c
// Purpose: Implements the logic for managing the QMC simulation state.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "state.h"
#include "hamiltonian.h"

bitset_t* lattice = NULL;
int q = 0;
int* Sq = NULL;
double* Energies = NULL;
DivDiff* weight_calculator = NULL;
double* factorials = NULL;

int state_init(const SimParams* params) {
    // ... (allocation block is unchanged) ...
    lattice = bitset_create(params->N);
    Sq = malloc(params->QMAX * sizeof(int));
    Energies = malloc((params->QMAX + 1) * sizeof(double));
    weight_calculator = divdiff_create(params->QMAX, 500);
    factorials = malloc((params->QMAX + 1) * sizeof(double)); // Allocate one extra slot

    if (!lattice || !Sq || !Energies || !weight_calculator || !factorials) {
        fprintf(stderr, "Error: Failed to allocate memory for QMC state.\n");
        if (lattice) bitset_free(lattice);
        free(Sq); free(Energies); free(factorials);
        divdiff_free(weight_calculator);
        return -1;
    }

    for (int i = 0; i < params->N; ++i) {
        if (rand() % 2) {
            bitset_set(lattice, i);
        }
    }

    q = 0;
    state_update_energy_list();
    divdiff_add_element(weight_calculator, Energies[0]);

    // for all possible values of q from 0 to QMAX.
    factorials[0] = 1.0;
    for (int i = 1; i <= params->QMAX; ++i) {
        factorials[i] = factorials[i - 1] * i;
    }

    return 0;
}

// ... (state_free and other functions are unchanged) ...
void state_free() {
    bitset_free(lattice);
    free(Sq);
    free(Energies);
    divdiff_free(weight_calculator);
    free(factorials);
}

double state_calculate_classical_energy(const bitset_t* config) {
    complex_t total_energy = 0.0;
    for (int i = 0; i < D0_size; ++i) {
        int overlap = 0;
        for(int bit = 0; bit < config->num_bits; ++bit) {
            if (bitset_get(config, bit) == 1 && bitset_get(D0_product[i], bit) == 1) {
                overlap++;
            }
        }
        if (overlap % 2 != 0) {
            total_energy -= D0_coeff[i];
        } else {
            total_energy += D0_coeff[i];
        }
    }
    return creal(total_energy);
}

void state_update_energy_list() {
    bitset_t* temp_lattice = bitset_create(lattice->num_bits);
    bitset_copy(temp_lattice, lattice);
    Energies[0] = state_calculate_classical_energy(temp_lattice);
    for (int i = 0; i < q; ++i) {
        int op_index = Sq[i];
        bitset_xor(temp_lattice, P_matrix[op_index]);
        Energies[i + 1] = state_calculate_classical_energy(temp_lattice);
    }
    bitset_free(temp_lattice);
}