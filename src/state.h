// File: src/state.h
// Purpose: Declares the dynamic state of the QMC simulation.
// (This file is correct and does not need changes)

#ifndef STATE_H
#define STATE_H

#include "datatypes.h"
#include "divdiff.h"
#include "hamiltonian.h"

// --- QMC State Structure ---
typedef struct {
    bitset_t* lattice;
    int q;
    int* Sq;
    Worm* worm;
    double* Energies;
    complex_t currD;
    DivDiff* weight_calculator;
    ExExFloat* beta_pow_factorial;
    double* factorials;
} QMCState;

// --- Public API Functions ---
QMCState* state_create(const SimParams* params);
void state_free(QMCState* state);
double state_calculate_classical_energy(const Hamiltonian* h, const bitset_t* config);
void state_recalculate_props(QMCState* state, const Hamiltonian* h);

/**
 * @brief Rebuilds the internal divided differences table from the state's Energies array.
 */
void rebuild_divdiff_from_energies(QMCState* state, const SimParams* params);

/**
 * @brief Calculates the full, properly scaled weight of the current configuration.
 */
ExExFloat get_full_weight(const QMCState* state);

#endif // STATE_H