// File: src/state.h
// Purpose: Declares the dynamic state of the QMC simulation.
// VERSION: Added beta_pow_factorial for correct weight calculation.

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
    double* Energies;
    complex_t currD;
    DivDiff* weight_calculator;
    ExExFloat* beta_pow_factorial; // CRITICAL: Stores (-beta)^k / k!
} QMCState;

// --- Public API Functions ---
QMCState* state_create(const SimParams* params);
void state_free(QMCState* state);
double state_calculate_classical_energy(const Hamiltonian* h, const bitset_t* config);
void state_recalculate_props(QMCState* state, const Hamiltonian* h);

#endif // STATE_H