// File: src/updates.h
// Purpose: Declares the functions that perform the Monte Carlo updates.
// VERSION: Extended to include all 7 update types via two main algorithms.

#ifndef UPDATES_H
#define UPDATES_H

#include "state.h"
#include "hamiltonian.h"

/**
 * @brief Performs one composite Monte Carlo update step.
 *        This is the standard update algorithm. It combines local moves
 *        (swap, pair ins/del, classical, block swap) with the global
 *        cycle completion move to ensure ergodicity.
 * @param state The current QMC state (will be modified).
 * @param h The Hamiltonian data.
 * @param params The simulation parameters.
 */
void do_composite_update(QMCState* state, const Hamiltonian* h, const SimParams* params);

/**
 * @brief Performs one full worm update sequence.
 *        This is an alternative global update algorithm that is often more
 *        efficient for systems with topological constraints. It starts by
 *        creating a "worm" (breaking the identity constraint) and propagates
 *        it until it "heals".
 * @param state The current QMC state (will be modified).
 * @param h The Hamiltonian data.
 * @param params The simulation parameters.
 */
void do_worm_update(QMCState* state, const Hamiltonian* h, const SimParams* params);

#endif // UPDATES_H