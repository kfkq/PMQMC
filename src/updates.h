// File: src/updates.h
// Purpose: Declares the functions that perform the Monte Carlo updates.

#ifndef UPDATES_H
#define UPDATES_H

#include "state.h"
#include "hamiltonian.h"

/**
 * @brief Performs one Monte Carlo update step.
 *        This function randomly chooses one of the available update moves
 *        and attempts to perform it.
 * @param state The current QMC state (will be modified).
 * @param h The Hamiltonian data.
 * @param params The simulation parameters.
 */
void do_update(QMCState* state, const Hamiltonian* h, const SimParams* params);

#endif // UPDATES_H