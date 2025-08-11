// File: src/state.h
// Purpose: Declares the dynamic state of the QMC simulation and the functions
//          to initialize and manage it.

#ifndef STATE_H
#define STATE_H

#include "datatypes.h"
#include "divdiff.h"

// --- Global Simulation State Variables ---
// These variables represent the current configuration of the Markov chain.
// They are defined in state.c and made accessible via 'extern'.

extern bitset_t* lattice;       // The current classical spin configuration
extern int q;                   // The current length of the operator sequence
extern int* Sq;                  // The sequence of operator indices {P_i, P_j, ...}
extern double* Energies;        // The list of energies {E_z0, E_z1, ...} for the current sequence
extern DivDiff* weight_calculator; // The calculator for the configuration weight
extern double* factorials;

// --- Public API Functions ---

/**
 * @brief Initializes the global simulation state.
 *        Allocates memory for all state variables and sets up a random initial configuration.
 * @param params A pointer to the loaded simulation parameters.
 * @return 0 on success, -1 on failure.
 */
int state_init(const SimParams* params);

/**
 * @brief Frees all memory associated with the global simulation state.
 */
void state_free();

/**
 * @brief Calculates the energy of a given classical spin configuration.
 * @param config The bitset representing the spin configuration.
 * @return The classical energy <config|D0|config>.
 */
double state_calculate_classical_energy(const bitset_t* config);

/**
 * @brief Populates the global 'Energies' array based on the current 'lattice' and 'Sq'.
 *        This is a key helper function called by update moves.
 */
void state_update_energy_list();

#endif // STATE_H