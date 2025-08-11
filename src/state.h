// File: src/state.h
// Purpose: Declares the dynamic state of the QMC simulation and the functions
//          to initialize and manage it.

#ifndef STATE_H
#define STATE_H

#include "datatypes.h"
#include "divdiff.h"
#include "hamiltonian.h"

// --- QMC State Structure ---
// This struct encapsulates the entire dynamic state of the simulation.
typedef struct {
    bitset_t* lattice;
    int q;
    int* Sq;
    double* Energies;
    DivDiff* weight_calculator;
    double* factorials;
} QMCState;

// --- Public API Functions ---

/**
 * @brief Creates and initializes the QMC simulation state.
 *        Allocates memory for all state variables and sets up a random initial configuration.
 * @param params A pointer to the loaded simulation parameters.
 * @return A pointer to the newly created QMCState, or NULL on failure.
 */
QMCState* state_create(const SimParams* params);

/**
 * @brief Frees all memory associated with the global simulation state.
 * @param state A pointer to the QMCState to be freed.
 */
void state_free(QMCState* state);

/**
 * @brief Calculates the energy of a given classical spin configuration.
 * @param h The Hamiltonian data.
 * @param config The bitset representing the spin configuration.
 * @return The classical energy <config|D0|config>.
 */
double state_calculate_classical_energy(const Hamiltonian* h, const bitset_t* config);

/**
 * @brief Populates the global 'Energies' array based on the current 'lattice' and 'Sq'.
 *        This is a key helper function called by update moves.
 * @param state The current QMC state to update.
 * @param h The Hamiltonian data.
 * @param params The simulation parameters.
 */
void state_update_energy_list(QMCState* state, const Hamiltonian* h, const SimParams* params);

#endif // STATE_H