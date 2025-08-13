// File: src/measurements.h
// Purpose: Declares functions to measure physical observables from the QMC state.

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "state.h"
#include "datatypes.h"

/**
 * @brief Measures the instantaneous value of the energy <H>.
 * @param state The current QMC state.
 * @param params The simulation parameters.
 * @return The instantaneous energy value.
 */
double measure_H(const QMCState* state, const SimParams* params);

/**
 * @brief Measures the instantaneous value of the energy squared <H^2>.
 * @param state The current QMC state.
 * @param params The simulation parameters.
 * @return The instantaneous energy squared value.
 */
double measure_H2(const QMCState* state, const SimParams* params);

#endif // MEASUREMENTS_H