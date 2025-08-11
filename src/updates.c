// File: src/updates.c
// Purpose: Implements the Monte Carlo update moves.
// VERSION: Corrected weight calculation and detailed balance factors.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "updates.h"
#include "utils.h"

// --- Private Helper: State Management ---

static void save_state(const QMCState* state, QMCState* backup) {
    backup->q = state->q;
    bitset_copy(backup->lattice, state->lattice);
    memcpy(backup->Sq, state->Sq, state->q * sizeof(int));
}

static void restore_state(QMCState* state, const QMCState* backup) {
    state->q = backup->q;
    bitset_copy(state->lattice, backup->lattice);
    memcpy(state->Sq, backup->Sq, backup->q * sizeof(int));
}

// Helper to get the full weight of a state.
// The weight is W = Re(D(z,Sq)) * ((-beta)^q / q!) * (q! * exp[...])
static ExExFloat get_full_weight(const QMCState* state) {
    if (state->q < 0 || state->q >= state->weight_calculator->current_len) return exex_from_double(0.0);
    
    ExExFloat divdiff_part = state->weight_calculator->results[state->q];
    ExExFloat beta_part = state->beta_pow_factorial[state->q];
    ExExFloat combined = exex_multiply(divdiff_part, beta_part);
    
    return exex_multiply_double(combined, creal(state->currD));
}

// Helper to rebuild the divdiff table from the state's Energies array
static void rebuild_divdiff(QMCState* state, const SimParams* params) {
    divdiff_free(state->weight_calculator);
    state->weight_calculator = divdiff_create(params->QMAX, 500);
    for (int i = 0; i <= state->q; ++i) {
        divdiff_add_element(state->weight_calculator, -params->BETA * state->Energies[i]);
    }
}

// --- Private Helper: Update Moves ---

static void update_classical_flip(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    QMCState* backup = state_create(params);
    save_state(state, backup);
    ExExFloat old_weight = get_full_weight(state);

    int spin_to_flip = rng_uniform_int(params->N);
    bitset_flip(state->lattice, spin_to_flip);

    state_recalculate_props(state, h);
    rebuild_divdiff(state, params);
    ExExFloat new_weight = get_full_weight(state);

    double acceptance_prob = 1.0;
    if (exex_to_double(old_weight) != 0.0) {
        double ratio = fabs(exex_to_double(exex_divide(new_weight, old_weight)));
        if (ratio < 1.0) acceptance_prob = ratio;
    }

    if (rng_uniform_double() >= acceptance_prob) { // REJECT
        restore_state(state, backup);
        state_recalculate_props(state, h);
        rebuild_divdiff(state, params);
    }
    state_free(backup);
}

static void update_pair_insertion_deletion(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    QMCState* backup = state_create(params);
    save_state(state, backup);
    ExExFloat old_weight = get_full_weight(state);

    double metropolis_factor = 1.0;
    int choice = rng_uniform_int(2);

    if (choice == 0 && state->q >= 2) { // --- ATTEMPT DELETION ---
        int m = rng_uniform_int(state->q - 1);
        if (state->Sq[m] == state->Sq[m+1]) {
            int old_q = state->q;
            state->q -= 2;
            memmove(&state->Sq[m], &state->Sq[m+2], (old_q - (m + 2)) * sizeof(int));
            metropolis_factor = (double)params->NOP;
        } else {
            state_free(backup); return;
        }
    } else if (choice == 1 && state->q + 2 < params->QMAX) { // --- ATTEMPT INSERTION ---
        int m = rng_uniform_int(state->q + 1);
        int op_to_insert = rng_uniform_int(params->NOP);
        int old_q = state->q;
        state->q += 2;
        memmove(&state->Sq[m+2], &state->Sq[m], (old_q - m) * sizeof(int));
        state->Sq[m] = op_to_insert;
        state->Sq[m+1] = op_to_insert;
        metropolis_factor = 1.0 / (double)params->NOP;
    } else {
        state_free(backup); return;
    }

    state_recalculate_props(state, h);
    rebuild_divdiff(state, params);
    ExExFloat new_weight = get_full_weight(state);

    double acceptance_prob = 1.0;
    if (exex_to_double(old_weight) != 0.0) {
        double ratio = fabs(exex_to_double(exex_divide(new_weight, old_weight)));
        if (ratio * metropolis_factor < 1.0) {
            acceptance_prob = ratio * metropolis_factor;
        }
    } else if (exex_to_double(new_weight) == 0.0) {
        acceptance_prob = 0.0;
    }

    if (rng_uniform_double() >= acceptance_prob) { // REJECT
        restore_state(state, backup);
        state_recalculate_props(state, h);
        rebuild_divdiff(state, params);
    }
    state_free(backup);
}

static void update_simple_swap(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    if (state->q < 2) return;

    int m = rng_uniform_int(state->q - 1);
    if (state->Sq[m] == state->Sq[m+1]) return;

    QMCState* backup = state_create(params);
    save_state(state, backup);
    ExExFloat old_weight = get_full_weight(state);

    int temp = state->Sq[m];
    state->Sq[m] = state->Sq[m+1];
    state->Sq[m+1] = temp;

    state_recalculate_props(state, h);
    rebuild_divdiff(state, params);
    ExExFloat new_weight = get_full_weight(state);

    double acceptance_prob = 1.0;
    if (exex_to_double(old_weight) != 0.0) {
        double ratio = fabs(exex_to_double(exex_divide(new_weight, old_weight)));
        if (ratio < 1.0) acceptance_prob = ratio;
    }

    if (rng_uniform_double() >= acceptance_prob) { // REJECT
        restore_state(state, backup);
        state_recalculate_props(state, h);
        rebuild_divdiff(state, params);
    }
    state_free(backup);
}


void do_update(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    double rand_val = rng_uniform_double();

    if (rand_val < 0.333) {
        update_classical_flip(state, h, params);
    } else if (rand_val < 0.666) {
        update_pair_insertion_deletion(state, h, params);
    } else {
        update_simple_swap(state, h, params);
    }
}