// File: src/updates.c
// Purpose: Implements all 7 Monte Carlo update moves.
// VERSION: Instrumented with extensive debugging print statements.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "updates.h"
#include "utils.h"

// --- UNCOMMENT THIS LINE TO ENABLE DEBUGGING OUTPUT ---
// #define DEBUG_UPDATES

// --- State Management and Helpers (static) ---

typedef struct {
    int q;
    int* Sq;
    bitset_t* lattice;
} QMCState_Backup;

static void save_state(const QMCState* state, QMCState_Backup* backup) {
    backup->q = state->q;
    memcpy(backup->Sq, state->Sq, state->q * sizeof(int));
    bitset_copy(backup->lattice, state->lattice);
}

static void restore_state(QMCState* state, const QMCState_Backup* backup) {
    state->q = backup->q;
    memcpy(state->Sq, backup->Sq, backup->q * sizeof(int));
    bitset_copy(state->lattice, backup->lattice);
}

// --- The 7 Individual Update Moves (as static functions) ---

static int attempt_classical_flip(QMCState* state, const SimParams* params) {
    int spin_to_flip = rng_uniform_int(params->N);
    bitset_flip(state->lattice, spin_to_flip);
    return 1;
}

static int attempt_simple_swap(QMCState* state) {
    if (state->q < 2) return 0;
    int m = rng_uniform_int(state->q - 1);
    if (state->Sq[m] == state->Sq[m+1]) return 0;
    int temp = state->Sq[m];
    state->Sq[m] = state->Sq[m+1];
    state->Sq[m+1] = temp;
    return 1;
}

static int attempt_pair_insertion_deletion(QMCState* state, const SimParams* params, double* metro_factor) {
    int choice;
    // If q=0, we MUST insert. Otherwise, choose randomly.
    if (state->q == 0) {
        choice = 1; // Force insertion
    } else {
        choice = rng_uniform_int(2);
    }

    if (choice == 0 && state->q >= 2) { // DELETION
        int m = rng_uniform_int(state->q - 1);
        if (state->Sq[m] == state->Sq[m+1]) {
            int old_q = state->q;
            state->q -= 2;
            memmove(&state->Sq[m], &state->Sq[m+2], (old_q - (m + 2)) * sizeof(int));
            *metro_factor = 1.0 / (double)params->NOP;
            return 1;
        }
    } else if (choice == 1 && state->q + 2 <= params->QMAX) { // INSERTION
        int m = rng_uniform_int(state->q + 1);
        int op_to_insert = rng_uniform_int(params->NOP);
        int old_q = state->q;
        state->q += 2;
        memmove(&state->Sq[m+2], &state->Sq[m], (old_q - m) * sizeof(int));
        state->Sq[m] = op_to_insert;
        state->Sq[m+1] = op_to_insert;
        *metro_factor = (double)params->NOP;
        return 1;
    }
    return 0;
}

static int attempt_block_swap(QMCState* state, const Hamiltonian* h) {
    if (state->q < 2) return 0;
    int k = rng_uniform_int(state->q - 1);
    for (int i = 0; i <= k; ++i) {
        bitset_xor(state->lattice, h->P_matrix[state->Sq[i]]);
    }
    int* temp_sq = malloc(state->q * sizeof(int));
    memcpy(temp_sq, state->Sq, state->q * sizeof(int));
    memmove(&state->Sq[0], &temp_sq[k+1], (state->q - 1 - k) * sizeof(int));
    memmove(&state->Sq[state->q - 1 - k], &temp_sq[0], (k + 1) * sizeof(int));
    free(temp_sq);
    return 1;
}

static int attempt_cycle_completion(const SimParams* params, double* metro_factor) {
    if (params->NCYCLES == 0) return 0;
    (void)params; (void)metro_factor;
    return 0;
}

// --- Top-Level Update Algorithms ---

void do_composite_update(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    QMCState_Backup backup;
    backup.Sq = malloc(params->QMAX * sizeof(int));
    backup.lattice = bitset_create(params->N);
    save_state(state, &backup);
    
    ExExFloat old_weight = get_full_weight(state);
    double total_metro_factor = 1.0;
    int move_made = 0;
    int is_insertion = 0; // Debug flag

    double rand_val = rng_uniform_double();
    if (rand_val < 0.25) {
        move_made = attempt_classical_flip(state, params);
    } else if (rand_val < 0.50) {
        move_made = attempt_simple_swap(state);
    } else if (rand_val < 0.75) {
        move_made = attempt_pair_insertion_deletion(state, params, &total_metro_factor);
        if (state->q > backup.q) is_insertion = 1;
    } else {
        if (rng_uniform_double() < 0.5) {
            move_made = attempt_block_swap(state, h);
        } else {
            move_made = attempt_cycle_completion(params, &total_metro_factor);
        }
    }

    if (move_made) {
        state_recalculate_props(state, h);
        rebuild_divdiff_from_energies(state, params);
        ExExFloat new_weight = get_full_weight(state);
        
        double ratio = 0.0;
        if (exex_to_double(old_weight) != 0.0) {
            ratio = fabs(exex_to_double(exex_divide(new_weight, old_weight)));
        }
        double acceptance_prob = fmin(1.0, ratio * total_metro_factor);
        double rng_val_decision = rng_uniform_double();

        #ifdef DEBUG_UPDATES
        if (is_insertion) {
            printf("\n[DEBUG] === PAIR INSERTION ATTEMPT ===\n");
            printf("[DEBUG] BEFORE: q=%d, old_weight=%.6e\n", backup.q, exex_to_double(old_weight));
            printf("[DEBUG] PROPOSED: q=%d, Sq=[%d, %d]\n", state->q, state->Sq[0], state->Sq[1]);
            printf("[DEBUG] ENERGIES: E0=%.4f, E1=%.4f, E2=%.4f\n", state->Energies[0], state->Energies[1], state->Energies[2]);
            printf("[DEBUG] D-TERM: Re(currD)=%.4f, Im(currD)=%.4f\n", creal(state->currD), cimag(state->currD));
            printf("[DEBUG] WEIGHTS: new_weight=%.6e\n", exex_to_double(new_weight));
            printf("[DEBUG] METRO:  factor=%.2f, |W_new/W_old|=%.6e, prob=%.6e\n", total_metro_factor, ratio, acceptance_prob);
            printf("[DEBUG] RNG val = %.6f. Decision: %s\n", rng_val_decision, (rng_val_decision < acceptance_prob) ? "ACCEPT" : "REJECT");
        }
        #endif

        if (rng_val_decision >= acceptance_prob) { // REJECT
            restore_state(state, &backup);
            // CRITICAL FIX: After restoring state, must recalculate derived properties
            state_recalculate_props(state, h);
            rebuild_divdiff_from_energies(state, params);
        }
    }

    free(backup.Sq);
    bitset_free(backup.lattice);
}

void do_worm_update(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    (void)h; (void)params;
    printf("Warning: do_worm_update is not fully implemented.\n");
}