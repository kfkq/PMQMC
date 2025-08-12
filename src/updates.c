// File: src/updates.c
// Purpose: Implements all 7 Monte Carlo update moves.

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
    backup->Sq = malloc(state->q * sizeof(int));  // Alloc on save
    memcpy(backup->Sq, state->Sq, state->q * sizeof(int));
    backup->lattice = bitset_create(state->lattice->num_bits);
    bitset_copy(backup->lattice, state->lattice);
}

static void restore_state(QMCState* state, const QMCState_Backup* backup) {
    state->q = backup->q;
    memcpy(state->Sq, backup->Sq, backup->q * sizeof(int));
    bitset_copy(state->lattice, backup->lattice);
}

static void free_backup(QMCState_Backup* backup) {
    free(backup->Sq);
    bitset_free(backup->lattice);
}

// --- The 7 Individual Update Moves (as static functions) ---

static int attempt_classical_flip(QMCState* state, const SimParams* params) {
    int spin_to_flip = rng_uniform_int(params->N);
    bitset_flip(state->lattice, spin_to_flip);
#ifdef DEBUG_UPDATES
    printf("[DEBUG] Classical flip on spin %d\n", spin_to_flip);
#endif
    return 1;  // Always succeeds
}

static int attempt_simple_swap(QMCState* state) {
    if (state->q < 2) return 0;
    int m = rng_uniform_int(state->q - 1);
    if (state->Sq[m] == state->Sq[m+1]) return 0;  // No swap if identical
    int temp = state->Sq[m];
    state->Sq[m] = state->Sq[m+1];
    state->Sq[m+1] = temp;
#ifdef DEBUG_UPDATES
    printf("[DEBUG] Simple swap at position %d\n", m);
#endif
    return 1;
}

static int attempt_pair_insertion_deletion(QMCState* state, const SimParams* params, double* metro_factor) {
    int choice;
    if (state->q == 0) {
        choice = 1; // Force insertion
    } else {
        choice = rng_uniform_int(2);
    }

    if (choice == 0 && state->q >= 2) { // DELETION
        int m = rng_uniform_int(state->q - 1);
        if (state->Sq[m] != state->Sq[m+1]) return 0;
        int old_q = state->q;
        state->q -= 2;
        memmove(&state->Sq[m], &state->Sq[m+2], (old_q - (m + 2)) * sizeof(int));
        *metro_factor = 1.0 / (double)params->NOP;
#ifdef DEBUG_UPDATES
        printf("[DEBUG] Pair deletion at %d (op %d), metro=%.2f\n", m, state->Sq[m], *metro_factor);
#endif
        return 1;
    } else if (choice == 1 && state->q + 2 <= params->QMAX) { // INSERTION
        int m = rng_uniform_int(state->q + 1);
        int op_to_insert = rng_uniform_int(params->NOP);
        int old_q = state->q;
        state->q += 2;
        memmove(&state->Sq[m+2], &state->Sq[m], (old_q - m) * sizeof(int));
        state->Sq[m] = op_to_insert;
        state->Sq[m+1] = op_to_insert;
        *metro_factor = (double)params->NOP;
#ifdef DEBUG_UPDATES
        printf("[DEBUG] Pair insertion at %d (op %d), metro=%.2f\n", m, op_to_insert, *metro_factor);
#endif
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
#ifdef DEBUG_UPDATES
    printf("[DEBUG] Block swap at %d\n", k);
#endif
    return 1;
}

static int attempt_cycle_completion(QMCState* state, const Hamiltonian* h, const SimParams* params, double* metro_factor) {
    if (params->NCYCLES == 0) return 0;
    
    // Pick random cycle
    int cycle_idx = rng_uniform_int(params->NCYCLES);
    bitset_t* cycle = h->cycles[cycle_idx];
    int cycle_len = bitset_count(cycle);
    if (cycle_len == 0) return 0;

    // Pick random position m in Sq to insert
    int m = rng_uniform_int(state->q + 1);

    // Insert the cycle's operators at m (shuffle order)
    if (state->q + cycle_len > params->QMAX) return 0;
    int old_q = state->q;
    state->q += cycle_len;

    // Shift existing Sq to make space
    memmove(&state->Sq[m + cycle_len], &state->Sq[m], (old_q - m) * sizeof(int));

    // Add cycle ops (random order)
    int insert_pos = m;
    for (int i = 0; i < params->NOP; ++i) {
        if (bitset_get(cycle, i)) {
            state->Sq[insert_pos++] = i;
        }
    }

    // Shuffle the inserted block (Fisher-Yates)
    for (int i = m + cycle_len - 1; i > m; --i) {
        int j = m + rng_uniform_int(i - m + 1);
        int temp = state->Sq[i];
        state->Sq[i] = state->Sq[j];
        state->Sq[j] = temp;
    }

    // Apply the inserted operators to lattice (for energy recalc)
    for (int i = m; i < m + cycle_len; ++i) {
        bitset_xor(state->lattice, h->P_matrix[state->Sq[i]]);
    }

    // Metropolis factor: From paper, depends on cycle prob (e.g., 1/Ncycles * other factors)
    *metro_factor = 1.0 / (double)params->NCYCLES;  // Simplified; adjust per detailed balance

#ifdef DEBUG_UPDATES
    printf("[DEBUG] Cycle completion: idx=%d, len=%d, at m=%d, metro=%.2f\n", cycle_idx, cycle_len, m, *metro_factor);
#endif
    return 1;
}

// Worm updates (two variants from paper)
// --- Top-Level Update Algorithms ---

void do_composite_update(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    QMCState_Backup backup;
    save_state(state, &backup);
    
    ExExFloat old_weight = get_full_weight(state);
    double total_metro_factor = 1.0;
    int move_made = 0;

    double rand_val = rng_uniform_double();
    if (rand_val < 0.2) {
        move_made = attempt_classical_flip(state, params);
    } else if (rand_val < 0.4) {
        move_made = attempt_simple_swap(state);
    } else if (rand_val < 0.6) {
        move_made = attempt_pair_insertion_deletion(state, params, &total_metro_factor);
    } else if (rand_val < 0.8) {
        move_made = attempt_block_swap(state, h);
    } else {
        move_made = attempt_cycle_completion(state, h, params, &total_metro_factor);
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
        printf("[DEBUG] Proposed move: acceptance_prob=%.4f, rng=%.4f\n", acceptance_prob, rng_val_decision);
#endif

        if (rng_val_decision >= acceptance_prob) { // REJECT
            restore_state(state, &backup);
            state_recalculate_props(state, h);
            rebuild_divdiff_from_energies(state, params);
        }
    }

    free_backup(&backup);
}

void do_worm_update(QMCState* state, const Hamiltonian* h, const SimParams* params) {
    QMCState_Backup backup;
    save_state(state, &backup);
    ExExFloat old_weight = get_full_weight(state);
    double metro_factor = 1.0;
    int move_made = 0;

    if (state->worm->active) {
        // Worm is active. Choose to propagate or heal.
        double p_heal = 0.5; // Simplified probability
        if (state->worm->i == state->worm->j && state->worm->k == state->worm->l && rng_uniform_double() < p_heal) {
            // --- ATTEMPT HEAL ---
            int m = -1;
            for (int idx = 0; idx < state->q; ++idx) {
                if (state->Sq[idx] == state->worm->k) { m = idx; break; }
            }
            if (m != -1) {
                memmove(&state->Sq[m], &state->Sq[m + 1], (state->q - 1 - m) * sizeof(int));
                state->q--;
                state->worm->active = 0;
                move_made = 1;
            }
        } else {
            // --- ATTEMPT PROPAGATION ---
            if (state->q + 1 <= params->QMAX) {
                int i_new = rng_uniform_int(params->N);
                int m = rng_uniform_int(state->q + 1);
                int k_new = rng_uniform_int(params->NOP);
                memmove(&state->Sq[m + 1], &state->Sq[m], (state->q - m) * sizeof(int));
                state->Sq[m] = k_new;
                state->q++;
                bitset_xor(state->lattice, h->P_matrix[k_new]);
                state->worm->i = i_new;
                state->worm->k = k_new;
                move_made = 1;
            }
        }
    } else {
        // --- ATTEMPT CREATE ---
        if (rng_uniform_double() < 0.5) {
            if (state->q + 1 <= params->QMAX) {
                int k = rng_uniform_int(params->NOP);
                int i = rng_uniform_int(params->N);
                int m = rng_uniform_int(state->q + 1);
                memmove(&state->Sq[m + 1], &state->Sq[m], (state->q - m) * sizeof(int));
                state->Sq[m] = k;
                state->q++;
                state->worm->active = 1;
                state->worm->k = k;
                state->worm->l = k;
                state->worm->i = i;
                state->worm->j = i;
                move_made = 1;
            }
        }
    }

    if (move_made) {
        state_recalculate_props(state, h);
        rebuild_divdiff_from_energies(state, params);
        ExExFloat new_weight = get_full_weight(state);
        
        double ratio = 0.0;
        if (exex_to_double(old_weight) > 1e-9) { // Avoid division by zero
            ratio = fabs(exex_to_double(exex_divide(new_weight, old_weight)));
        }
        double acceptance_prob = fmin(1.0, ratio * metro_factor);

        if (rng_uniform_double() >= acceptance_prob) { // REJECT
            restore_state(state, &backup);
            // After restoring, we must also restore the derived properties
            state_recalculate_props(state, h);
            rebuild_divdiff_from_energies(state, params);
        }
    }
    
    free_backup(&backup);
}
