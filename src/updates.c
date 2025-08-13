// File: src/updates.c
// Purpose: Implements all 7 Monte Carlo update moves.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "updates.h"
#include "utils.h"

// Warning function for when QMAX is reached
static void warn_qmax_reached(const char* update_name) {
    static int warning_count = 0;
    static const int max_warnings = 10; // Limit number of warnings printed
    
    if (warning_count < max_warnings) {
        fprintf(stderr, "Warning: QMAX reached during %s update.\n", update_name);
        warning_count++;
        if (warning_count == max_warnings) {
            fprintf(stderr, "Warning: Suppressing further QMAX warnings.\n");
        }
    }
}

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

// --- NEW HELPER FUNCTIONS FOR CYCLE COMPLETION ---

// Checks for the absence of repetitions in a sequence of length r
static int no_repetition_check(const int* sequence, int r) {
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < i; ++j) {
            if (sequence[j] == sequence[i]) {
                return 0; // Repetition found
            }
        }
    }
    return 1; // No repetitions
}

// Finds cycles containing all operators in a subsequence and returns one randomly.
// Also returns the total count of found cycles via the found_cycles_count pointer.
static int find_random_cycle(const int* subseq, int r, const Hamiltonian* h, const SimParams* params, int lmin, int lmax, int* found_cycles_count) {
    if (params->NCYCLES == 0) {
        *found_cycles_count = 0;
        return -1;
    }

    bitset_t* candidates = bitset_create(params->NCYCLES);
    // Start with all cycles as candidates
    for (int i = 0; i < params->NCYCLES; ++i) bitset_set(candidates, i);

    // Filter by ANDing with the P_in_cycles for each operator in the subsequence
    // NOTE: This requires a P_in_cycles table, which is not in your current hamiltonian.h
    // We will build it on the fly here for now. A more optimized approach would pre-compute it.
    for (int i = 0; i < r; ++i) {
        int op_idx = subseq[i];
        bitset_t* p_in_cycles_mask = bitset_create(params->NCYCLES);
        for (int c_idx = 0; c_idx < params->NCYCLES; ++c_idx) {
            if (bitset_get(h->cycles[c_idx], op_idx)) {
                bitset_set(p_in_cycles_mask, c_idx);
            }
        }
        bitset_and(candidates, p_in_cycles_mask);
        bitset_free(p_in_cycles_mask);
    }

    // Further filter by cycle length
    for (int i = 0; i < params->NCYCLES; ++i) {
        if (bitset_get(candidates, i)) {
            int len = bitset_count(h->cycles[i]);
            if (len < lmin || len > lmax) {
                // Flip the bit to 0 to remove it from candidates
                bitset_flip(candidates, i);
            }
        }
    }

    *found_cycles_count = bitset_count(candidates);
    if (*found_cycles_count == 0) {
        bitset_free(candidates);
        return -1;
    }

    // Pick a random one from the valid candidates
    int choice = rng_uniform_int(*found_cycles_count);
    int selected_cycle_idx = -1;
    int count = 0;
    for (int i = 0; i < params->NCYCLES; ++i) {
        if (bitset_get(candidates, i)) {
            if (count == choice) {
                selected_cycle_idx = i;
                break;
            }
            count++;
        }
    }

    bitset_free(candidates);
    return selected_cycle_idx;
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
    } else if (choice == 1 && state->q + 2 > params->QMAX) {
        // Warn when QMAX would be exceeded
        warn_qmax_reached("pair insertion");
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
    // This is a simplified version of the C++ code's cycle completion with gaps (u=0).
    // The full version uses a geometric distribution for 'u' (gaps).
    const int u = 0; 
    
    // Define cycle search parameters (as in mainQMC.hpp)
    // A more robust implementation would pre-calculate these and store them.
    const int cycle_min_len = 2; // Simplification
    const int cycle_max_len = params->NOP; // Simplification
    const int rmin = (cycle_min_len - 1) / 2;
    const int rmax = (cycle_max_len + 1) / 2;

    if (state->q < u + rmin) return 0;

    // 1. Pick subsequence length 'r'
    int inv_pr_A = fmin(rmax, state->q - u) - rmin + 1;
    if (inv_pr_A <= 0) return 0;
    int r = rng_uniform_int(inv_pr_A) + rmin;

    // 2. Pick subsequence S
    int m = rng_uniform_int(state->q - (r + u) + 1);
    int* Sq_subseq = &state->Sq[m]; // This is our S

    // For simplicity, we are not shuffling here as the C++ code does.
    // The shuffle is mainly for the 'gaps' part (u>0).
    if (!no_repetition_check(Sq_subseq, r)) return 0;

    // 3. Find a suitable cycle
    int lmin = 2 * r - 1;
    int lmax = 2 * r + 1;
    int found_cycles_A = 0;
    int cycle_idx = find_random_cycle(Sq_subseq, r, h, params, lmin, lmax, &found_cycles_A);

    if (cycle_idx == -1) return 0;

    // 4. Propose the move
    bitset_t* S_mask = bitset_create(params->NOP);
    for (int i = 0; i < r; ++i) bitset_set(S_mask, Sq_subseq[i]);

    bitset_t* P_mask = bitset_create(params->NOP);
    bitset_copy(P_mask, h->cycles[cycle_idx]);
    bitset_xor(P_mask, S_mask); // P = C XOR S

    int p = bitset_count(P_mask);
    if (state->q + p - r > params->QMAX) {
        bitset_free(S_mask);
        bitset_free(P_mask);
        warn_qmax_reached("cycle completion");
        return 0;
    }

    // Find cycles for the reverse move (for the Metropolis factor)
    int* S_prime = malloc(p * sizeof(int));
    int s_prime_idx = 0;
    for(int i=0; i<params->NOP; ++i) if(bitset_get(P_mask, i)) S_prime[s_prime_idx++] = i;
    
    int lmin_rev = 2 * p - 1;
    int lmax_rev = 2 * p + 1;
    int found_cycles_B = 0;
    find_random_cycle(S_prime, p, h, params, lmin_rev, lmax_rev, &found_cycles_B);
    if (found_cycles_B == 0) { // Reverse move is impossible
        free(S_prime);
        bitset_free(S_mask);
        bitset_free(P_mask);
        return 0; // Cannot guarantee detailed balance
    }

    // 5. Apply the move to Sq
    int old_q = state->q;
    int new_q = old_q + p - r;
    int* new_Sq = malloc(new_q * sizeof(int));
    
    // Copy part before subsequence
    memcpy(new_Sq, state->Sq, m * sizeof(int));
    // Copy the new operators (S')
    memcpy(&new_Sq[m], S_prime, p * sizeof(int));
    // Copy part after subsequence
    memcpy(&new_Sq[m + p], &state->Sq[m + r], (old_q - m - r) * sizeof(int));
    
    // Replace old Sq with new one
    memcpy(state->Sq, new_Sq, new_q * sizeof(int));
    state->q = new_q;
    free(new_Sq);

    // 6. Calculate the full Metropolis factor
    double wfactor = (double)found_cycles_A / (double)found_cycles_B;
    wfactor *= state->factorials[p] / state->factorials[r];
    
    int inv_pr_B = fmin(rmax, state->q - u) - rmin + 1;
    wfactor *= (double)inv_pr_A / (double)inv_pr_B;
    
    *metro_factor = wfactor;

    // Cleanup
    free(S_prime);
    bitset_free(S_mask);
    bitset_free(P_mask);

#ifdef DEBUG_UPDATES
    printf("[DEBUG] Cycle completion: r=%d, p=%d, wfactor=%.4f\n", r, p, wfactor);
#endif
    return 1; // Move was successfully proposed
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
            } else {
                warn_qmax_reached("worm propagation");
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
            } else {
                warn_qmax_reached("worm creation");
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
