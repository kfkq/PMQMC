// File: src/hamiltonian.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hamiltonian.h"

static char* trim_whitespace(char *str) {
    while (isspace((unsigned char)*str)) str++;
    if (*str == 0) return str;
    char *end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return str;
}

Hamiltonian* hamiltonian_create_and_load(const char* filename, SimParams* params) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening input file");
        return NULL;
    }

    // Allocate the main Hamiltonian struct
    Hamiltonian* h = (Hamiltonian*)calloc(1, sizeof(Hamiltonian));
    if (!h) {
        perror("Failed to allocate Hamiltonian struct");
        fclose(fp);
        return NULL;
    }

    char line_buffer[8192];

    // First Pass: Read SIMULATION_PARAMS
    while (fgets(line_buffer, sizeof(line_buffer), fp)) {
        if (strncmp(line_buffer, "SIMULATION_PARAMS_BEGIN", 23) == 0) {
            while (fgets(line_buffer, sizeof(line_buffer), fp) && strncmp(line_buffer, "SIMULATION_PARAMS_END", 21) != 0) {
                char key[100];
                if (sscanf(line_buffer, " %s", key) != 1) continue;
                if (strcmp(key, "N") == 0) sscanf(line_buffer, "%*s %d", &params->N);
                if (strcmp(key, "NOP") == 0) sscanf(line_buffer, "%*s %d", &params->NOP);
                if (strcmp(key, "NCYCLES") == 0) sscanf(line_buffer, "%*s %d", &params->NCYCLES);
                if (strcmp(key, "BETA") == 0) sscanf(line_buffer, "%*s %lf", &params->BETA);
                if (strcmp(key, "TSTEPS") == 0) sscanf(line_buffer, "%*s %lld", &params->TSTEPS);
                if (strcmp(key, "STEPS") == 0) sscanf(line_buffer, "%*s %lld", &params->STEPS);
                if (strcmp(key, "STEPS_PER_MEASUREMENT") == 0) sscanf(line_buffer, "%*s %d", &params->STEPS_PER_MEASUREMENT);
                if (strcmp(key, "QMAX") == 0) sscanf(line_buffer, "%*s %d", &params->QMAX);
                if (strcmp(key, "NBINS") == 0) sscanf(line_buffer, "%*s %d", &params->NBINS);
                if (strcmp(key, "WORM") == 0) {
                    char worm_str[10];
                    sscanf(line_buffer, "%*s %s", worm_str);
                    params->WORM = (strcmp(worm_str, "True") == 0 || strcmp(worm_str, "true") == 0 || strcmp(worm_str, "1") == 0) ? 1 : 0;
                }
            }
            break;
        }
    }

    // Allocate Memory
    if (params->N <= 0) {
        fprintf(stderr, "Error: N must be > 0.\n");
        hamiltonian_free(h, params); // Use the new free function for cleanup
        fclose(fp);
        return NULL;
    }
    if (params->NOP > 0) h->P_matrix = calloc(params->NOP, sizeof(bitset_t*));
    if (params->NCYCLES > 0) h->cycles = calloc(params->NCYCLES, sizeof(bitset_t*));
    if (params->NOP > 0) {
        h->D_sizes = calloc(params->NOP, sizeof(int));
        h->D_coeffs = calloc(params->NOP, sizeof(complex_t*));
        h->D_products = calloc(params->NOP, sizeof(bitset_t**));
    }

    // Second Pass: Read data
    rewind(fp);
    while (fgets(line_buffer, sizeof(line_buffer), fp)) {
        char* line = trim_whitespace(line_buffer);
        if (line[0] == '#' || line[0] == '\0') continue;

        if (strcmp(line, "P_MATRIX_BEGIN") == 0) {
            for (int i = 0; i < params->NOP; ++i) {
                if (!fgets(line_buffer, sizeof(line_buffer), fp)) {
                    fprintf(stderr, "Error: Unexpected end of file in P_MATRIX block.\n");
                    break;
                }
                h->P_matrix[i] = bitset_create_from_string(trim_whitespace(line_buffer));
            }
        } else if (strcmp(line, "CYCLES_BEGIN") == 0) {
            for (int i = 0; i < params->NCYCLES; ++i) {
                if (!fgets(line_buffer, sizeof(line_buffer), fp)) {
                     fprintf(stderr, "Error: Unexpected end of file in CYCLES block.\n");
                    break;
                }
                h->cycles[i] = bitset_create_from_string(trim_whitespace(line_buffer));
            }
        }
        // ... (rest of the parsing logic for DIAGONAL and OFF_DIAGONAL is the same as before) ...
        else if (strcmp(line, "DIAGONAL_TERM_BEGIN") == 0) {
            while (fgets(line_buffer, sizeof(line_buffer), fp) && strncmp(trim_whitespace(line_buffer), "DIAGONAL_TERM_END", 17) != 0) {
                char* current_line = trim_whitespace(line_buffer);
                if (strncmp(current_line, "SIZE", 4) == 0) {
                    sscanf(current_line, "%*s %d", &h->D0_size);
                    if (h->D0_size > 0) {
                        h->D0_coeff = malloc(h->D0_size * sizeof(complex_t));
                        h->D0_product = malloc(h->D0_size * sizeof(bitset_t*));
                    }
                } else if (strcmp(current_line, "COEFFS_BEGIN") == 0) {
                    for (int i = 0; i < h->D0_size; ++i) {
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        double real, imag;
                        sscanf(line_buffer, "%lf %lf", &real, &imag);
                        h->D0_coeff[i] = real + imag * I;
                    }
                } else if (strcmp(current_line, "PRODUCTS_BEGIN") == 0) {
                    for (int i = 0; i < h->D0_size; ++i) {
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        h->D0_product[i] = bitset_create_from_string(trim_whitespace(line_buffer));
                    }
                }
            }
        } else if (strcmp(line, "OFF_DIAGONAL_TERMS_BEGIN") == 0) {
             while (fgets(line_buffer, sizeof(line_buffer), fp) && strncmp(trim_whitespace(line_buffer), "OFF_DIAGONAL_TERMS_END", 22) != 0) {
                char* current_line = trim_whitespace(line_buffer);
                if (strcmp(current_line, "SIZES_BEGIN") == 0) {
                    fgets(line_buffer, sizeof(line_buffer), fp);
                    char* token = strtok(line_buffer, " \t\n");
                    for (int i = 0; i < params->NOP; ++i) {
                        h->D_sizes[i] = atoi(token);
                        token = strtok(NULL, " \t\n");
                    }
                } else if (strcmp(current_line, "COEFFS_BEGIN") == 0) {
                    for (int i = 0; i < params->NOP; ++i) {
                        h->D_coeffs[i] = malloc(h->D_sizes[i] * sizeof(complex_t));
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        char* token_real = strtok(line_buffer, " \t\n");
                        for (int j = 0; j < h->D_sizes[i]; ++j) {
                            char* token_imag = strtok(NULL, " \t\n");
                            h->D_coeffs[i][j] = atof(token_real) + atof(token_imag) * I;
                            token_real = strtok(NULL, " \t\n");
                        }
                    }
                } else if (strcmp(current_line, "PRODUCTS_BEGIN") == 0) {
                     for (int i = 0; i < params->NOP; ++i) {
                        h->D_products[i] = malloc(h->D_sizes[i] * sizeof(bitset_t*));
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        char* token = strtok(line_buffer, ";\n");
                        for (int j = 0; j < h->D_sizes[i]; ++j) {
                            h->D_products[i][j] = bitset_create_from_string(trim_whitespace(token));
                            token = strtok(NULL, ";\n");
                        }
                    }
                }
            }
        }
    }

    fclose(fp);
    return h;
}

void hamiltonian_free(Hamiltonian* h, const SimParams* params) {
    if (!h || !params) return;

    if (h->P_matrix) {
        for (int i = 0; i < params->NOP; ++i) bitset_free(h->P_matrix[i]);
        free(h->P_matrix);
    }
    if (h->cycles) {
        for (int i = 0; i < params->NCYCLES; ++i) bitset_free(h->cycles[i]);
        free(h->cycles);
    }
    if (h->D0_product) {
        for (int i = 0; i < h->D0_size; ++i) bitset_free(h->D0_product[i]);
        free(h->D0_product);
    }
    free(h->D0_coeff);

    if (h->D_sizes) {
        for (int i = 0; i < params->NOP; ++i) {
            if (h->D_coeffs && h->D_coeffs[i]) free(h->D_coeffs[i]);
            if (h->D_products && h->D_products[i]) {
                for (int j = 0; j < h->D_sizes[i]; ++j) {
                    bitset_free(h->D_products[i][j]);
                }
                free(h->D_products[i]);
            }
        }
    }
    free(h->D_sizes);
    free(h->D_coeffs);
    free(h->D_products);

    // Finally, free the main struct
    free(h);
}