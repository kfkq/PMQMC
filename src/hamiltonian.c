// File: src/hamiltonian.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hamiltonian.h"

// ... (global variable definitions and trim_whitespace are unchanged) ...
bitset_t** P_matrix = NULL;
bitset_t** cycles = NULL;
complex_t* D0_coeff = NULL;
bitset_t** D0_product = NULL;
int D0_size = 0;
int* D_sizes = NULL;
complex_t** D_coeffs = NULL;
bitset_t*** D_products = NULL;

static char* trim_whitespace(char *str) {
    while (isspace((unsigned char)*str)) str++;
    if (*str == 0) return str;
    char *end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return str;
}


int hamiltonian_load(const char* filename, SimParams* params) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening input file");
        return -1;
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
            }
            break;
        }
    }

    // Allocate Memory
    if (params->N <= 0) {
        fprintf(stderr, "Error: N must be > 0.\n");
        fclose(fp);
        return -1;
    }
    if (params->NOP > 0) P_matrix = calloc(params->NOP, sizeof(bitset_t*));
    if (params->NCYCLES > 0) cycles = calloc(params->NCYCLES, sizeof(bitset_t*));
    if (params->NOP > 0) {
        D_sizes = calloc(params->NOP, sizeof(int));
        D_coeffs = calloc(params->NOP, sizeof(complex_t*));
        D_products = calloc(params->NOP, sizeof(bitset_t**));
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
                P_matrix[i] = bitset_create_from_string(trim_whitespace(line_buffer));
            }
        } else if (strcmp(line, "CYCLES_BEGIN") == 0) {
            for (int i = 0; i < params->NCYCLES; ++i) {
                if (!fgets(line_buffer, sizeof(line_buffer), fp)) {
                     fprintf(stderr, "Error: Unexpected end of file in CYCLES block.\n");
                    break;
                }
                cycles[i] = bitset_create_from_string(trim_whitespace(line_buffer));
            }
        }
        // ... (rest of the parsing logic for DIAGONAL and OFF_DIAGONAL is the same as before) ...
        else if (strcmp(line, "DIAGONAL_TERM_BEGIN") == 0) {
            while (fgets(line_buffer, sizeof(line_buffer), fp) && strncmp(trim_whitespace(line_buffer), "DIAGONAL_TERM_END", 17) != 0) {
                char* current_line = trim_whitespace(line_buffer);
                if (strncmp(current_line, "SIZE", 4) == 0) {
                    sscanf(current_line, "%*s %d", &D0_size);
                    if (D0_size > 0) {
                        D0_coeff = malloc(D0_size * sizeof(complex_t));
                        D0_product = malloc(D0_size * sizeof(bitset_t*));
                    }
                } else if (strcmp(current_line, "COEFFS_BEGIN") == 0) {
                    for (int i = 0; i < D0_size; ++i) {
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        double real, imag;
                        sscanf(line_buffer, "%lf %lf", &real, &imag);
                        D0_coeff[i] = real + imag * I;
                    }
                } else if (strcmp(current_line, "PRODUCTS_BEGIN") == 0) {
                    for (int i = 0; i < D0_size; ++i) {
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        D0_product[i] = bitset_create_from_string(trim_whitespace(line_buffer));
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
                        D_sizes[i] = atoi(token);
                        token = strtok(NULL, " \t\n");
                    }
                } else if (strcmp(current_line, "COEFFS_BEGIN") == 0) {
                    for (int i = 0; i < params->NOP; ++i) {
                        D_coeffs[i] = malloc(D_sizes[i] * sizeof(complex_t));
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        char* token_real = strtok(line_buffer, " \t\n");
                        for (int j = 0; j < D_sizes[i]; ++j) {
                            char* token_imag = strtok(NULL, " \t\n");
                            D_coeffs[i][j] = atof(token_real) + atof(token_imag) * I;
                            token_real = strtok(NULL, " \t\n");
                        }
                    }
                } else if (strcmp(current_line, "PRODUCTS_BEGIN") == 0) {
                     for (int i = 0; i < params->NOP; ++i) {
                        D_products[i] = malloc(D_sizes[i] * sizeof(bitset_t*));
                        fgets(line_buffer, sizeof(line_buffer), fp);
                        char* token = strtok(line_buffer, ";\n");
                        for (int j = 0; j < D_sizes[i]; ++j) {
                            D_products[i][j] = bitset_create_from_string(trim_whitespace(token));
                            token = strtok(NULL, ";\n");
                        }
                    }
                }
            }
        }
    }

    fclose(fp);
    return 0;
}

void hamiltonian_free(const SimParams* params) {
    if (!params) return;

    if (P_matrix) {
        for (int i = 0; i < params->NOP; ++i) bitset_free(P_matrix[i]);
        free(P_matrix);
    }
    if (cycles) {
        for (int i = 0; i < params->NCYCLES; ++i) bitset_free(cycles[i]);
        free(cycles);
    }
    if (D0_product) {
        for (int i = 0; i < D0_size; ++i) bitset_free(D0_product[i]);
        free(D0_product);
    }
    free(D0_coeff);

    if (D_sizes) {
        for (int i = 0; i < params->NOP; ++i) {
            if (D_coeffs && D_coeffs[i]) free(D_coeffs[i]);
            if (D_products && D_products[i]) {
                for (int j = 0; j < D_sizes[i]; ++j) {
                    bitset_free(D_products[i][j]);
                }
                free(D_products[i]);
            }
        }
    }
    free(D_sizes);
    free(D_coeffs);
    free(D_products);
}