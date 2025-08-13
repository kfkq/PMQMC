import numpy as np
import cmath
import sys
from typing import List, Tuple, Union

class OpSum:
    """Builds a Hamiltonian as a sum of Pauli operator terms. Uses 0-based site indexing."""
    def __init__(self):
        self.terms: List[Tuple[float, List[Tuple[int, str]]]] = []
    
    def add(self, coeff: float, *ops: Union[int, str]):
        """Adds a term to the Hamiltonian. Sites are 0-indexed."""
        if len(ops) % 2 != 0:
            raise ValueError("Operators must come in site-pauli pairs")
        term = []
        for i in range(0, len(ops), 2):
            site = ops[i]
            pauli = ops[i+1].upper()
            if not isinstance(site, int) or site < 0:
                raise ValueError(f"Invalid site: {site}. Sites must be 0-indexed non-negative integers.")
            if pauli not in ['X', 'Y', 'Z']:
                raise ValueError(f"Invalid Pauli: {pauli}")
            term.append((site, pauli))
        self.terms.append((coeff, term))
    
    def __iadd__(self, other):
        """Allows adding terms with the += operator."""
        if isinstance(other, tuple):
            self.add(*other)
        else:
            raise TypeError("Can only add tuples of (coeff, site0, pauli0, ...)")
        return self

def Hamiltonian(op_sum: OpSum):
    """Converts the OpSum object to the specialized PMR-QMC format."""
    data = [(term[0], [item for sublist in term[1] for item in (sublist[0], sublist[1])]) for term in op_sum.terms]
    
    max_site = max((max(site for site, pauli in term_ops) if term_ops else -1 for _, term_ops in op_sum.terms), default=-1)
    n_qubit = max_site + 1
    
    if n_qubit == 0:
        return 0, [], [], data

    Ps = {}
    for coeff, flat_terms in data:
        terms = [(flat_terms[i], flat_terms[i+1]) for i in range(0, len(flat_terms), 2)]
        p_bit = [0] * n_qubit
        z_qubits = []
        curr_coeff = coeff
        for q, p in terms:
            if p == 'X':
                p_bit[q] = 1
            elif p == 'Y':
                p_bit[q] = 1
                curr_coeff *= -1j
                z_qubits.append(q)
            elif p == 'Z':
                z_qubits.append(q)
        p_tuple = tuple(p_bit)
        if p_tuple not in Ps:
            Ps[p_tuple] = []
        Ps[p_tuple].append((curr_coeff, sorted(z_qubits)))
    
    def bit_tuple_to_int(t):
        val = 0
        for i, bit in enumerate(t):
            if bit:
                val |= 1 << i
        return val

    sorted_ps = sorted(Ps.items(), key=lambda item: bit_tuple_to_int(item[0]))
    
    ps_list = [list(k) for k, v in sorted_ps]
    coeff_zs = [v for k, v in sorted_ps]
    return n_qubit, ps_list, coeff_zs, data

def gf2_null(A):
    """
    Finds the null space of the rows of matrix A over GF(2).
    """
    if A.size == 0:
        return np.array([])

    matrix = A.T.tolist() # Transpose to work on rows (operators)
    num_rows = len(matrix)
    if num_rows == 0:
        return np.array([])
    num_cols = len(matrix[0])
    
    matrix_RE = [row[:] for row in matrix]
    marked_rows = []

    # Gaussian Elimination part from prepare.cpp
    for j in range(num_cols):
        for i in range(num_rows):
            if matrix_RE[i][j] == 1:
                marked_rows.append(i)
                for k in range(num_cols):
                    if k != j and matrix_RE[i][k] == 1:
                        for l in range(num_rows):
                            matrix_RE[l][k] = matrix_RE[l][j] ^ matrix_RE[l][k]
                break

    marked_rows = sorted(list(set(marked_rows)))
    
    marked_matrix = [matrix_RE[i] for i in marked_rows]
    
    nullspace_basis = []
    unmarked_rows = [i for i in range(num_rows) if i not in marked_rows]

    for i in unmarked_rows:
        row_i = matrix_RE[i]
        null_vector = [0] * num_rows
        null_vector[i] = 1
        for j in range(num_cols):
            if row_i[j] == 1:
                for k_idx, k in enumerate(marked_rows):
                    if marked_matrix[k_idx][j] == 1:
                        null_vector[k] = 1
        nullspace_basis.append(null_vector)
        
    return np.array(nullspace_basis).T if nullspace_basis else np.array([])

def minimize_cycles(basis):
    """Finds a basis for the null space with the shortest possible vectors."""
    if basis.size == 0:
        return basis
    basis = basis.T
    
    # Sort by length initially
    basis = sorted(basis, key=sum)
    
    n = len(basis)
    changed = True
    while changed:
        changed = False
        for i in range(n):
            for j in range(i + 1, n):
                new_vec = [(basis[i][k] + basis[j][k]) % 2 for k in range(len(basis[i]))]
                if 1 < sum(new_vec) < sum(basis[j]):
                    basis[j] = new_vec
                    changed = True
        # Re-sort after each pass
        basis = sorted(basis, key=sum)
        
    return np.array(basis).T

def z_to_bit(z_list, n_qubit):
    """Converts a list of 0-indexed Z operator locations to a bitstring."""
    bit = [0] * n_qubit
    for q in z_list:
        bit[q] = 1
    return bit

def preprocess(op_sum: OpSum, params: dict):
    """
    Takes an OpSum object and simulation parameters, and writes them to a
    robust, block-structured input file named 'hamiltonian.in'.
    """
    # The output filename is now fixed.
    output_file = 'hamiltonian.in'
    
    # --- Core logic remains the same ---
    n_qubit, ps_list, coeff_zs, _ = Hamiltonian(op_sum)
    
    d0 = ps_list and not any(ps_list[0])
    
    d0_coeff = [c for c,z in coeff_zs[0]] if d0 else []
    d0_z = [z_to_bit(z, n_qubit) for c,z in coeff_zs[0]] if d0 else []
    
    d_coeffs = [[c for c,z in op] for op in (coeff_zs[1:] if d0 else coeff_zs)]
    d_zs = [[z_to_bit(z, n_qubit) for c,z in op] for op in (coeff_zs[1:] if d0 else coeff_zs)]
    
    ps_mat = np.array(ps_list[1:] if d0 else ps_list, dtype=int).T
    cycles = gf2_null(ps_mat)
    min_cycles = minimize_cycles(cycles)
    
    nop = len(ps_list) - (1 if d0 else 0)
    ncycles = cycles.shape[1] if cycles.size else 0

    # --- New block-structured file writing ---
    with open(output_file, 'w') as f:
        f.write("# PMR-QMC Simulation Input File\n")
        f.write("# Generated by Python preprocessor script.\n\n")

        f.write("SIMULATION_PARAMS_BEGIN\n")
        f.write(f"  N         {n_qubit}\n")
        f.write(f"  NOP       {nop}\n")
        f.write(f"  NCYCLES   {ncycles}\n")
        f.write(f"  BETA      {params['beta']}\n")
        f.write(f"  STEPS     {params['steps']}\n")
        f.write(f"  STEPS_PER_MEASUREMENT {params['steps_per_measurement']}\n")
        f.write(f"  QMAX      {params['qmax']}\n")
        f.write(f"  WORM      {params.get('worm_updates', False)}\n")
        f.write("SIMULATION_PARAMS_END\n\n")

        if 'default_measurements' in params and params['default_measurements']:
            f.write("DEFAULT_MEASUREMENTS_BEGIN\n")
            for obs in params['default_measurements']:
                f.write(f"  {obs}\n")
            f.write("DEFAULT_MEASUREMENTS_END\n\n")

        f.write("P_MATRIX_BEGIN\n")
        for p in ps_list[1:] if d0 else ps_list:
            f.write(''.join(map(str, reversed(p))) + '\n')
        f.write("P_MATRIX_END\n\n")

        f.write("CYCLES_BEGIN\n")
        if cycles.size > 0:
            for i in range(ncycles):
                f.write(''.join(map(str, reversed(min_cycles[:, i]))) + '\n')
        f.write("CYCLES_END\n\n")

        f.write("DIAGONAL_TERM_BEGIN\n")
        f.write(f"  SIZE {len(d0_coeff)}\n")
        f.write("  COEFFS_BEGIN\n")
        for c in d0_coeff:
            f.write(f"    {c.real} {c.imag}\n")
        f.write("  COEFFS_END\n")
        f.write("  PRODUCTS_BEGIN\n")
        for z in d0_z:
            f.write("    " + ''.join(map(str, reversed(z))) + '\n')
        f.write("  PRODUCTS_END\n")
        f.write("DIAGONAL_TERM_END\n\n")

        f.write("OFF_DIAGONAL_TERMS_BEGIN\n")
        f.write("  SIZES_BEGIN\n    ")
        f.write(' '.join(map(str, [len(op) for op in d_coeffs])) + '\n')
        f.write("  SIZES_END\n")
        f.write("  COEFFS_BEGIN\n")
        for op in d_coeffs:
            f.write("    " + ' '.join(f"{c.real} {c.imag}" for c in op) + '\n')
        f.write("  COEFFS_END\n")
        f.write("  PRODUCTS_BEGIN\n")
        for op_z in d_zs:
            f.write("    " + ';'.join(''.join(map(str, reversed(z))) for z in op_z) + '\n')
        f.write("  PRODUCTS_END\n")
        f.write("OFF_DIAGONAL_TERMS_END\n")

    print(f"Successfully created input file: {output_file}")