### Pseudocode and Algorithm Structure of the PMR-QMC Spin-1/2 Implementation

Based on the provided merged codebase (from the GitHub repo LevBarash/PMRQMC, corresponding to the 2024 paper by Barash et al.), I'll outline the algorithm structure and pseudocode. This is the core spin-1/2 PMR-QMC implementation. The 2025 paper by Ezzell and Hen extends it for advanced measurements (e.g., custom/dynamic observables), but the base structure (PMR decomposition, cycle computation, updates, weights) remains the same. Their code [41] is not explicitly linked in the search results, but it's described as building on this repo (e.g., adding estimators to existing PMR-QMC codes). If it's public by August 2025, it would likely be at a repo like github.com/NicEzzell/PMR-QMC or similar—short runs suggest focusing on custom observable logic in observables.hpp (generated similarly to hamiltonian.hpp).

The codebase is C++ (not pure C), modular, with MPI for parallelism. Key flow: Preprocess Hamiltonian → Simulate MC chain → Post-process stats.

#### Overall Algorithm Structure
1. **Preprocessing (prepare.cpp)**: Convert input Hamiltonian (Pauli strings in text file) to PMR form. Compute fundamental cycles via GF(2) nullspace. Generate `hamiltonian.hpp` with arrays (e.g., permutations, diagonals, cycles).
   
2. **Initialization (mainQMC.hpp, PMRQMC.cpp/PMRQMC_mpi.cpp)**:
   - Load parameters (beta, steps, etc.) from `parameters.hpp`.
   - Init RNG, divided differences (divdiff.hpp for weights).
   - Load or generate initial config (random spin state |z>).
   - Compute initial weight.

3. **Simulation Loop**:
   - Thermalization: Tsteps updates.
   - Measurements: steps updates, measure every stepsPerMeasurement.
   - Updates: Composite (swaps, pair ins/del, cycle completion) or worm (single ops + healing).
   - Weights: Updated via divided differences on energy multiset.

4. **Measurements**: During loop, compute observables (e.g., <H>, <H^2>) using current config/weight.

5. **Post-Processing (datasummary.cpp or in main)**: Binning/jackknife for means/stddevs. Derived observables (e.g., specific heat).

6. **Parallelism (PMRQMC_mpi.cpp)**: MPI for independent runs, reduce stats.

7. **Checkpointing**: Save/load state on SIGTERM or resume.

Edge cases: Sign problem (use abs weights), qmax hit (warn), unused cycles (warn).

#### High-Level Pseudocode
```
# Preprocess (prepare.bin)
Read H.txt (Pauli strings: coeff q1 sigma1 q2 sigma2 ...)
Group by X-strings (permutations P), diagonals D (Z-strings)
Compute fundamental cycles: Gaussian elim over GF(2) on P matrix, minimize lengths
Generate hamiltonian.hpp: Defines N (qubits), Nop (off-diag terms), Ncycles, arrays P_matrix[], cycles[], D_coeff[][], etc.

# Main Simulation (PMRQMC.bin or MPI version)
Include hamiltonian.hpp, parameters.hpp

Init:
  RNG seed
  Divided diff objects (for exp(-beta [E0..Eq]))
  Random spin state |z> (bitset)
  q = 0
  Compute initial Energies, weight = exp(-beta E_z) * factorial terms

Thermalization:
  For Tsteps:
    Update()  # See below

Measurements:
  For measurements:
    For stepsPerMeasurement:
      Update()
    Measure()  # Compute observables, add to bins (with sign)

Update() (composite version):
  While True (until break prob or non-zero weight):
    Random choice: swap, pair insert/delete, cycle completion
    Propose change to Sq (operator sequence)
    Update Energies, currD (matrix elements)
    New weight = UpdateWeight() via divided diff add/remove
    Accept via Metropolis (detailed balance, factors for cycle probs)
  If reject all, revert

Cycle Completion:
  Pick random r (subseq length), u (gaps)
  Pick subseq from Sq, shuffle
  Find cycles containing subseq (bitset intersect)
  Pick random cycle, insert complement (S'), shuffle
  Adjust q, update weight

Worm Update (alternative):
  Introduce "worm head" (add/remove single op, break identity)
  While not identity:
    Local updates or single op moves
    Heal back to identity
  Accept intermediates via detailed balance

Measure():
  For each observable (H, H2, etc.):
    R = formula using divided diffs, Energies, currD
    Add R * sgn(weight) to bin

Post-Process:
  Binning: Means, variances from bins
  Jackknife for derived (e.g., Cv = beta^2 (<H2> - <H>^2))
  Error estimation with sign problem correction (<O> = <O sgn> / <sgn>)

Output: Means, stddevs, warnings (qmax, unused cycles)
```

#### Detailed Pseudocode for Key Components

##### 1. Preprocessing (prepare.cpp)
```
def data_extract(file):
  Parse lines: coeff, [q1 sigma1 q2 sigma2 ...]
  Return list of (coeff, [ints])  # sigma: 1=X,2=Y,3=Z

def PZcomp(data):
  For each term:
    coeff = data.coeff
    Extract Zs (qubits with Z or Y)
    Bitstring for Ps (1 at qubit for X/Y)
    If Y: coeff *= 1j, add Z qubit
  Group by Ps bitstring
  Sum coeffs for same Z under same P
  Return Ps, coeffs, Zs, Z_track (map Z to P)

def Null2(matrix):  # GF(2) nullspace
  Gaussian elim to row echelon
  Find basis for nullspace
  Minimize lengths via additions

Main:
  data = data_extract(H.txt)
  PZ = PZcomp(data)
  Ps_binary = bit_to_intvec(PZ.Ps[1:])  # Skip identity if present
  cycles = Null2(Ps_binary)
  While cycle_minimize(cycles): pass  # Add to shorten

  Output to hamiltonian.hpp:
    #define N no_qubit
    #define Nop len(Ps_nontrivial)
    #define Ncycles len(cycles)
    bitset<N> P_matrix[Nop] = {...}
    bitset<Nop> cycles[Ncycles] = {...}
    D0_size, D0_coeff[], D0_product[] (diagonals)
    D_size[Nop], D_coeff[Nop][D_maxsize], D_product[Nop][D_maxsize] (off-diags)
```

##### 2. Initialization (mainQMC.hpp)
```
init_rng(): Seed mt19937, uniform/geometric dists

init_basic():
  Precompute beta_pow_factorial[q] = (-beta)^q / q!
  factorial[q]

init():
  Precompute cycle_len[], P_in_cycles[]
  Random lattice (bitset for |z>)
  q=0, currWeight = GetWeight()  # exp(-beta E_z)
  cycles_used[] = 0 or 1 (worm)
  valid_observable[] for enabled measures
```

##### 3. Weight Computation (divdiff.hpp, mainQMC.hpp)
```
class divdiff:  # Extended float for precision
  AddElement(znew): Update Newton form divided diffs for exp(-beta [z0..zq])
  RemoveElement(): Reverse add
  UpdateWeight(): Recompute Energies[] via ApplyOperator, calc_d
    Adjust divdiffs add/remove
    Return divdiffs[q] * beta_pow_factorial[q] * Re/Abs(currD)
```

##### 4. Update (mainQMC.hpp)
Composite (default):
```
While not break (prob 0.9) and weight !=0:
  Random v:
    if v<0.25 and q>=2: Swap random adjacent ops
    elif v<0.5 and q>=2: Delete pair if equal
    elif v<0.75 and q+2<qmax: Insert pair
    else: Cycle completion
      Pick r (rmin..rmax), u (geometric)
      Pick/shuffle subseq of r+u
      Find cycles containing subseq
      Pick random cycle, insert complement S', shuffle
  Propose newWeight
Accept via Metropolis(R new/old)
Else revert
```

Cycle Completion Detail:
```
PickSubsequence(r+u): Random contiguous from Sq
FindCycles(r): Bitset intersect P_in_cycles for subseq ops
  Filter length lmin..lmax
Insert: Shift Sq, add S' (cycle - subseq), add gaps, shuffle
```

Worm (alternative):
```
While not identity:
  Local update or single op add/del
  Accept intermediate via detailed balance
  Reject whole with prob p_f=0.1
Heal to identity
```

##### 5. Measurements (mainQMC.hpp)
```
measure():
  GetWeight()  # Ensure current
  sgn = sign(Re/Abs(currD))
  For each observable:
    R = formula (e.g., H: z_q/-beta + (div[q-1]/div[q])*q/-beta)
    Add R*sgn to bin
```

##### 6. Post-Processing
```
process_single_run():
  Bin means/vars for sgn, O
  Correct for sign: <O> = <O sgn>/<sgn> with cov terms
  Jackknife for derived (e.g., Cv)

In MPI: Gather/reduce bins, compute global
```

This covers the core. For 2025 extensions: Add custom observables via similar preprocessing (PZcomp for O.txt), measure via MD0, MD calc in loop. If needed, I can refine.