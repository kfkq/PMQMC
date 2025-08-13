# PMQMC in C

### Introduction

A high-performance C implementation of the Permutation Matrix Representation Quantum Monte Carlo (PMR-QMC) algorithm for simulating arbitrary spin-1/2 systems.

The project uses a three-stage workflow for robust and reproducible scientific analysis:
1.  **Preprocessing:** A Python script converts a user-defined Hamiltonian into the required C input format (`hamiltonian.in`).
2.  **Simulation:** A fast C program (`pmqmc`) runs the simulation and saves raw data to a compressed HDF5 file (`raw_data.h5`).
3.  **Post-Processing:** A flexible Python script (`scripts/postprocess.py`) inspects the raw data, automatically determines optimal binning from the autocorrelation time, performs a Jackknife analysis, and writes the final results to `default_obs.dat`.

### Requirements

*   **C Compiler:** GCC or Clang.
*   **GNU Make:** For building the C executable.
*   **Python 3 & Libraries:** For preprocessing and analysis.
    ```bash
    pip install numpy h5py matplotlib
    ```
*   **HDF5 C Library:** Required for compiling the C simulator.
    *   On Debian/Ubuntu: `sudo apt-get install libhdf5-dev`
    *   On macOS (Homebrew): `brew install hdf5`

### Usage

**1. Configure and Preprocess**

Edit a Python script (e.g., `examples/1d_chain.py`) to define your model, simulation parameters, and which observables to measure. Then, run it:
```bash
python3 examples/1d_chain.py
```
This creates the `hamiltonian.in` file.

**2. Compile**

**Important:** Before compiling for the first time, you must edit the `Makefile` and set the `HDF5_INCLUDE_DIR` and `HDF5_LIB_DIR` variables to match the paths of your HDF5 installation.

Then, compile the C simulator:
```bash
make
```
This creates the `pmqmc` executable.

**3. Run Simulation**

Execute the simulator inside the folder where `hamiltonian.in` exist. It will automatically find `hamiltonian.in`.
```bash
../pmqmc
```
This runs the simulation and generates a compressed `raw_data.h5` file.

**4. Analyze Results**

Execute the Python analysis script inside the folder where `raw_data.h5` exist. It automatically finds the necessary files.
```bash
python3 ../scripts/postprocess.py
```

This script will print the final means and standard errors to the console and save them in `default_obs.dat`.

### HDF5 Data Structure

The output data is stored in `raw_data.h5`. You can inspect this file with tools like `h5dump` or load it easily in Python, MATLAB, etc. The internal structure is as follows:

*   **/measurements** `(Dataset)`: A table containing the time-series data. Each row has the following fields (if measured):
    *   `sgn` (double): The sign of the configuration weight.
    *   `H` (double): The instantaneous value of the energy, `<H>`.
    *   `H2` (double): The instantaneous value of `<H^2>`.
    *   `Z_mag` (double): The instantaneous value of the Z-magnetization.
    *   `q` (int): The length of the operator sequence.

*   **/measurements** `(Attributes)`: The dataset is annotated with metadata, storing the key simulation parameters used for the run (e.g., `BETA`, `N`, `STEPS`, `SKIP_MEASUREMENTS`, etc.).