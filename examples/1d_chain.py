import sys
from pathlib import Path

# --- Add src directory to Python path for local imports ---
project_root = Path(__file__).resolve().parent.parent
src_path = project_root / 'scripts'
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Now that 'src' is in the path, we can import from preprocess.py.
from preprocess import OpSum, preprocess

def main():
    """Example: 8-qubit Transverse-Field Heisenberg spin chain pre-processing."""
    print("--- Building Hamiltonian for 4-qubit Ising Model ---")
    
    # Define model parameters.
    N_QUBITS = 8
    J_coupling = 1.0
    h_field = 0.5

    # Build the Hamiltonian using OpSum.
    op_sum = OpSum()
    # Add interaction terms (-J * S_i . S_{i+1}) for sites 0, 1, 2
    for i in range(N_QUBITS):
        op_sum.add(-J_coupling, i, 'X', (i + 1) % N_QUBITS, 'X')
        op_sum.add(-J_coupling, i, 'Y', (i + 1) % N_QUBITS, 'Y')
        op_sum.add(-J_coupling, i, 'Z', (i + 1) % N_QUBITS, 'Z')
        op_sum.add(-h_field, i, 'X')

    # Define simulation parameters for the backend solver.
    simulation_params = {
    'beta': 1.0,
    'steps': 1_000_000,
    'steps_per_measurement': 10,
    'qmax': 30,
    'worm_updates': False,
    
    'default_measurements': [
        "H",
        "H2",
        "Z_MAGNETIZATION"
    ]}

    preprocess(op_sum, simulation_params)

if __name__ == "__main__":
    main()