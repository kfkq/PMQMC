import sys
from pathlib import Path

# --- Add src directory to Python path for local imports ---
# (This assumes the script is run from the project root, e.g., python examples/2d_xy_model.py)
project_root = Path(__file__).resolve().parent.parent
src_path = project_root / 'scripts'
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from preprocess import OpSum, preprocess

def main():
    """Example: 4x4 2D XY Model with a Z-field pre-processing."""
    print("--- Building Hamiltonian for 4x4 2D XY Model ---")
    
    # --- Define model parameters ---
    L = 6  # The lattice is L x L
    N_QUBITS = L * L
    J_coupling = 1.0  # Ferromagnetic XY interaction
    h_field = -0.8     # Magnetic field strength in the Z direction

    # --- Build the Hamiltonian using OpSum ---
    op_sum = OpSum()

    # Iterate over every site in the 2D lattice
    for i in range(N_QUBITS):
        # Map the 1D index 'i' to 2D coordinates (x, y)
        x = i % L
        y = i // L

        # 1. Add interaction terms with neighbors (with periodic boundaries)
        # To avoid double counting, we only add bonds to the "right" and "down" neighbors.
        
        # Right neighbor
        j_right = y * L + (x + 1) % L
        op_sum.add(-J_coupling, i, 'X', j_right, 'X')
        op_sum.add(-J_coupling, i, 'Y', j_right, 'Y')

        # Down neighbor
        j_down = ((y + 1) % L) * L + x
        op_sum.add(-J_coupling, i, 'X', j_down, 'X')
        op_sum.add(-J_coupling, i, 'Y', j_down, 'Y')

        # 2. Add the magnetic field term for the current site
        op_sum.add(-h_field, i, 'Z')

    # --- Define simulation parameters for the backend solver ---
    simulation_params = {
        'beta': 2.0,
        'steps': 1_000_000,
        'steps_per_measurement': 10,
        'qmax': 500,   
        'worm_updates': False,
        
        # Specify which default observables to measure
        'default_measurements': [
            "H",
            "H2",
            "Z_MAGNETIZATION"
        ]
    }

    # Generate the hamiltonian.in file
    preprocess(op_sum, simulation_params)

if __name__ == "__main__":
    main()