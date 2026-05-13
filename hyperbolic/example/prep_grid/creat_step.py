import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_step
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    flag_printf = 1
    nz = 201
    num_of_step = nz - 1
    dh = -15.0
    step_vals = np.full(num_of_step, dh)

    file_name = '../step_file_3d.txt'
    export_step(step_vals, num_of_step, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(range(1, num_of_step + 1), step_vals, 'b-', linewidth=1.5)
        ax.set_xlabel('Step index', fontweight='bold')
        ax.set_ylabel('Step value', fontweight='bold')
        ax.set_title('Step distribution', fontweight='bold')
        ax.grid(True)
        plt.tight_layout()
        plt.savefig('../plot/hyper_step.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/hyper_step.png')

if __name__ == "__main__":
    main()
