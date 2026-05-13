import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_step

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    flag_printf = 1
    nz = 201
    num_of_step = nz - 1
    flag_flip = 1

    # Gradient Grid
    incre_layer = 12
    max_ratio = 3
    incre_ratio = np.exp((np.log(max_ratio) / incre_layer))

    step = np.zeros(num_of_step)
    for i in range(10):  # i 0-9 (MATLAB 1-10)
        step[i] = 1
    for i in range(10, 10 + incre_layer):  # i 10-21 (MATLAB 11-22)
        step[i] = step[i-1] * incre_ratio
    for i in range(10 + incre_layer, num_of_step):  # MATLAB 23 to num_of_step
        step[i] = step[10 + incre_layer - 1]

    if flag_flip:
        step = np.flip(step)

    sum_step = np.sum(step)
    # normalization step
    step_nor = step / sum_step
    sum_step_nor = np.sum(step_nor)
    if (sum_step_nor - 1) > 1e-8:
        raise ValueError("step set is error, please check and reset")

    # export step file
    file_name = '../step_file_3d.txt'
    export_step(step_nor, num_of_step, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(range(1, num_of_step + 1), step_nor, 'b-', linewidth=1.5)
        ax.set_xlabel('Step index', fontweight='bold')
        ax.set_ylabel('Normalized step size', fontweight='bold')
        ax.set_title('Step distribution', fontweight='bold')
        ax.grid(True)
        plt.tight_layout()
        plt.savefig('../plot/para_step.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/para_step.png')


if __name__ == "__main__":
    main()
