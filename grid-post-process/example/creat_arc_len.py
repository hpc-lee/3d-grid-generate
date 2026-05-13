import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'plotting'))
from export_c import export_arc_len
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

flag_printf = 1
nz = 201
num_of_step = nz - 1
flag_flip = 1

step = np.zeros(num_of_step)
for i in range(num_of_step):
    step[i] = 1.0

# # Gradient Grid
# incre_layer = 12
# max_ratio = 3
# incre_ratio = np.exp((np.log(max_ratio)/incre_layer))
#
# for i in range(10):
#     step[i] = 1.0
#
# for i in range(10, 10+incre_layer):
#     step[i] = step[i-1]*incre_ratio
#
# for i in range(10+incre_layer, num_of_step):
#     step[i] = step[10+incre_layer-1]

if flag_flip:
    step = np.flip(step)

sum_step = np.sum(step)

# normalization step
step_nor = step / sum_step

arc_len = np.zeros(nz)
arc_len[0] = 0.0
for i in range(1, nz):
    arc_len[i] = arc_len[i-1] + step_nor[i-1]

if (arc_len[nz-1] - 1) > 1e-8:
    raise ValueError("step set is error, please check and reset")

# create step file
file_name = './arc_len_file1.txt'
export_arc_len(arc_len, nz, file_name)

if flag_printf:
    os.makedirs('../plot', exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(1, nz + 1), arc_len, 'b-', linewidth=1.5)
    ax.set_xlabel('Point index', fontweight='bold')
    ax.set_ylabel('Arc length', fontweight='bold')
    ax.set_title('Arc length distribution', fontweight='bold')
    ax.grid(True)
    plt.tight_layout()
    plt.savefig('../plot/arc_len.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved: ../plot/arc_len.png')
