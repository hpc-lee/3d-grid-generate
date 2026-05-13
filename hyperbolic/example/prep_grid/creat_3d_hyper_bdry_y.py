import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_by
from bdry_operations import arc_stretch_2d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    flag_printf = 1
    flag_topo_y = 1

    nx = 400
    nz = 300

    dx = 10.0
    dz = 10.0
    origin_x = 0.0
    origin_y = 0.0
    origin_z = 0.0

    # by shape is (nz, nx, 3) as per MATLAB
    by = np.zeros((nz, nx, 3))

    for k in range(nz):
        for i in range(nx):
            by[k, i, 0] = origin_x + i * dx
            by[k, i, 1] = origin_y
            by[k, i, 2] = origin_z + k * dz

    if flag_topo_y:
        point_x = origin_x + (nx // 2) * dx
        point_z = origin_z + (nz // 2) * dz
        L = 0.3 * nx * dx
        H = 0.2 * nx * dx
        for k in range(nz):
            for i in range(nx):
                r1 = np.sqrt((by[k, i, 0] - point_x)**2 + (by[k, i, 2] - point_z)**2)
                topo = 0.0
                if r1 < L:
                    topo = 0.5 * H * (1 + np.cos(np.pi * r1 / L))
                by[k, i, 1] = by[k, i, 1] + topo

    A = 0.00001
    by = arc_stretch_2d(A, by)

    file_name = '../data_file_3d.txt'
    export_bdry_1_by(by, nx, nz, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nz, nx) // 40)
        # MATLAB: plot3 with permuted axes + original
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        for k in range(0, nz, stride):
            ax.plot(by[k, :, 0], by[k, :, 1], by[k, :, 2], 'b-', linewidth=0.5)
        for i in range(0, nx, stride):
            ax.plot(by[:, i, 0], by[:, i, 1], by[:, i, 2], 'b-', linewidth=0.5)
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('by (both directions)')
        plt.tight_layout()
        plt.savefig('../plot/hyper_bdry_y.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/hyper_bdry_y.png')

if __name__ == "__main__":
    main()
