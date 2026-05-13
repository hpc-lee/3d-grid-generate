import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_by
from bdry_operations import extend_abs_layer, arc_stretch_2d

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    flag_printf = 1
    flag_topo_y = 1
    num_pml = 0
    nx1 = 400
    nz1 = 200
    nx = nx1 + 2 * num_pml
    nz = nz1 + 2 * num_pml
    dx = 100
    dz = 100
    ny = 300
    dy = 100
    origin_x = 0
    origin_y = 0
    origin_z = 0

    # by1 and by2 are (nz, nx, 3)
    by1 = np.zeros((nz, nx, 3))
    by2 = np.zeros((nz, nx, 3))

    for k in range(nz1):  # k 0..nz1-1 (MATLAB 1..nz1)
        for i in range(nx1):  # i 0..nx1-1 (MATLAB 1..nx1)
            by2[k + num_pml, i + num_pml, 0] = origin_x + i * dx
            by2[k + num_pml, i + num_pml, 1] = origin_y
            by2[k + num_pml, i + num_pml, 2] = origin_z + k * dz

    if flag_topo_y:
        point_x = origin_x + (nx1 // 2) * dx
        point_z = origin_z + (nz1 // 2) * dz
        L = 0.2 * nx * dx
        H = 0.2 * ny * dy
        for k in range(nz1):
            for i in range(nx1):
                r1 = np.sqrt(
                    (by2[k + num_pml, i + num_pml, 0] - point_x)**2 +
                    (by2[k + num_pml, i + num_pml, 2] - point_z)**2
                )
                topo = 0.0
                if r1 < L:
                    topo = 0.5 * H * (1 + np.cos(np.pi * r1 / L))
                by2[k + num_pml, i + num_pml, 1] -= topo

    # Extend PML
    # Note: extend_abs_layer expects (nx, ny, 3), but here we have (nz, nx, 3)
    # Wait, MATLAB's extend_abs_layer is called with (by2, dx, dz, nx, nz, num_pml)
    # So in Python, let's swap the parameters to match the dimensions:
    # In bdry_operations.py, extend_abs_layer's parameters are (bz, dx, dy, nx, ny, num_pml)
    # So here, bz is (nz, nx, 3), dx=dz, dy=dx, nx=nz, ny=nx, num_pml=num_pml
    by2 = extend_abs_layer(by2, dz, dx, nz, nx, num_pml)

    for k in range(nz):
        for i in range(nx):
            by1[k, i, 0] = by2[k, i, 0]
            by1[k, i, 1] = origin_y - (ny - 1) * dy
            by1[k, i, 2] = by2[k, i, 2]

    A = 0.000001
    by2 = arc_stretch_2d(A, by2)

    # Export by2 using export_bdry_1_by (shape (nz, nx, 3), nx and nz as parameters)
    file_name = '../data_file_3d.txt'
    export_bdry_1_by(by2, nx, nz, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nz, nx) // 40)
        # Figure 1: by1 + by2 (plot3 style)
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        for k in range(0, nz, stride):
            ax.plot(by1[k, :, 0], by1[k, :, 1], by1[k, :, 2], 'b-', linewidth=0.5)
            ax.plot(by2[k, :, 0], by2[k, :, 1], by2[k, :, 2], 'b-', linewidth=0.5)
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('by1 + by2')
        plt.tight_layout()
        plt.savefig('../plot/para_bdry_y_plot3.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        # Figure 2: by2 surf
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        X = by2[:, :, 0]; Y = by2[:, :, 1]; Z = by2[:, :, 2]
        ax.plot_surface(X[::stride, ::stride], Y[::stride, ::stride],
                        Z[::stride, ::stride], cmap='jet', edgecolor='none')
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('by2 surface')
        plt.tight_layout()
        plt.savefig('../plot/para_bdry_y_surf.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/para_bdry_y_*.png')


if __name__ == "__main__":
    main()
