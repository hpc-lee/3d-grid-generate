import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_bx
from bdry_operations import extend_abs_layer, arc_stretch_2d

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    flag_printf = 1
    flag_topo_x = 1
    num_pml = 0
    ny1 = 300
    nz1 = 200
    ny = ny1 + 2 * num_pml
    nz = nz1 + 2 * num_pml
    dy = 100
    dz = 100
    nx = 400
    dx = 100
    origin_x = 0
    origin_y = 0
    origin_z = 0

    # bx1 and bx2 are (nz, ny, 3)
    bx1 = np.zeros((nz, ny, 3))
    bx2 = np.zeros((nz, ny, 3))

    for k in range(nz1):  # k 0..nz1-1
        for j in range(ny1):  # j 0..ny1-1
            bx2[k + num_pml, j + num_pml, 0] = origin_x
            bx2[k + num_pml, j + num_pml, 1] = origin_y + j * dy
            bx2[k + num_pml, j + num_pml, 2] = origin_z + k * dz

    if flag_topo_x:
        point_y = origin_y + (ny1 // 2) * dy
        point_z = origin_z + (nz1 // 2) * dz
        L = 0.2 * ny * dy
        H = 0.15 * nx * dx
        for k in range(nz1):
            for j in range(ny1):
                r1 = np.sqrt(
                    (bx2[k + num_pml, j + num_pml, 1] - point_y)**2 +
                    (bx2[k + num_pml, j + num_pml, 2] - point_z)**2
                )
                topo = 0.0
                if r1 < L:
                    topo = 0.5 * H * (1 + np.cos(np.pi * r1 / L))
                bx2[k + num_pml, j + num_pml, 0] -= topo

    # Extend PML
    # Similar to creat_bdry_y: bz is (nz, ny, 3), dx=dz, dy=dy, nx=nz, ny=ny
    bx2 = extend_abs_layer(bx2, dz, dy, nz, ny, num_pml)

    for k in range(nz):
        for j in range(ny):
            bx1[k, j, 0] = origin_x - (nx - 1) * dx
            bx1[k, j, 1] = bx2[k, j, 1]
            bx1[k, j, 2] = bx2[k, j, 2]

    A = 0.000001
    bx2 = arc_stretch_2d(A, bx2)

    # Export using export_bdry_1_bx (shape (nz, ny, 3), ny and nz as parameters)
    file_name = '../data_file_3d.txt'
    export_bdry_1_bx(bx2, ny, nz, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nz, ny) // 40)
        # Figure 1: bx1 + bx2 (plot3 style)
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        for k in range(0, nz, stride):
            ax.plot(bx1[k, :, 0], bx1[k, :, 1], bx1[k, :, 2], 'b-', linewidth=0.5)
            ax.plot(bx2[k, :, 0], bx2[k, :, 1], bx2[k, :, 2], 'b-', linewidth=0.5)
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('bx1 + bx2')
        plt.tight_layout()
        plt.savefig('../plot/para_bdry_x_plot3.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        # Figure 2: bx2 surf
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        X = bx2[:, :, 0]; Y = bx2[:, :, 1]; Z = bx2[:, :, 2]
        ax.plot_surface(X[::stride, ::stride], Y[::stride, ::stride],
                        Z[::stride, ::stride], cmap='jet', edgecolor='none')
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('bx2 surface')
        plt.tight_layout()
        plt.savefig('../plot/para_bdry_x_surf.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/para_bdry_x_*.png')


if __name__ == "__main__":
    main()
