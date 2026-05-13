import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_6
from bdry_operations import extend_abs_layer, arc_stretch_1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    # Parameters
    flag_printf = 1
    flag_topo_x = 0
    flag_topo_y = 0
    flag_topo_z = 0

    num_pml = 0
    nx1 = 101
    ny1 = 71

    nx = nx1 + 2 * num_pml
    ny = ny1 + 2 * num_pml
    nz = 51

    dx = 10.0
    dy = 10.0
    dz = 10.0

    origin_x = 0.0
    origin_y = 0.0
    origin_z = 0.0

    # Initialize boundaries
    bz1 = np.zeros((nx, ny, 3))
    bz2 = np.zeros((nx, ny, 3))
    by1 = np.zeros((nx, nz, 3))
    by2 = np.zeros((nx, nz, 3))
    bx1 = np.zeros((ny, nz, 3))
    bx2 = np.zeros((ny, nz, 3))

    # Fill z boundaries (bz1 bottom, bz2 top)
    for j in range(ny1):
        for i in range(nx1):
            ii = i + num_pml
            jj = j + num_pml
            bz1[ii, jj, 0] = origin_x + i * dx
            bz1[ii, jj, 1] = origin_y + j * dy
            bz1[ii, jj, 2] = origin_z - (nz - 1) * dz

            bz2[ii, jj, 0] = origin_x + i * dx
            bz2[ii, jj, 1] = origin_y + j * dy
            bz2[ii, jj, 2] = origin_z

    # Add Gaussian topography if flag_topo_z is set
    if flag_topo_z:
        point_x = origin_x + (nx1 // 2) * dx
        point_y = origin_y + (ny1 // 2) * dy
        L = 0.3 * nx * dx
        H = 0.1 * nx * dx
        for j in range(ny1):
            for i in range(nx1):
                ii = i + num_pml
                jj = j + num_pml
                r1 = np.sqrt(
                    (bz2[ii, jj, 0] - point_x)**2 +
                    (bz2[ii, jj, 1] - point_y)**2
                )
                topo = 0.0
                if r1 < L:
                    topo = 0.5 * H * (1 + np.cos(np.pi * r1 / L))
                bz2[ii, jj, 2] += topo

    # Extend PML layers
    bz1 = extend_abs_layer(bz1, dx, dy, nx, ny, num_pml)
    bz2 = extend_abs_layer(bz2, dx, dy, nx, ny, num_pml)

    # Fill x boundaries (bx1 left, bx2 right)
    for j in range(ny):
        dz1 = (bz2[0, j, 2] - bz1[0, j, 2]) / (nz - 1)
        dz2 = (bz2[nx-1, j, 2] - bz1[nx-1, j, 2]) / (nz - 1)
        for k in range(nz):
            bx1[j, k, 0] = bz1[0, j, 0]
            bx1[j, k, 1] = bz1[0, j, 1]
            bx1[j, k, 2] = bz1[0, j, 2] + k * dz1

            bx2[j, k, 0] = bz1[nx-1, j, 0]
            bx2[j, k, 1] = bz1[nx-1, j, 1]
            bx2[j, k, 2] = bz1[nx-1, j, 2] + k * dz2

    # Fill y boundaries (by1 front, by2 back)
    for i in range(nx):
        dz1 = (bz2[i, 0, 2] - bz1[i, 0, 2]) / (nz - 1)
        dz2 = (bz2[i, ny-1, 2] - bz1[i, ny-1, 2]) / (nz - 1)
        for k in range(nz):
            by1[i, k, 0] = bz1[i, 0, 0]
            by1[i, k, 1] = bz1[i, 0, 1]
            by1[i, k, 2] = bz1[i, 0, 2] + k * dz1

            by2[i, k, 0] = bz1[i, ny-1, 0]
            by2[i, k, 1] = bz1[i, ny-1, 1]
            by2[i, k, 2] = bz1[i, ny-1, 2] + k * dz2

    A = 0.0001
    # Uncomment these if you need arc stretching
    # by1 = arc_stretch_1d(A, by1)
    # by2 = arc_stretch_1d(A, by2)
    # bz1 = arc_stretch_1d(A, bz1)
    # bz2 = arc_stretch_1d(A, bz2)
    #
    # bx1 = arc_stretch_1d(A, bx1)
    # bx2 = arc_stretch_1d(A, bx2)
    # bz1 = arc_stretch_1d(A, bz1.swapaxes(0,1)).swapaxes(0,1)
    # bz2 = arc_stretch_1d(A, bz2.swapaxes(0,1)).swapaxes(0,1)
    #
    # bx1 = arc_stretch_1d(A, bx1.swapaxes(0,1)).swapaxes(0,1)
    # bx2 = arc_stretch_1d(A, bx2.swapaxes(0,1)).swapaxes(0,1)
    # by1 = arc_stretch_1d(A, by1.swapaxes(0,1)).swapaxes(0,1)
    # by2 = arc_stretch_1d(A, by2.swapaxes(0,1)).swapaxes(0,1)

    # Export boundary file
    file_dir = os.path.join(os.path.dirname(__file__), '../data_file_3d.txt')
    export_bdry_6(bx1, bx2, by1, by2, bz1, bz2, nx, ny, nz, file_dir)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nx, ny, nz) // 30)

        # Figure 1: bx1 + bx2 (plot3 style, like MATLAB)
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        for k in range(0, bx1.shape[1], stride):
            ax.plot(bx1[:, k, 0], bx1[:, k, 1], bx1[:, k, 2], 'b-', linewidth=0.5)
            ax.plot(bx2[:, k, 0], bx2[:, k, 1], bx2[:, k, 2], 'b-', linewidth=0.5)
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('bx1 + bx2')
        ax.set_box_aspect([1, 1, 0.5])
        plt.tight_layout()
        plt.savefig('../plot/elli_gauss_bx.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        # Figure 2: by1 + by2
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        for k in range(0, by1.shape[1], stride):
            ax.plot(by1[:, k, 0], by1[:, k, 1], by1[:, k, 2], 'b-', linewidth=0.5)
            ax.plot(by2[:, k, 0], by2[:, k, 1], by2[:, k, 2], 'b-', linewidth=0.5)
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('by1 + by2')
        ax.set_box_aspect([1, 1, 0.5])
        plt.tight_layout()
        plt.savefig('../plot/elli_gauss_by.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        # Figure 3: bz1 + bz2
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        for j in range(0, bz1.shape[1], stride):
            ax.plot(bz1[:, j, 0], bz1[:, j, 1], bz1[:, j, 2], 'b-', linewidth=0.5)
            ax.plot(bz2[:, j, 0], bz2[:, j, 1], bz2[:, j, 2], 'b-', linewidth=0.5)
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('bz1 + bz2')
        ax.set_box_aspect([1, 1, 0.5])
        plt.tight_layout()
        plt.savefig('../plot/elli_gauss_bz.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        # Figure 4: bz2 surf (like MATLAB's surf with 'edgecolor','none')
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        X = bz2[:, :, 0]; Y = bz2[:, :, 1]; Z = bz2[:, :, 2]
        ax.plot_surface(X[::stride, ::stride], Y[::stride, ::stride],
                        Z[::stride, ::stride], cmap='jet', edgecolor='none')
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('bz2 (top surface)')
        ax.set_box_aspect([1, 1, 0.5])
        plt.tight_layout()
        plt.savefig('../plot/elli_gauss_bz2_surf.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        print('Saved: ../plot/elli_gauss_*.png')


if __name__ == "__main__":
    main()
