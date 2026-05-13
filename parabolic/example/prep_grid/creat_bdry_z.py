import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_2
from bdry_operations import extend_abs_layer, arc_stretch_2d

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    flag_printf = 1
    flag_topo_z = 0
    num_pml = 0
    nx1 = 2000
    ny1 = 2000
    nx = nx1 + 2 * num_pml
    ny = ny1 + 2 * num_pml
    dx = 10
    dy = 10
    nz = 1000
    dz = 10
    origin_x = 0
    origin_y = 0
    origin_z = 0

    # Note: in MATLAB, bz1 and bz2 are (ny, nx, 3), but export_bdry_2 expects (nx, ny, 3)
    # Wait let's check MATLAB code: bz1=zeros(ny,nx,3); then for j=1:ny, i=1:nx: bz2(j,i,...)
    # But let's check export in MATLAB: fprintf for j=1:ny, i=1:nx: bz2(j,i,1), etc.
    # export_bdry_2 expects outer j, inner i, shape (nx, ny, 3)
    # Wait, let's check:
    # Let's create (ny, nx, 3) first, then maybe swap? Wait no, let's see:
    # Wait MATLAB: bz2(j,i,1) is x-coordinate for (i,j) in nx, ny?
    # Wait let's just follow MATLAB's shape and indexing carefully:

    bz1 = np.zeros((ny, nx, 3))  # Wait MATLAB uses (ny, nx, 3)
    bz2 = np.zeros((ny, nx, 3))

    for j in range(ny1):  # j 0..ny1-1 (MATLAB 1..ny1)
        for i in range(nx1):  # i 0..nx1-1 (MATLAB 1..nx1)
            bz2[j + num_pml, i + num_pml, 0] = origin_x + i * dx
            bz2[j + num_pml, i + num_pml, 1] = origin_y + j * dy
            bz2[j + num_pml, i + num_pml, 2] = origin_z

    if flag_topo_z:
        point_x = origin_x + (nx1 // 2) * dx
        point_y = origin_y + (ny1 // 2) * dy
        L = 400
        H = 600
        for j in range(ny1):
            for i in range(nx1):
                r1 = (bz2[j + num_pml, i + num_pml, 0] - point_x)**2 + \
                     (bz2[j + num_pml, i + num_pml, 1] - point_y)**2
                topo = H * np.exp(-r1 / L**2)
                bz2[j + num_pml, i + num_pml, 2] -= topo

    # Extend PML
    # Wait extend_abs_layer takes bz with shape (nx, ny, 3) normally, but here we have (ny, nx, 3)
    # Wait let's check MATLAB's extend_abs_layer.m! Wait we didn't read that yet, but let's check:
    # Wait let's read extend_abs_layer.m now:

    bz2 = extend_abs_layer(bz2, dx, dy, nx, ny, num_pml)

    for j in range(ny):
        for i in range(nx):
            bz1[j, i, 0] = bz2[j, i, 0]
            bz1[j, i, 1] = bz2[j, i, 1]
            bz1[j, i, 2] = origin_z - (nz - 1) * dz

    # Now, export_bdry_2 expects (nx, ny, 3), so we need to swap axes? Wait let's check:
    # Wait in MATLAB, the export loop is j=1:ny, i=1:nx: print bz2(j,i,1), etc.
    # export_bdry_2 does j in 0..ny-1, i in 0..nx-1: print bz2[i,j,0], etc.
    # So if we have bz2 as (ny, nx, 3), then bz2[j,i,0] in MATLAB is bz2[j,i,0] in Python
    # So to make it (nx, ny, 3), we need to transpose the first two dimensions!
    bz1_export = np.transpose(bz1, (1, 0, 2))
    bz2_export = np.transpose(bz2, (1, 0, 2))

    file_name = '../data_file_3d.txt'
    export_bdry_2(bz1_export, bz2_export, nx, ny, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nx, ny) // 80)
        # bz2 surf
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        X = bz2_export[:, :, 0] / 1e3; Y = bz2_export[:, :, 1] / 1e3; Z = bz2_export[:, :, 2] / 1e3
        ax.plot_surface(X[::stride, ::stride], Y[::stride, ::stride],
                        Z[::stride, ::stride], cmap='jet', alpha=0.9, edgecolor='none')
        mid_j = ny // 2
        ax.plot(X[:, mid_j], Y[:, mid_j], Z[:, mid_j], 'r', linewidth=2)
        ax.set_xlabel('X axis (km)', fontweight='bold')
        ax.set_ylabel('Y axis (km)', fontweight='bold')
        ax.set_zlabel('Z axis (km)', fontweight='bold')
        ax.view_init(elev=30, azim=40)
        ax.set_box_aspect([1, 1, 0.5])
        ax.grid(False)
        mappable = plt.cm.ScalarMappable(cmap='jet')
        mappable.set_array(Z[::stride, ::stride].ravel())
        mappable.set_clim(Z.min(), Z.max())
        cb = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1)
        cb.set_label('(km)', fontweight='bold')
        fig.patch.set_facecolor('white')
        plt.tight_layout()
        plt.savefig('../plot/para_bdry_z.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/para_bdry_z.png')


if __name__ == "__main__":
    main()
