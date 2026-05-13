import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_bz
from bdry_operations import extend_abs_layer, arc_stretch_2d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    flag_printf = 1
    flag_topo_z = 1

    num_pml = 0
    nx1 = 1000
    nx = nx1 + 2 * num_pml

    ny1 = 1000
    ny = ny1 + 2 * num_pml

    dx = 1.0
    dy = 1.0
    origin_x = 0.0
    origin_y = 0.0
    origin_z = 0.0

    # bz shape is (ny, nx, 3) as per MATLAB's creat_3d_hyper_bdry_z
    bz = np.zeros((ny, nx, 3))

    # Fill the main region (with PML offset)
    for j in range(ny1):
        for i in range(nx1):
            bz[j + num_pml, i + num_pml, 0] = origin_x + i * dx  # i starts at 0 in Python
            bz[j + num_pml, i + num_pml, 1] = origin_y + j * dy  # j starts at 0 in Python
            bz[j + num_pml, i + num_pml, 2] = origin_z

    if flag_topo_z:
        point_x = origin_x + (nx1 // 2) * dx
        point_y = origin_y + (ny1 // 2) * dy
        L = 0.2 * nx * dx
        H = 0.25 * nx * dx
        for j in range(ny1):
            for i in range(nx1):
                jj = j + num_pml
                ii = i + num_pml
                r1 = np.sqrt((bz[jj, ii, 0] - point_x)**2 + (bz[jj, ii, 1] - point_y)**2)
                topo = 0.0
                if r1 < L:
                    topo = 0.5 * H * (1 + np.cos(np.pi * r1 / L))
                bz[jj, ii, 2] = bz[jj, ii, 2] - topo

    # Apply PML extension
    bz = extend_abs_layer(bz, dx, dy, nx, ny, num_pml)

    # Arc stretch (commented out in MATLAB, so we leave it commented here too)
    # A = 0.00001
    # bz = arc_stretch_2d(A, bz)

    file_name = '../data_file_3d.txt'
    export_bdry_1_bz(bz, nx, ny, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nx, ny) // 80)
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        X = bz[:, :, 0] / 1e3; Y = bz[:, :, 1] / 1e3; Z = bz[:, :, 2] / 1e3
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
        plt.savefig('../plot/hyper_bdry_z.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/hyper_bdry_z.png')

if __name__ == "__main__":
    main()
