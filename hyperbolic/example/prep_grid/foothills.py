import sys, os
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import griddata
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_bz
from bdry_operations import extend_abs_layer
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    flag_printf = 1
    flag_sample = 0

    # Load topo data
    topo_data = loadmat("/data/lihl/code/foothills_3d/topo_new.mat")
    # Assume the variable name in the mat file is 'topo' (common for importdata)
    # Wait, let's check: MATLAB's importdata stores in .data if it's a matrix
    # But let's just get the first variable if possible
    topo = None
    for key in topo_data:
        if not key.startswith('__'):
            topo = topo_data[key]
            break
    if topo is None:
        raise ValueError("Could not find topo data in mat file")

    num_pml = 20
    nx1 = 501
    ny1 = 501
    nx = nx1 + 2 * num_pml
    ny = ny1 + 2 * num_pml

    i1 = 0          # 1 in MATLAB, 0-based here
    i2 = i1 + nx1   # i1 + nx1 in 0-based (since i1:i2 in MATLAB is i1 to i2 inclusive, which is i2 - i1 + 1 elements)
    j1 = 500        # 501 in MATLAB, 0-based here
    j2 = j1 + ny1
    topo_local = topo[i1:i2, j1:j2]

    nz = 301
    dx = 10.0
    dy = 10.0
    dz = 10.0
    org_x = i1 * dx  # (i1 - 1)*dx in MATLAB is i1*dx in Python since i1 is 0-based
    org_y = j1 * dy

    # bz2 shape is (ny, nx, 3) for export_bdry_1_bz
    bz2 = np.zeros((ny, nx, 3))

    for j in range(ny1):
        for i in range(nx1):
            jj = j + num_pml
            ii = i + num_pml
            bz2[jj, ii, 0] = org_x + i * dx
            bz2[jj, ii, 1] = org_y + j * dy
            bz2[jj, ii, 2] = topo_local[i, j]

    # Extend PML (note: extend_abs_layer expects bz shape (ny, nx, 3) here)
    bz2 = extend_abs_layer(bz2, dx, dy, nx, ny, num_pml)

    if flag_sample:
        nx_new = 641
        ny_new = 641
        bz2_new = np.zeros((ny_new, nx_new, 3))
        x_line = np.linspace(bz2[0, 0, 0], bz2[0, nx-1, 0], nx_new)
        y_line = np.linspace(bz2[0, 0, 1], bz2[ny-1, 0, 1], ny_new)
        xx, yy = np.meshgrid(x_line, y_line)
        # Get all points for interpolation
        points = np.vstack((bz2[:, :, 0].ravel(), bz2[:, :, 1].ravel())).T
        values = bz2[:, :, 2].ravel()
        # Interpolate
        bz2_new[:, :, 0] = xx
        bz2_new[:, :, 1] = yy
        bz2_new[:, :, 2] = griddata(points, values, (xx, yy), method='cubic')
    else:
        bz2_new = bz2
        nx_new = nx
        ny_new = ny

    file_name = '../data_file_3d.txt'
    export_bdry_1_bz(bz2_new, nx_new, ny_new, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nx_new, ny_new) // 80)
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        X = bz2_new[:, :, 0] / 1e3; Y = bz2_new[:, :, 1] / 1e3; Z = bz2_new[:, :, 2] / 1e3
        ax.plot_surface(X[::stride, ::stride], Y[::stride, ::stride],
                        Z[::stride, ::stride], cmap='jet', alpha=0.9, edgecolor='none')
        # MATLAB: receiver markers at specific indices + source marker
        receivers = [(121, 401), (201, 401), (281, 401), (361, 401), (441, 401)]
        for idx, (ri, rj) in enumerate(receivers):
            ri2 = min(ri - 1, nx_new - 1)
            rj2 = min(rj - 1, ny_new - 1)
            ax.plot([X[ri2, rj2]], [Y[ri2, rj2]], [Z[ri2, rj2] + 0.1],
                    'rv', markersize=10, markeredgewidth=2)
            ax.text(X[ri2, rj2] - 0.1, Y[ri2, rj2], Z[ri2, rj2] + 0.15,
                    f'R{idx+1}', color='r', fontsize=12, fontweight='bold')
        # Source marker
        si, sj = 221, 221
        si2 = min(si - 1, nx_new - 1)
        sj2 = min(sj - 1, ny_new - 1)
        ax.plot([X[si2, sj2]], [Y[si2, sj2]], [Z[si2, sj2] + 0.25],
                'rp', markersize=12, markeredgewidth=2)
        ax.set_xlabel('x-axis (km)', fontweight='bold', fontsize=15)
        ax.set_ylabel('y-axis (km)', fontweight='bold', fontsize=15)
        ax.set_zlabel('z-axis (km)', fontweight='bold', fontsize=15)
        ax.view_init(elev=30, azim=60)
        ax.set_box_aspect([1, 1, 0.5])
        ax.grid(False)
        mappable = plt.cm.ScalarMappable(cmap='jet')
        mappable.set_array(Z[::stride, ::stride].ravel())
        mappable.set_clim(0.2, 1.5)
        cb = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1)
        cb.set_label('(km)', fontweight='bold')
        fig.patch.set_facecolor('white')
        plt.tight_layout()
        plt.savefig('../plot/foothills.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/foothills.png')

if __name__ == "__main__":
    main()
