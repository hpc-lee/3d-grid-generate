import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_bz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    flag_printf = 1

    nx = 401
    ny = 301

    p1 = 0      # 1 in MATLAB, 0-based here
    p2 = 80      # 81 in MATLAB
    p3 = 130     # 131 in MATLAB
    p4 = 270     # 271 in MATLAB
    p5 = 320     # 321 in MATLAB
    p6 = 400     # 401 in MATLAB

    dh = 10.0
    origin_x = 0.0
    origin_y = 0.0
    origin_z = 0.0

    # bz shape is (ny, nx, 3)
    bz = np.zeros((ny, nx, 3))

    for j in range(ny):
        # First segment
        for i in range(p1, p2 + 1):
            bz[j, i, 0] = origin_x + i * dh
            bz[j, i, 1] = origin_y + j * dh
            bz[j, i, 2] = origin_z

        # Second segment
        for i in range(p2 + 1, p3 + 1):
            bz[j, i, 0] = origin_x + bz[j, p2, 0] + np.cos(0.5 * np.pi * 90 / 90) * (i - p2) * dh
            bz[j, i, 1] = origin_y + j * dh
            bz[j, i, 2] = origin_z + bz[j, p2, 2] + np.sin(0.5 * np.pi * 90 / 90) * (i - p2) * dh

        # Third segment
        for i in range(p3 + 1, p4 + 1):
            bz[j, i, 0] = origin_x + bz[j, p3, 0] + (i - p3) * dh
            bz[j, i, 1] = origin_y + j * dh
            bz[j, i, 2] = origin_z + bz[j, p3, 2]

        # Fourth segment
        for i in range(p4 + 1, p5 + 1):
            bz[j, i, 0] = origin_x + bz[j, p4, 0] + np.cos(0.5 * np.pi * 90 / 90) * (i - p4) * dh
            bz[j, i, 1] = origin_y + j * dh
            bz[j, i, 2] = origin_z + bz[j, p4, 2] - np.sin(0.5 * np.pi * 90 / 90) * (i - p4) * dh

        # Fifth segment
        for i in range(p5 + 1, p6 + 1):
            bz[j, i, 0] = origin_x + bz[j, p5, 0] + (i - p5) * dh
            bz[j, i, 1] = origin_y + j * dh
            bz[j, i, 2] = origin_z + bz[j, p5, 2]

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
        # MATLAB: plot3(bz(:,151,1)/1e3,...)
        mid_j = 151
        ax.plot(X[:, mid_j], Y[:, mid_j], Z[:, mid_j], 'r', linewidth=3)
        ax.set_xlabel('X axis (km)', fontweight='bold')
        ax.set_ylabel('Y axis (km)', fontweight='bold')
        ax.set_zlabel('Z axis (km)', fontweight='bold')
        ax.view_init(elev=30, azim=20)
        ax.set_box_aspect([1, 1, 0.5])
        ax.grid(False)
        mappable = plt.cm.ScalarMappable(cmap='jet')
        mappable.set_array(Z[::stride, ::stride].ravel())
        mappable.set_clim(Z.min(), Z.max())
        cb = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1)
        cb.set_label('(km)', fontweight='bold')
        fig.patch.set_facecolor('white')
        plt.tight_layout()
        plt.savefig('../plot/hyper_model2.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/hyper_model2.png')

if __name__ == "__main__":
    main()
