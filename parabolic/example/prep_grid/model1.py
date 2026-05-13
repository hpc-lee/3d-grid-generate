import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_2

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def smooth_1d(arr, span):
    """Reimplementation of MATLAB smooth(y, span) with 'moving' method.

    Even span is reduced by 1 (MATLAB convention). At boundaries, uses
    anchored windows (left-anchored at start, right-anchored at end)
    rather than centered truncated windows.
    """
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    if span % 2 == 0:
        span = span - 1
    half = span // 2
    result = np.zeros(n)
    for i in range(n):
        if i < half:
            lo = 0
            hi = 2 * i
        elif i >= n - half:
            lo = 2 * i - n + 1
            hi = n - 1
        else:
            lo = i - half
            hi = i + half
        result[i] = np.mean(arr[lo:hi+1])
    return result


def main():
    flag_printf = 1
    flag_topo_z = 1

    nx = 401
    ny = 401
    dx = 10
    dy = 10

    origin_x = 0
    origin_y = 0
    origin_z = 0

    bz1 = np.zeros((nx, ny, 3))
    bz2 = np.zeros((nx, ny, 3))

    for j in range(ny):
        for i in range(nx):
            bz2[i, j, 0] = origin_x + i * dx
            bz2[i, j, 1] = origin_y + j * dy
            bz2[i, j, 2] = origin_z

    if flag_topo_z:
        x1 = 1000
        y1 = 2000
        x2 = 3000
        y2 = 2000
        L = 200
        H = 500
        for j in range(ny):
            for i in range(nx):
                r1 = (bz2[i, j, 0] - x1)**2 + (bz2[i, j, 1] - y1)**2
                r2 = (bz2[i, j, 0] - x2)**2 + (bz2[i, j, 1] - y2)**2
                topo = H * np.exp(-r1 / L**2) - H * np.exp(-r2 / L**2)
                bz2[i, j, 2] += topo

    # Smooth along two directions
    window_len = 20
    for i in range(nx):
        bz2[i, :, 0] = smooth_1d(bz2[i, :, 0], window_len)
        bz2[i, :, 1] = smooth_1d(bz2[i, :, 1], window_len)
        bz2[i, :, 2] = smooth_1d(bz2[i, :, 2], window_len)
    for j in range(ny):
        bz2[:, j, 0] = smooth_1d(bz2[:, j, 0], window_len)
        bz2[:, j, 1] = smooth_1d(bz2[:, j, 1], window_len)
        bz2[:, j, 2] = smooth_1d(bz2[:, j, 2], window_len)

    # bz1 gets x,y from bz2, z=-2000
    bz1[:, :, 0] = bz2[:, :, 0]
    bz1[:, :, 1] = bz2[:, :, 1]
    bz1[:, :, 2] = -2000

    file_name = '../data_file_3d.txt'
    export_bdry_2(bz1, bz2, nx, ny, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nx, ny) // 80)
        fig = plt.figure(figsize=(11, 6))
        ax = fig.add_subplot(111, projection='3d')
        X = bz2[:, :, 0] / 1e3; Y = bz2[:, :, 1] / 1e3; Z = bz2[:, :, 2] / 1e3
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
        plt.savefig('../plot/para_model1.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/para_model1.png')


if __name__ == "__main__":
    main()
