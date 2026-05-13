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

    nx = 401
    ny = 301
    dh = 10
    origin_x = 0
    origin_y = 0
    origin_z = 0

    bz1 = np.zeros((nx, ny, 3))
    bz2 = np.zeros((nx, ny, 3))

    # Convert to 0-based indices
    p1 = 0
    p2 = 100
    p3 = 150
    p4 = 250
    p5 = 300
    p6 = 400

    for j in range(ny):
        # i from p1 to p2 (inclusive)
        for i in range(p1, p2 + 1):
            bz2[i, j, 0] = origin_x + i * dh
            bz2[i, j, 1] = origin_y + j * dh
            bz2[i, j, 2] = origin_z
        # i from p2+1 to p3
        for i in range(p2 + 1, p3 + 1):
            bz2[i, j, 0] = origin_x + bz2[p2, j, 0] + np.cos(0.5 * np.pi * 70 / 90) * (i - p2) * dh
            bz2[i, j, 1] = origin_y + j * dh
            bz2[i, j, 2] = origin_z + bz2[p2, j, 2] - np.sin(0.5 * np.pi * 70 / 90) * (i - p2) * dh
        # i from p3+1 to p4
        for i in range(p3 + 1, p4 + 1):
            bz2[i, j, 0] = origin_x + bz2[p3, j, 0] + (i - p3) * dh
            bz2[i, j, 1] = origin_y + j * dh
            bz2[i, j, 2] = origin_z + bz2[p3, j, 2]
        # i from p4+1 to p5
        for i in range(p4 + 1, p5 + 1):
            bz2[i, j, 0] = origin_x + bz2[p4, j, 0] + np.cos(0.5 * np.pi * 70 / 90) * (i - p4) * dh
            bz2[i, j, 1] = origin_y + j * dh
            bz2[i, j, 2] = origin_z + bz2[p4, j, 2] + np.sin(0.5 * np.pi * 70 / 90) * (i - p4) * dh
        # i from p5+1 to p6
        for i in range(p5 + 1, p6 + 1):
            bz2[i, j, 0] = origin_x + bz2[p5, j, 0] + (i - p5) * dh
            bz2[i, j, 1] = origin_y + j * dh
            bz2[i, j, 2] = origin_z + bz2[p5, j, 2]

    # Smooth along two directions
    window_len = 40
    for i in range(nx):
        bz2[i, :, 0] = smooth_1d(bz2[i, :, 0], window_len)
        bz2[i, :, 1] = smooth_1d(bz2[i, :, 1], window_len)
        bz2[i, :, 2] = smooth_1d(bz2[i, :, 2], window_len)
    for j in range(ny):
        bz2[:, j, 0] = smooth_1d(bz2[:, j, 0], window_len)
        bz2[:, j, 1] = smooth_1d(bz2[:, j, 1], window_len)
        bz2[:, j, 2] = smooth_1d(bz2[:, j, 2], window_len)

    # bz1 x: linspace(-500, 3800, 401), y: linspace(0, 3000, 301), z=-2000
    x = np.linspace(-500, 3800, nx)
    y = np.linspace(0, 3000, ny)
    for j in range(ny):
        for i in range(nx):
            bz1[i, j, 0] = x[i]
            bz1[i, j, 1] = y[j]
            bz1[i, j, 2] = -2000

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
        plt.savefig('../plot/para_model2.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/para_model2.png')


if __name__ == "__main__":
    main()
