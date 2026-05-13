import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_bz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bdry_operations import arc_stretch_2d

def smooth1d(arr, span):
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

    dx = 10.0
    dy = 10.0

    origin_x = 0.0
    origin_y = 0.0
    origin_z = 0.0

    # MATLAB's bz is (nx, ny, 3), but we need (ny, nx, 3) for export_bdry_1_bz
    # Wait, let's think: MATLAB writes j (1:ny) then i (1:nx), and for each, writes bz(i,j,:)
    # So in Python, to have bz[j,i,:] correspond to MATLAB's bz(i,j,:), we need to make
    # bz shape (ny, nx, 3), and fill j (0:ny-1), i (0:nx-1)
    bz = np.zeros((ny, nx, 3))

    for j in range(ny):
        for i in range(nx):
            bz[j, i, 0] = origin_x + i * dx
            bz[j, i, 1] = origin_y + j * dy
            bz[j, i, 2] = origin_z

    if flag_topo_z:
        point_x = origin_x + (nx // 2) * dx
        point_y = origin_y + (ny // 2) * dy
        L = 200.0
        H = 2000.0
        for j in range(ny):
            for i in range(nx):
                r1 = (bz[j, i, 0] - point_x)**2 + (bz[j, i, 1] - point_y)**2
                topo = H * np.exp(-r1 / (L**2))
                if topo > 500.0:
                    topo = 500.0
                bz[j, i, 2] = bz[j, i, 2] - topo

    # Smooth along both directions (replicate MATLAB's smooth)
    # First, smooth along i (second dimension) for each j
    window_len = 20
    for j in range(ny):
        bz[j, :, 0] = smooth1d(bz[j, :, 0], window_len)
        bz[j, :, 1] = smooth1d(bz[j, :, 1], window_len)
        bz[j, :, 2] = smooth1d(bz[j, :, 2], window_len)

    # Then smooth along j (first dimension) for each i
    for i in range(nx):
        bz[:, i, 0] = smooth1d(bz[:, i, 0], window_len)
        bz[:, i, 1] = smooth1d(bz[:, i, 1], window_len)
        bz[:, i, 2] = smooth1d(bz[:, i, 2], window_len)

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
        # MATLAB: plot3(bz(:,201,1)/1e3,...)
        mid_j = 201
        ax.plot(X[:, mid_j], Y[:, mid_j], Z[:, mid_j], 'r', linewidth=3)
        ax.set_xlabel('X axis (km)', fontweight='bold')
        ax.set_ylabel('Y axis (km)', fontweight='bold')
        ax.set_zlabel('Z axis (km)', fontweight='bold')
        ax.view_init(elev=70, azim=40)
        ax.set_box_aspect([1, 1, 0.5])
        ax.grid(False)
        mappable = plt.cm.ScalarMappable(cmap='jet')
        mappable.set_array(Z[::stride, ::stride].ravel())
        mappable.set_clim(Z.min(), Z.max())
        cb = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1)
        cb.set_label('(km)', fontweight='bold')
        fig.patch.set_facecolor('white')
        plt.tight_layout()
        plt.savefig('../plot/hyper_model3.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/hyper_model3.png')

if __name__ == "__main__":
    main()
