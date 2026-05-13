import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'plotting'))
from export_c import export_bdry_1_bx
from bdry_operations import arc_stretch_2d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    flag_printf = 1
    flag_topo_x = 1

    ny = 400
    nz = 300

    dy = 10.0
    dz = 10.0
    origin_x = 0.0
    origin_y = 0.0
    origin_z = 0.0

    # bx shape is (nz, ny, 3) as per MATLAB
    bx = np.zeros((nz, ny, 3))

    for k in range(nz):
        for j in range(ny):
            bx[k, j, 0] = origin_x
            bx[k, j, 1] = origin_y + j * dy
            bx[k, j, 2] = origin_z + k * dz

    if flag_topo_x:
        point_y = origin_y + (ny // 2) * dy
        point_z = origin_z + (nz // 2) * dz
        L = 0.3 * ny * dy
        H = 0.2 * ny * dy
        for k in range(nz):
            for j in range(ny):
                r1 = np.sqrt((bx[k, j, 1] - point_y)**2 + (bx[k, j, 2] - point_z)**2)
                topo = 0.0
                if r1 < L:
                    topo = 0.5 * H * (1 + np.cos(np.pi * r1 / L))
                bx[k, j, 0] = bx[k, j, 0] + topo

    A = 0.00001
    bx = arc_stretch_2d(A, bx)

    file_name = '../data_file_3d.txt'
    export_bdry_1_bx(bx, ny, nz, file_name)

    if flag_printf:
        os.makedirs('../plot', exist_ok=True)
        stride = max(1, max(nz, ny) // 40)
        # MATLAB: plot3 with permuted axes + original
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        for k in range(0, nz, stride):
            ax.plot(bx[k, :, 0], bx[k, :, 1], bx[k, :, 2], 'b-', linewidth=0.5)
        for j in range(0, ny, stride):
            ax.plot(bx[:, j, 0], bx[:, j, 1], bx[:, j, 2], 'b-', linewidth=0.5)
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title('bx (both directions)')
        plt.tight_layout()
        plt.savefig('../plot/hyper_bdry_x.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print('Saved: ../plot/hyper_bdry_x.png')

if __name__ == "__main__":
    main()
