"""Export functions for C-project 3D boundary file format.

C format includes integer count lines on the first non-comment line,
which io_get_nextline() reads.
"""
import numpy as np


def export_bdry_6(bx1: np.ndarray, bx2: np.ndarray,
                  by1: np.ndarray, by2: np.ndarray,
                  bz1: np.ndarray, bz2: np.ndarray,
                  nx: int, ny: int, nz: int, file_dir: str) -> None:
    """Export 6-boundary format (elliptic).

    bx1, bx2: shape (ny, nz, 3) — x-direction boundaries
    by1, by2: shape (nx, nz, 3) — y-direction boundaries
    bz1, bz2: shape (nx, ny, 3) — z-direction boundaries
    """
    with open(file_dir, 'w') as fd:
        fd.write("# nx number\n")
        fd.write(f"{nx}\n")
        fd.write("# ny number\n")
        fd.write(f"{ny}\n")
        fd.write("# nz number\n")
        fd.write(f"{nz}\n")

        fd.write("# bx1 coords\n")
        for k in range(nz):
            for j in range(ny):
                fd.write(f'{bx1[j,k,0]:.9e} {bx1[j,k,1]:.9e} {bx1[j,k,2]:.9e}\n')

        fd.write("# bx2 coords\n")
        for k in range(nz):
            for j in range(ny):
                fd.write(f'{bx2[j,k,0]:.9e} {bx2[j,k,1]:.9e} {bx2[j,k,2]:.9e}\n')

        fd.write("# by1 coords\n")
        for k in range(nz):
            for i in range(nx):
                fd.write(f'{by1[i,k,0]:.9e} {by1[i,k,1]:.9e} {by1[i,k,2]:.9e}\n')

        fd.write("# by2 coords\n")
        for k in range(nz):
            for i in range(nx):
                fd.write(f'{by2[i,k,0]:.9e} {by2[i,k,1]:.9e} {by2[i,k,2]:.9e}\n')

        fd.write("# bz1 coords\n")
        for j in range(ny):
            for i in range(nx):
                fd.write(f'{bz1[i,j,0]:.9e} {bz1[i,j,1]:.9e} {bz1[i,j,2]:.9e}\n')

        fd.write("# bz2 coords\n")
        for j in range(ny):
            for i in range(nx):
                fd.write(f'{bz2[i,j,0]:.9e} {bz2[i,j,1]:.9e} {bz2[i,j,2]:.9e}\n')


def export_bdry_2(bz1: np.ndarray, bz2: np.ndarray,
                  nx: int, ny: int, file_dir: str) -> None:
    """Export 2-boundary format (parabolic).

    bz1, bz2: shape (nx, ny, 3)
    """
    with open(file_dir, 'w') as fd:
        fd.write("# nx ny number\n")
        fd.write(f"{nx} {ny}\n")

        fd.write("# bz1 coords\n")
        for j in range(ny):
            for i in range(nx):
                fd.write(f'{bz1[i,j,0]:.9e} {bz1[i,j,1]:.9e} {bz1[i,j,2]:.9e}\n')

        fd.write("# bz2 coords\n")
        for j in range(ny):
            for i in range(nx):
                fd.write(f'{bz2[i,j,0]:.9e} {bz2[i,j,1]:.9e} {bz2[i,j,2]:.9e}\n')


def export_bdry_1_bz(bz: np.ndarray, nx: int, ny: int,
                     file_dir: str) -> None:
    """Export 1-boundary format — z-direction marching (hyperbolic).

    bz: shape (nx, ny, 3) or (ny, nx, 3) depending on convention.
    Written as outer loop j=0:ny-1, inner loop i=0:nx-1.
    """
    with open(file_dir, 'w') as fd:
        fd.write("# nx number\n")
        fd.write(f"{nx}\n")
        fd.write("# ny number\n")
        fd.write(f"{ny}\n")
        fd.write("# bz coords\n")
        for j in range(ny):
            for i in range(nx):
                fd.write(f'{bz[j,i,0]:.9e} {bz[j,i,1]:.9e} {bz[j,i,2]:.9e}\n')


def export_bdry_1_bx(bx: np.ndarray, ny: int, nz: int,
                     file_dir: str) -> None:
    """Export 1-boundary format — x-direction marching (hyperbolic).

    bx: shape (nz, ny, 3).
    Written as outer loop k=0:nz-1, inner loop j=0:ny-1.
    """
    with open(file_dir, 'w') as fd:
        fd.write("# ny number\n")
        fd.write(f"{ny}\n")
        fd.write("# nz number\n")
        fd.write(f"{nz}\n")
        fd.write("# bx coords\n")
        for k in range(nz):
            for j in range(ny):
                fd.write(f'{bx[k,j,0]:.9e} {bx[k,j,1]:.9e} {bx[k,j,2]:.9e}\n')


def export_bdry_1_by(by: np.ndarray, nx: int, nz: int,
                     file_dir: str) -> None:
    """Export 1-boundary format — y-direction marching (hyperbolic).

    by: shape (nz, nx, 3).
    Written as outer loop k=0:nz-1, inner loop i=0:nx-1.
    """
    with open(file_dir, 'w') as fd:
        fd.write("# nx number\n")
        fd.write(f"{nx}\n")
        fd.write("# nz number\n")
        fd.write(f"{nz}\n")
        fd.write("# by coords\n")
        for k in range(nz):
            for i in range(nx):
                fd.write(f'{by[k,i,0]:.9e} {by[k,i,1]:.9e} {by[k,i,2]:.9e}\n')


def export_step(step_vals: np.ndarray, num_of_step: int,
                file_dir: str) -> None:
    """Export step file."""
    with open(file_dir, 'w') as fd:
        fd.write("# number of step\n")
        fd.write(f"{num_of_step}\n")
        fd.write("# step\n")
        for i in range(num_of_step):
            fd.write(f'{step_vals[i]:.9e}\n')


def export_arc_len(arc_vals: np.ndarray, npoints: int,
                   file_dir: str) -> None:
    """Export arc length file (grid post-process stretch)."""
    with open(file_dir, 'w') as fd:
        fd.write("# number of points\n")
        fd.write(f"{npoints}\n")
        fd.write("# arc_len\n")
        for i in range(npoints):
            fd.write(f'{arc_vals[i]:.9e}\n')
