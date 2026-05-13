import numpy as np


def extend_abs_layer(bz: np.ndarray, dx: float, dy: float,
                     nx: int, ny: int, num_pml: int) -> np.ndarray:
    """Extend PML absorbing layer on a 3D boundary surface.

    bz: shape (nx, ny, 3) for z-direction boundaries (bz1, bz2)
        or shape (ny, nx, 3) as used in hyperbolic creat_3d_hyper_bdry_z
    """
    coef = 0.7

    # x1 left
    for i in range(num_pml - 1, -1, -1):
        for j in range(ny):
            dz = bz[i + 2, j, 2] - bz[i + 1, j, 2]
            slope = coef * dz / dx
            bz[i, j, 0] = bz[i + 1, j, 0] - dx
            bz[i, j, 1] = bz[i + 1, j, 1]
            bz[i, j, 2] = bz[i + 1, j, 2] - slope * dx

    # x2 right
    for i in range(nx - num_pml, nx):
        for j in range(ny):
            dz = bz[i - 1, j, 2] - bz[i - 2, j, 2]
            slope = coef * dz / dx
            bz[i, j, 0] = bz[i - 1, j, 0] + dx
            bz[i, j, 1] = bz[i - 1, j, 1]
            bz[i, j, 2] = bz[i - 1, j, 2] + slope * dx

    # y1 front
    for i in range(nx):
        for j in range(num_pml - 1, -1, -1):
            dz = bz[i, j + 2, 2] - bz[i, j + 1, 2]
            slope = coef * dz / dy
            bz[i, j, 0] = bz[i, j + 1, 0]
            bz[i, j, 1] = bz[i, j + 1, 1] - dy
            bz[i, j, 2] = bz[i, j + 1, 2] - slope * dy

    # y2 back
    for i in range(nx):
        for j in range(ny - num_pml, ny):
            dz = bz[i, j - 1, 2] - bz[i, j - 2, 2]
            slope = coef * dz / dy
            bz[i, j, 0] = bz[i, j - 1, 0]
            bz[i, j, 1] = bz[i, j - 1, 1] + dy
            bz[i, j, 2] = bz[i, j - 1, 2] + slope * dy

    return bz


def arc_stretch_1d(A: float, bdry: np.ndarray) -> np.ndarray:
    """Arc-length redistribution along the second dimension of a 3D boundary.

    bdry: shape (n1, n2, 3) — redistributed along n2 axis.
    """
    n2 = bdry.shape[1]
    n1 = bdry.shape[0]

    for k in range(n1):
        x_0 = bdry[k, :, 0].copy()
        y_0 = bdry[k, :, 1].copy()
        z_0 = bdry[k, :, 2].copy()

        # arc length along second dimension
        s = np.zeros(n2)
        for i in range(1, n2):
            s[i] = s[i - 1] + np.sqrt(
                (x_0[i] - x_0[i - 1]) ** 2 +
                (y_0[i] - y_0[i - 1]) ** 2 +
                (z_0[i] - z_0[i - 1]) ** 2)

        # normalized
        u = np.zeros(n2)
        for i in range(1, n2):
            u[i] = s[i] / s[-1]

        for i in range(1, n2 - 1):
            xi = i / (n2 - 1)
            r = (np.exp(A * xi) - 1) / (np.exp(A) - 1)

            n = 0
            for m in range(n2 - 1):
                if r >= u[m] and r < u[m + 1]:
                    n = m
                    break

            bdry[k, i, 0] = x_0[n] + (x_0[n + 1] - x_0[n]) * (r - u[n]) / (u[n + 1] - u[n])
            bdry[k, i, 1] = y_0[n] + (y_0[n + 1] - y_0[n]) * (r - u[n]) / (u[n + 1] - u[n])
            bdry[k, i, 2] = z_0[n] + (z_0[n + 1] - z_0[n]) * (r - u[n]) / (u[n + 1] - u[n])

    return bdry


def arc_stretch_2d(A: float, bdry: np.ndarray) -> np.ndarray:
    """Arc-length redistribution along both dimensions of a 3D boundary.

    bdry: shape (n1, n2, 3) — first along n2, then along n1.
    """
    n1 = bdry.shape[0]
    n2 = bdry.shape[1]

    x_0 = bdry[:, :, 0].copy()
    y_0 = bdry[:, :, 1].copy()
    z_0 = bdry[:, :, 2].copy()
    x = bdry[:, :, 0].copy()
    y = bdry[:, :, 1].copy()
    z = bdry[:, :, 2].copy()

    # along second dimension
    s = np.zeros((n1, n2))
    for j in range(n1):
        for i in range(1, n2):
            s[j, i] = s[j, i - 1] + np.sqrt(
                (x_0[j, i] - x_0[j, i - 1]) ** 2 +
                (y_0[j, i] - y_0[j, i - 1]) ** 2 +
                (z_0[j, i] - z_0[j, i - 1]) ** 2)

    u = np.zeros((n1, n2))
    for j in range(n1):
        for i in range(1, n2):
            u[j, i] = s[j, i] / s[j, -1]

    for j in range(n1):
        for i in range(1, n2 - 1):
            xi = i / (n2 - 1)
            r = (np.exp(A * xi) - 1) / (np.exp(A) - 1)
            n = 0
            for m in range(n2 - 1):
                if r >= u[j, m] and r < u[j, m + 1]:
                    n = m
                    break
            bdry[j, i, 0] = x_0[j, n] + (x_0[j, n + 1] - x_0[j, n]) * (r - u[j, n]) / (u[j, n + 1] - u[j, n])
            bdry[j, i, 1] = y_0[j, n] + (y_0[j, n + 1] - y_0[j, n]) * (r - u[j, n]) / (u[j, n + 1] - u[j, n])
            bdry[j, i, 2] = z_0[j, n] + (z_0[j, n + 1] - z_0[j, n]) * (r - u[j, n]) / (u[j, n + 1] - u[j, n])

    # update x, y, z for second pass
    x[:] = bdry[:, :, 0]
    y[:] = bdry[:, :, 1]
    z[:] = bdry[:, :, 2]

    # along first dimension
    s = np.zeros((n1, n2))
    for j in range(1, n1):
        for i in range(n2):
            s[j, i] = s[j - 1, i] + np.sqrt(
                (x[j, i] - x[j - 1, i]) ** 2 +
                (y[j, i] - y[j - 1, i]) ** 2 +
                (z[j, i] - z[j - 1, i]) ** 2)

    u = np.zeros((n1, n2))
    for j in range(1, n1):
        for i in range(n2):
            u[j, i] = s[j, i] / s[-1, i]

    for j in range(1, n1 - 1):
        for i in range(n2):
            et = j / (n1 - 1)
            r = (np.exp(A * et) - 1) / (np.exp(A) - 1)
            n = 0
            for m in range(n1 - 1):
                if r >= u[m, i] and r < u[m + 1, i]:
                    n = m
                    break
            bdry[j, i, 0] = x[n, i] + (x[n + 1, i] - x[n, i]) * (r - u[n, i]) / (u[n + 1, i] - u[n, i])
            bdry[j, i, 1] = y[n, i] + (y[n + 1, i] - y[n, i]) * (r - u[n, i]) / (u[n + 1, i] - u[n, i])
            bdry[j, i, 2] = z[n, i] + (z[n + 1, i] - z[n, i]) * (r - u[n, i]) / (u[n + 1, i] - u[n, i])

    return bdry
