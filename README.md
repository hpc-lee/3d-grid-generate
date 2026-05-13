# 3D Curvilinear Grid Generation (C)

C implementation of 3D curvilinear grid generation for seismic wave numerical simulation. Provides three grid generation methods — parabolic, hyperbolic, and elliptic — with grid quality evaluation, post-processing, and Python-based visualization.

## Project Structure

```
3d-grid-generate/
├── lib/                      # Shared C libraries
│   ├── cJSON.c / cJSON.h     # JSON parser
│   ├── lib_math.c / lib_math.h  # 3x3 matrix operations + geometry
│   └── lib_mem.c / lib_mem.h    # Memory allocation utilities
├── plotting/                 # Grid visualization (Python)
│   ├── draw_grid.py          # Plot 3D grid coordinates (2D slices)
│   ├── draw_quality.py       # Plot quality metrics (2D slices)
│   ├── draw_grid_func.py     # NetCDF I/O utilities (3D)
│   ├── export_c.py           # C-format boundary export functions
│   └── bdry_operations.py    # 3D PML extension + arc-length stretching
├── elliptic/                 # Elliptic grid generation (MPI parallel)
│   ├── src/                  # C source code
│   ├── Makefile
│   └── example/
│       ├── grid_generate.sh  # Run script
│       └── prep_grid/        # Python boundary generators
├── hyperbolic/               # Hyperbolic grid generation
│   ├── src/
│   ├── Makefile
│   └── example/
│       ├── model*.sh
│       ├── grid_generation.sh
│       └── prep_grid/
├── parabolic/                # Parabolic grid generation
│   ├── src/
│   ├── Makefile
│   └── example/
│       ├── model*.sh
│       ├── grid_generation.sh
│       └── prep_grid/
├── grid-post-process/        # Grid post-processing (merge, sample, stretch, PML)
│   ├── src/
│   ├── Makefile
│   └── example/
│       ├── grid_post_proc.sh
│       └── creat_arc_len.py
└── docs/                     # Documentation
    ├── user_manual.tex
    └── user_manual.pdf
```

## Features

- **Three grid generation methods** with different algorithmic trade-offs
- **Grid quality evaluation**: orthogonality (3 components), Jacobian, smoothness, step sizes
- **Post-processing**: multi-grid merging, bilinear sampling, arc-length stretching, PML layer support
- **MPI parallel execution** for elliptic and parabolic solvers (3D Cartesian decomposition)
- **Three marching directions** for hyperbolic and parabolic methods (x, y, z)
- **NetCDF I/O** for interoperability with seismic simulation codes
- **Python-based boundary generation** and visualization

## Prerequisites

### C Compilation
- GCC or MPI compiler (mpicc)
- NetCDF library
- Environment variables: `GNUHOME`, `MPIHOME` (for elliptic/parabolic), `NETCDF`

### Python (for prep_grid and plotting)
- Python >= 3.10
- NumPy
- Matplotlib
- netCDF4
- SciPy (for model5 topography import)

## Installation

### Build C Programs

```bash
# Hyperbolic
cd hyperbolic && make && cd ..

# Parabolic (requires MPI)
cd parabolic && make && cd ..

# Grid Post-Process
cd grid-post-process && make && cd ..

# Elliptic (requires MPI)
cd elliptic && make && cd ..
```

## Quick Start

### Step 1: Generate Boundary Data

```bash
cd elliptic/example/prep_grid
python3 creat_3d_gauss_boundary.py   # generates data_file_3d.txt
```

### Step 2: Run Grid Generation

```bash
cd elliptic/example
bash grid_generate.sh
```

### Step 3: Visualize Results

```bash
cd plotting
# Edit cfs_file in draw_grid.py to point to output config.json
python3 draw_grid.py
python3 draw_quality.py
```

## Configuration

All parameters are specified in a JSON file. Keys prefixed with `#` are treated as comments and ignored.

### Common Parameters

| Key | Type | Description |
|-----|------|-------------|
| `number_of_grid_points_x` | int | Number of grid points in x (xi) direction |
| `number_of_grid_points_y` | int | Number of grid points in y (eta) direction |
| `number_of_grid_points_z` | int | Number of grid points in z (zeta) direction |
| `check_orth_xiet` | int (0/1) | Check orthogonality xi-eta |
| `check_orth_xizt` | int (0/1) | Check orthogonality xi-zeta |
| `check_orth_etzt` | int (0/1) | Check orthogonality eta-zeta |
| `check_jac` | int (0/1) | Check Jacobian determinant |
| `check_smooth_xi` | int (0/1) | Check smoothness in xi direction |
| `check_smooth_et` | int (0/1) | Check smoothness in eta direction |
| `check_smooth_zt` | int (0/1) | Check smoothness in zeta direction |
| `check_step_xi` | int (0/1) | Check step size in xi direction |
| `check_step_et` | int (0/1) | Check step size in eta direction |
| `check_step_zt` | int (0/1) | Check step size in zeta direction |
| `geometry_input_file` | string | Path to boundary geometry file |
| `grid_export_dir` | string | Output directory for grid files |

### Parabolic Parameters

| Key | Type | Description |
|-----|------|-------------|
| `step_input_file` | string | Path to step length file |
| `coef` | int | Clustering coefficient |
| `parabolic.direction` | string | Marching direction ("x", "y", "z") |

### Hyperbolic Parameters

| Key | Type | Description |
|-----|------|-------------|
| `step_input_file` | string | Path to step length file |
| `flag_stretch` | int (0/1) | Enable arc-length stretching |
| `hyperbolic.coef` | int | Smoothing/dissipation coefficient |
| `hyperbolic.bdry_type` | [int, int] | Boundary types |
| `hyperbolic.epsilon` | [float, float] | Damping parameters |
| `hyperbolic.direction` | string | Marching direction ("x", "y", "z") |

### Elliptic Parameters

| Key | Type | Description |
|-----|------|-------------|
| `number_of_mpiprocs_x` | int | MPI processes in x direction |
| `number_of_mpiprocs_y` | int | MPI processes in y direction |
| `number_of_mpiprocs_z` | int | MPI processes in z direction |
| `grid_method.elli_diri.coef` | [int x4] | Source term coefficients |
| `grid_method.elli_diri.weight` | [float x2] | Exponential decay weights |
| `grid_method.elli_diri.iter_err` | float | Convergence tolerance |
| `grid_method.elli_diri.max_iter` | int | Maximum SOR iterations |

### Post-Processing Parameters

| Key | Type | Description |
|-----|------|-------------|
| `input_grid_number` | int | Number of input grids to merge |
| `merge_direction` | string | "x", "y", or "z" (required if input_grid_number >= 2) |
| `flag_sample` | int (0/1) | Enable bilinear sampling |
| `sample_factor` | [int, int, int] | Upsampling factors [x, y, z] |
| `pml_weight_2x` | int (0/1) | Double weight for PML layers |

## Grid Quality Metrics

| Metric | Variable | Description | Ideal Value |
|--------|----------|-------------|-------------|
| Orthogonality xi-eta | `orth_xiet` | Angle deviation between xi and eta grid lines | 90 deg |
| Orthogonality xi-zeta | `orth_xizt` | Angle deviation between xi and zeta grid lines | 90 deg |
| Orthogonality eta-zeta | `orth_etzt` | Angle deviation between eta and zeta grid lines | 90 deg |
| Jacobian | `jacobi` | Determinant of coordinate transformation | > 0 |
| Smooth xi | `smooth_xi` | Ratio of adjacent step lengths in xi | Close to 1 |
| Smooth eta | `smooth_et` | Ratio of adjacent step lengths in eta | Close to 1 |
| Smooth zeta | `smooth_zt` | Ratio of adjacent step lengths in zeta | Close to 1 |
| Step xi | `step_xi` | Physical step length in xi direction | — |
| Step eta | `step_et` | Physical step length in eta direction | — |
| Step zeta | `step_zt` | Physical step length in zeta direction | — |

## Data Flow

```
prep_grid/ (Python)        C binary              plotting/ (Python)
creat_3d_gauss_     -->    ./main config.json  -->  draw_grid.py
boundary.py                |                       draw_quality.py
(data_file_3d.txt)         v
                      output/ (NetCDF)
                      coord_*.nc
                      orth_*_*.nc, jacobi_*.nc, ...
```

## Input File Formats

### Geometry File (`data_file_3d.txt`)

C-format with integer count lines:

**6-boundary (elliptic):**
```
# nx number
NX
# ny number
NY
# nz number
NZ
# bx1 coords
x y z
x y z
... (ny*nz lines)
# bx2 coords
... (ny*nz lines)
# by1 coords
... (nx*nz lines)
# by2 coords
... (nx*nz lines)
# bz1 coords
... (nx*ny lines)
# bz2 coords
... (nx*ny lines)
```

**2-boundary (parabolic):**
```
# nx ny number
NX NY
# bz1 coords
x y z
... (nx*ny lines)
# bz2 coords
... (nx*ny lines)
```

**1-boundary (hyperbolic, z-direction):**
```
# nx number
NX
# ny number
NY
# bz coords
x y z
... (nx*ny lines)
```

**1-boundary (hyperbolic, x-direction):**
```
# ny number
NY
# nz number
NZ
# bx coords
x y z
... (ny*nz lines)
```

**1-boundary (hyperbolic, y-direction):**
```
# nx number
NX
# nz number
NZ
# by coords
x y z
... (nx*nz lines)
```

### Step File (`step_file_3d.txt`)

```
# number of step
N
# step
value
value
... (nz-1 lines)
```

### Arc-Length File (`arc_len_file*.txt`)

```
# number of points
N
# arc_len
value
value
... (N lines)
```

## Plotting

The `plotting/` directory contains Python scripts for visualization. Edit the `cfs_file` variable at the top of each script to point to the desired output `config.json`:

```python
# In draw_grid.py or draw_quality.py:
cfs_file = '../elliptic/output/config.json'
#cfs_file = '../hyperbolic/output/config.json'
#cfs_file = '../parabolic/output/config.json'
#cfs_file = '../grid-post-process/output/config.json'
```

For 3D visualization, scripts extract 2D slices. Set `subc[n]=1` to select a plane:
- `subc[0]=1` → YOZ plane
- `subc[1]=1` → XOZ plane
- `subc[2]=1` → XOY plane

```python
# Example: plot XOZ plane at y-index 125
subs = [0, 0, 0]
subc = [-1, 1, -1]   # YOZ=-1(all), XOZ=1(single), XOY=-1(all)
subt = [1, 1, 1]
```

For quality plots, also edit `varnm`:
```python
varnm = 'step_zt'  # 'orth_xiet', 'orth_xizt', 'orth_etzt', 'jacobi',
                    # 'smooth_xi', 'smooth_et', 'smooth_zt',
                    # 'step_xi', 'step_et', 'step_zt'
```

## License

See LICENSE file.
