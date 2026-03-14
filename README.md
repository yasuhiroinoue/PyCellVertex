# PyCellVertex

PyCellVertex is a 2D vertex model simulation package written in Python. It provides a robust, physically grounded framework for simulating multicellular dynamics, tissue morphogenesis, and cellular rearrangements. 

This project is a faithful, modernized Python port of the original C++ codebase [`celldynamics`](https://github.com/yasuhiroinoue/celldynamics). It achieves mathematical and topological parity with the original C++ implementation while offering the flexibility, readability, and rich ecosystem integration of Python.

## Features

- **2D Vertex Model**: Simulates epithelial tissue dynamics using a mechanically consistent vertex model.
- **Cell Division**: Simulates autonomous cell proliferation with staggered timing and proper topological inheritance.
- **Topological Surgery**: Accurately handles cell intercalation (T1 transitions) and complex boundary intersections (T2/T3-like transitions) while strictly preserving graph topology.
- **PCP (Planar Cell Polarity) Oscillation**: Supports anisotropic line tensions and active oscillatory mechanics.
- **VTK Visualization**: Outputs simulation states in VTK format for direct visualization and analysis in tools like ParaView.

## Requirements

- Python 3.11+
- [uv](https://github.com/astral-sh/uv) (recommended for dependency management)

## Installation & Usage

Clone the repository and run the simulation using `uv`:

```bash
git clone https://github.com/yasuhiroinoue/PyCellVertex.git
cd PyCellVertex
```

### Basic Run (All mechanics enabled)
```bash
uv run python main.py --out-dir output/basic_run --steps 100000 --dump-vtk --vtk-mode alive
```

### Example: PCP Oscillation Only (No Division)
To observe pure tissue deformation and cell intercalations driven by Planar Cell Polarity (PCP) active forces, disable cell division:
```bash
uv run python main.py \
  --out-dir output/pcp_only \
  --steps 5000000 \
  --k1-pcp 0.2 --power-pcp 2.0 --pulse-t 55.0 \
  --no-enable-division \
  --dump-vtk --vtk-mode alive --vtk-step 10000
```

### Example: Cell Proliferation Only (No PCP)
To observe tissue growth through cell division without anisotropic active forces:
```bash
uv run python main.py \
  --out-dir output/division_only \
  --steps 5000000 \
  --k1-pcp 0.0 \
  --enable-division --division-time 1000.0 --division-stagger-frac 0.25 \
  --dump-vtk --vtk-mode alive --vtk-step 10000
```

### Key Parameters
- `--out-dir`: Directory where simulation logs and VTK artifacts are saved (default: `output/`).
- `--steps`: Total number of simulation steps. The default time step (`DELTA_TIME`) is `1e-4` in simulation time. Thus, `100000` steps corresponds to `10.0` units of simulation time.
- `--num-x`, `--num-y`: Initial grid size of cells.
- `--k-area`, `--area-eq`: Area elasticity constant and equilibrium area.
- `--k1-pcp`, `--pulse-t`: Parameters governing PCP-driven mechanical oscillations.
- `--enable-rearrange`, `--enable-intersection`: Toggle topological modifications.
- `--enable-division`: Enable cell proliferation.
- `--division-time`: The baseline physical time required for a cell to divide.
- `--division-stagger-frac`: A fraction `[0, 1]` that defines the maximum random offset for the cell cycle clock. The initial cell time is randomized uniformly in the range `[0, division-time * division-stagger-frac]`. This offset is also freshly re-assigned to both daughter cells every time a cell divides, ensuring division timings remain staggered (asynchronous) throughout the tissue.

## Visualization

Simulation outputs are generated in the `artifacts/` directory by default (or the path specified by `--out-dir`).

### VTK Output Modes (`--vtk-mode`)
The simulation provides dual-mode VTK outputs to separate physical visualization from internal topological debugging:
- `--vtk-mode alive` (Recommended): Outputs only active, topologically valid structures. Use this mode for all standard visualizations and analysis.
- `--vtk-mode both`: Outputs both `alive` and `raw` data into separate subdirectories (`vtk_alive/` and `vtk_raw/`).
- `--vtk-mode raw`: Includes logically deleted "dead" vertices and edges. Strictly for internal topology debugging.

### Loading in ParaView
Load the resulting `.vtk` files from the **`vtk_alive/`** directory into [ParaView](https://www.paraview.org/):
- **`2dv_line*.vtk`**: Contains cell edge (boundary) information.
- **`2dv_face*.vtk`**: Contains cell polygon (face) information and cell phase dynamics.
