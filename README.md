# PyCellVertex

PyCellVertex is a 2D vertex model simulation package written in Python. It provides a robust, physically grounded framework for simulating multicellular dynamics, tissue morphogenesis, and cellular rearrangements. 

This project is a faithful, modernized Python port of the original C++ codebase [`celldynamics`](https://github.com/yasuhiroinoue/celldynamics). It achieves mathematical and topological parity with the original C++ implementation while offering the flexibility, readability, and rich ecosystem integration of Python.

## 🚀 Dual-Engine Architecture (Educational vs Practical)

PyCellVertex embraces two distinct audiences: researchers who want to **understand the math/algorithms** and those who want to **run large-scale simulations fast**. To support both, the codebase features a Dual-Engine Architecture:

* **`--engine basic` (Source: `src_basic/`) - Default**
  * **Purpose**: Education, algorithmic understanding, and debugging.
  * **Characteristics**: Pure Object-Oriented Python. Mathematical formulas from the vertex model literature are translated directly into readable, clean code. Highly recommended for students and first-time contributors.

* **`--engine numpy` (Source: `src_fast/`)**
  * **Purpose**: High-performance practical simulations.
  * **Characteristics**: Replaces the force calculations and center updates with C-level SIMD operations via NumPy (Flattened Arrays, `np.add.at`, slicing). Vectorization yields a **~6x speedup** on long simulations while retaining the exact same graph topology surgeries.

You can switch engines transparently without altering your configuration:
```bash
# Run educational engine
uv run python main.py --steps 5000 --engine basic

# Run high-performance vectorized engine
uv run python main.py --steps 5000 --engine numpy
```

---

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
  --enable-division --division-time 50.0 --division-stagger-frac 0.5 \
  --dump-vtk --vtk-mode alive --vtk-step 10000
```

### Key Parameters

- `--out-dir`: Directory where simulation logs and VTK artifacts are saved (default: `output/`).
- `--steps`: Total number of simulation steps. The default time step (`DELTA_TIME`) is `1e-4` in simulation time. Thus, `100000` steps corresponds to `10.0` units of simulation time.
- `--num-x`, `--num-y`: Initial grid size of cells.

The simulation behavior is governed by several command-line arguments corresponding to the mathematical and topological properties of the vertex model.

#### Base Mechanics
| Argument | Description | Default |
| :--- | :--- | :--- |
| `--k-area` | Area elasticity constant ($K_a$). Governs resistance to area changes. | `1.0` |
| `--area-eq` | Equilibrium cell area ($A_0$). | `2.60` |
| `--k1-length` | Baseline active line tension ($K_1$). Often negative to promote adhesion. | `-6.0` |
| `--k2-length` | Edge elasticity constant ($K_2$). Prevents edges from becoming infinitely long. | `1.0` |
| `--length-eq` | Equilibrium edge length ($L_0$). | `0.0` |

#### PCP Oscillation
| Argument | Description | Default |
| :--- | :--- | :--- |
| `--k1-pcp` | Amplitude of Planar Cell Polarity (PCP) driven active line tension oscillation. | `0.0` |
| `--power-pcp` | Exponent $n$ for the PCP orientation term, controlling the sharpness of anisotropy. | `2.0` |
| `--pulse-t` | Period of the mechanical oscillation ($T$). | `55.0` |
| `--phase-x`, `--phase-y` | Spatial phase shift components along the X and Y axes. | `0.7854`, `1.5708` |

#### Topology & Division
| Argument | Description | Default |
| :--- | :--- | :--- |
| `--step-reconnect` | Interval (in steps) for checking and executing topological changes (T1/T2). | `100000` |
| `--enable-division` | Flag to enable autonomous cell proliferation. | `True` |
| `--division-time` | Base physical time required for a cell cycle. (Note: $\Delta t = 10^{-4}$). | `1000.0` |
| `--division-stagger-frac` | Fraction `[0, 1]` defining the maximum random phase offset assigned to daughter cells. | `0.25` |

For deep customization (such as modifying hardcoded physical constants like $\Delta t$ or the core force equations), please see [ADVANCED.md](ADVANCED.md).

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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
