# PyCellVertex

PyCellVertex is a 2D vertex model simulation package written in Python. It provides a robust, physically grounded framework for simulating multicellular dynamics, tissue morphogenesis, and cellular rearrangements.

## Features

- **2D Vertex Model**: Simulates epithelial tissue dynamics using a mechanically consistent vertex model.
- **Topological Surgery**: Accurately handles cell division, cell intercalation (T1 transitions), and complex boundary intersections (T2/T3-like transitions) while strictly preserving graph topology.
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
uv run python main.py --steps 100000 --dump-vtk
```

### Key Parameters
- `--steps`: Total number of simulation steps.
- `--num-x`, `--num-y`: Initial grid size of cells.
- `--k-area`, `--area-eq`: Area elasticity constant and equilibrium area.
- `--k1-pcp`, `--pulse-t`: Parameters governing PCP-driven mechanical oscillations.
- `--enable-rearrange`, `--enable-intersection`: Toggle topological modifications.

## Visualization

Simulation outputs are generated in the `artifacts/` directory by default. 
Load the resulting `.vtk` files into [ParaView](https://www.paraview.org/) to visualize cell boundaries and phase dynamics.
