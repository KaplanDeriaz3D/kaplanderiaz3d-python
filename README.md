# KaplanDeriaz3D — Python App

<img width="1920" height="994" alt="image" src="https://github.com/user-attachments/assets/dede5028-0404-4fb4-bb42-24aaf390ea2d" />

<img width="1918" height="948" alt="image" src="https://github.com/user-attachments/assets/46a2c685-c767-4ac7-8a26-f4d6e44b1afe" />

A Python/PySide6 + PyVista application for designing runner blades of **Kaplan** (axial-flow) and **Deriaz** (diagonal-flow) hydraulic turbines — no MATLAB license required.

Built as part of a Master's Thesis (TFM) on hydraulic machinery design. See [kaplanderiaz3d-theory](https://github.com/KaplanDeriaz3D/kaplanderiaz3d-theory) for the underlying theory and thesis documentation, and [kaplanderiaz3d-matlab](https://github.com/KaplanDeriaz3D/kaplanderiaz3d-matlab) for the original MATLAB version.

## Features

**Turbine hydraulic design**
- Supports both **Kaplan** (axial-flow) and **Deriaz** (diagonal-flow) runner types, switchable from a single interface.
- Computes the full blade backbone from operating parameters: rotational speed (RPM), rotation direction, flow rate, net head, gravity, fluid density, hydraulic/volumetric/mechanical efficiencies, and solidity coefficient (σ).
- Configurable interpolation law for the blade angle β along the chord (cubic, linear, cosine, inlet-loaded or outlet-loaded).
- Geometry inputs specific to each turbine type (Kaplan hub/tip radius and axial length; Deriaz spherical radii and cone angles).

**Hydrofoil & 3D mesh generation**
- Generates NACA-based (among others) hydrofoil profiles per radial station, with configurable thickness/chord ratio, hub-to-tip thickness taper and mesh resolution.
- Builds a full **solid 3D blade** (not just a mean surface). Thickness is applied along each streamline's own surface of revolution (cylinder for Kaplan, sphere for Deriaz), so even very thick blades stay flush with the hub and shroud.

**Interactive 3D visualization**
- Dual PyVista viewers: mean-surface (multi-blade) view and solid-runner (airfoil profile) view, matching the MATLAB app's layout.
- Light and Dark visualization themes, togglable from the menu.
- Reset 3D view and camera controls.

**Export**
- **STEP (default)**: parametric export with exact B-spline (NURBS) surfaces that pass through every computed point — no faceting or loss of quality. With OpenCASCADE installed (`cadquery-ocp`), the blade is written as a **closed B-Rep solid** ready for any CAD package.
- **IGES**: exact B-spline surfaces (IGES 5.3, entity 128).
- **STL**: watertight triangulated mesh with consistent outward normals, for 3D printing or meshing tools.
- Choose between the full solid blade or the mean surface only. STEP/IGES are written in millimetres; STL in metres.
- **Excel export**: inlet/outlet blade angles (β1, β2) per radial station, exported to a `.xlsx` workbook.
- **Copy Results Table**: copy the computed results table directly to the clipboard.

**Configuration management**
- **Save/Load Configuration**: store all current design parameters to a `.json` file and reload them later, to reproduce or share a specific design.
- **Reset to Defaults**: restore all parameters to their default values in one click.

## Requirements

- Python 3 (tested with Python 3.12)
- Dependencies listed in [`requirements.txt`](requirements.txt): PySide6, PyVista, PyVistaQt, NumPy, openpyxl and cadquery-ocp
- `cadquery-ocp` (OpenCASCADE) is only needed to export STEP as a closed solid. Without it, the app still runs and STEP/IGES are exported as exact B-spline surfaces.

> **Note:** the code is written against PySide6, but is also compatible with PyQt6 with minimal changes (see the note at the top of the source file).

## How to run

1. Clone or download this repository.
2. (Recommended) Create a virtual environment:
```bash
   python -m venv venv
   source venv/bin/activate      # on Windows: venv\Scripts\activate
```
3. Install the dependencies:
```bash
   pip install -r requirements.txt
```
4. Run the app:
```bash
   python kaplan_deriaz3d.py
```

## License

Licensed under [GPLv3](LICENSE).

## Part of the KaplanDeriaz3D project

[Organization overview](https://github.com/KaplanDeriaz3D) · [MATLAB app](https://github.com/KaplanDeriaz3D/kaplanderiaz3d-matlab) · [Python app](https://github.com/KaplanDeriaz3D/kaplanderiaz3d-python) · [Theory & TFM](https://github.com/KaplanDeriaz3D/kaplanderiaz3d-theory)
