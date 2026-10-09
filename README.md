# MOC Nozzle Contour Export

MATLAB tool that turns the wall contour computed by the method of characteristics program **MOC_Grid_BDE** (`summary.out`) into a complete, physically scaled supersonic nozzle wall. It adds a convergent section and a circular throat, applies a turbulent boundary layer displacement correction (Edenfield, 1968) to the divergent wall, and writes coordinate files ready for:

- **meshing / CFD** (the assembled wall points, in metres), and
- **CNC machining** (uniformly re-sampled wall, in metres or millimetres),

each in an **inviscid** and a **boundary-layer-corrected (BLC)** version. Axisymmetric (`AXI`) and planar (`TWOD`) nozzles are supported.

The MOC program itself is in the companion repository
[Three-Dimensional-Nozzle-Design-Code](https://github.com/TheWa1kingDead/Three-Dimensional-Nozzle-Design-Code) (corrected, standalone-build edition of NASA LEW-20180).

---

## Contents

- [Workflow](#workflow)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Input file](#input-file)
- [User parameters](#user-parameters)
- [What the script does](#what-the-script-does)
- [Outputs](#outputs)
- [Console summary](#console-summary)
- [Figures](#figures)
- [Checks and warnings](#checks-and-warnings)
- [Assumptions and limitations](#assumptions-and-limitations)
- [Repository layout](#repository-layout)
- [How to cite](#how-to-cite)
- [References](#references)
- [License](#license)

Equations and the full derivation of every step are in [`docs/THEORY.md`](docs/THEORY.md).

---

## Workflow

```
MOC_Grid_BDE  (design Mach, gamma, nozzle type, throat curvature)
     |
     v
summary.out   "Calculated wall contour" block: J, X/R*, R/R*, Mach, ...
     |
     v
verify_and_export_coords_full.m
     |-- scale to physical size (throat radius R* from De and Me)
     |-- convergent cubic + circular throat arc + MOC divergent wall
     |-- Edenfield turbulent BL displacement thickness, wall-normal shift
     |-- re-sample at a fixed axial step for CNC
     v
<name>_mesh_inviscid.txt   <name>_mesh_blc.txt     (m)
<name>_cnc_inviscid.txt    <name>_cnc_blc.txt      (m)
<name>_cnc_blc_mm.txt                              (mm, rounded)
```

A typical design loop is:

1. Run MOC_Grid_BDE for the design exit Mach number.
2. Run this script with the required exit diameter `De`.
3. Read the *Suggested De for next MOC run* in the console. If the BL-corrected exit is too large for the test section or the hardware, change `De` and repeat step 2.
4. Export (`prt = 4`) and send the `cnc` file to the workshop and the `mesh` file to the mesh generator.

---

## Requirements

- MATLAB **R2019b** or newer (`tiledlayout`, `writematrix`, `isfile`).
- **Aerospace Toolbox**, for `flowisentropic` (isentropic area ratio).
- Optimization of the cubic connector uses `fzero`, `fminbnd` and `integral` from base MATLAB.

---

## Quick start

1. Clone or download this repository.
2. Open `verify_and_export_coords_full.m` in MATLAB and press **Run**.

With no changes it reads the bundled example
[`examples/M3.5_perfect_axi/summary.out`](examples/M3.5_perfect_axi/) (Mach 3.5 perfect axisymmetric nozzle), builds a nozzle with a 100 mm exit diameter and a 3 inch inlet, prints the BL summary and draws the four plots. Nothing is written to disk until `prt` is set.

To use your own MOC result, either set `pname` and `fname` at the top of the script, or set `fname = ''` and pick the file in the dialog that opens.

---

## Input file

The script accepts either of these:

| File | What is read |
|---|---|
| `summary.out` written by MOC_Grid_BDE | The block that starts at the line `Calculated wall contour`. The column header is skipped and rows are read until the first blank line. |
| A plain table (any extension) | Only the wall block, copied out of `summary.out`. Header lines are allowed. |

Columns used (all from the MOC wall block):

| Column | Quantity | Units |
|---|---|---|
| 1 | J (point index) | - |
| 2 | X/R* | throat radii |
| 3 | R/R* (half-height h/h* for `TWOD`) | throat radii |
| 4 | Mach number on the wall | - |

The remaining columns (angle, P, T, density, mass flow) are not needed. Coordinates are non-dimensional, so the same MOC run can be scaled to any size.

---

## User parameters

All inputs are in the *Control parameters* section at the top of the script.

| Variable | Default | Units | Meaning |
|---|---|---|---|
| `prt` | `0` | - | Export selector: `0` none, `1` BLC (CNC + mesh), `2` inviscid (CNC + mesh), `3` BLC CNC in mm, rounded, `4` everything |
| `geom` | `'AXI'` | - | `'AXI'` axisymmetric or `'TWOD'` planar |
| `Me` | `3.5` | - | Design exit Mach number, must match the MOC run |
| `g` | `1.4` | - | Ratio of specific heats |
| `Rgas` | `287` | J/(kg K) | Gas constant |
| `Di` | `3*25.4e-3` | m | Inlet diameter (inlet height for `TWOD`) |
| `De` | `100e-3` | m | Exit diameter (exit height for `TWOD`) |
| `T0` | `300` | K | Stagnation temperature (for the BL calculation) |
| `P0` | `10e5` | Pa | Stagnation pressure (for the BL calculation) |
| `dx_cnc` | `10e-5` | m | Axial step of the CNC files (0.1 mm) |
| `cnc_decimals` | `2` | - | Decimal places (mm) kept in the `prt = 3` file |
| `scale_mode` | `'isentropic'` | - | How R* is chosen, see below |
| `x0_mode` | `'tangent'` | - | Where the boundary layer starts, see below |

Geometry parameters further down the script:

| Variable | Default | Meaning |
|---|---|---|
| `th_min`, `th_max`, `d_th` | `225`, `270`, `1` deg | Start, end and step of the circular throat arc |
| `R_arc` | `1` | Arc radius in throat radii (upstream radius of curvature of the throat) |
| `L_spline` | `0.5*Di/Dt_r` | Length of the convergent cubic, in throat radii |

**`scale_mode`**

- `'isentropic'`: the throat diameter is computed from the isentropic area ratio at `Me`, `Dt = De/sqrt(Ae/A*)` (AXI) or `Dt = De/(Ae/A*)` (TWOD). The exit of the scaled contour is then close to, but not exactly, `De`, because the MOC exit radius differs slightly from the one-dimensional value. Both values are printed.
- `'contour'`: the throat diameter is chosen so that the last MOC wall point lands exactly on `De/2`.

**`x0_mode`**

- `'tangent'`: virtual origin from a straight-line fit to the first 20% of the divergent wall, extrapolated to the axis.
- `'inlet'`: the boundary layer starts at the beginning of the convergent section.

---

## What the script does

1. **Read** the MOC wall contour (`X/R*`, `R/R*`, `M`) and check that the exit Mach number matches `Me`.
2. **Throat size.** Compute the physical throat diameter `Dt` from `De` (see `scale_mode`).
3. **Circular throat arc.** Arc of radius `R_arc*R*` centred at `(0, 2R*)`, from `th_min` to `th_max`. At 270 deg it meets the MOC wall at the throat `(0, R*)` with zero slope.
4. **Convergent section.** Cubic `y(x)` from the inlet radius `Di/2` (zero slope) to the start of the arc, matching the arc slope `tan(-(th_max - th_min))`, so the wall is smooth (C1) at both ends.
5. **Assemble** convergent cubic, arc and MOC divergent wall; sort, remove coincident points and scale by `R* = Dt/2`.
6. **Boundary layer.** Along the divergent wall, edge conditions come from the isentropic relations with the MOC wall Mach number. The Edenfield correlation with Eckert reference temperature and Sutherland viscosity gives the displacement thickness `delta*`. The wall is moved outward along its local normal by `delta*`. In the convergent section the radial shift is ramped linearly from zero at the inlet to `delta*` at the throat.
7. **CNC re-sampling.** Both walls are interpolated onto a uniform axial grid with step `dx_cnc`.
8. **Report, plot and export.**

Full equations: [`docs/THEORY.md`](docs/THEORY.md).

---

## Outputs

Files are written next to the input file and are named after it (`summary.out` gives `summary_*.txt`). All are tab-separated, three columns `x  y  z` with `z = 0`, so they import directly as a 3D curve in CAD, CAM and meshing tools.

| `prt` | File | Content | Units | Typical use |
|---|---|---|---|---|
| 1, 4 | `<name>_cnc_blc.txt` | BLC wall, uniform `dx_cnc` step | m | CAM / CNC |
| 1, 4 | `<name>_mesh_blc.txt` | BLC wall, assembled points | m | Mesh generator (Gambit, ICEM, Pointwise, Gmsh) |
| 2, 4 | `<name>_cnc_inviscid.txt` | Inviscid wall, uniform `dx_cnc` step | m | CAM / CNC, inviscid CFD |
| 2, 4 | `<name>_mesh_inviscid.txt` | Inviscid wall, assembled points | m | Mesh generator, inviscid (Euler) CFD |
| 3, 4 | `<name>_cnc_blc_mm.txt` | BLC wall, uniform step, rounded to `cnc_decimals`, repeated rows removed | mm | Direct import into CAM software |

Which wall to use:

- **Viscous CFD or the real hardware**: use the **BLC** wall. Its inner surface is displaced outward so that the inviscid core sees the MOC contour.
- **Inviscid (Euler) CFD**: use the **inviscid** wall.

`x = 0` is at the throat; the convergent section has negative `x`. `y` is the radius (AXI) or the half-height (TWOD).

---

## Console summary

```
 Boundary Layer Correction Summary (AXI, Me = 3.50)
  Input file                       :  ...\examples\M3.5_perfect_axi\summary.out
  Scaling mode                     :  isentropic
  Throat diameter (inviscid)       :  ... mm
  BL virtual origin  x0 (tangent)  : ... mm
  delta* at throat                 :  ... mm
  delta* at exit                   :  ... mm
  Requested exit diameter (De)     :  ... mm
  Inviscid exit diameter (contour) :  ... mm
  BLC-corrected exit diameter      :  ... mm
  Suggested De for next MOC run    :  ... mm

  BLC throat diameter    :  ... mm
  Total nozzle length    :  ... mm
```

*Suggested De* is `De + 2*delta*_exit`: the exit size the hardware will have once the boundary layer correction is added.

---

## Figures

One window with four tiles:

1. Assembled inviscid wall (convergent, arc and MOC divergent sections), equal axes.
2. Inviscid wall and BLC wall together.
3. `delta*` (left axis) and wall Mach number (right axis) along the divergent section.
4. The re-sampled CNC walls, inviscid and BLC.

---

## Checks and warnings

| Message | Cause | What to do |
|---|---|---|
| `Exit Mach in the MOC file ... differs from Me` | `Me` does not match the MOC run | Set `Me` to the MOC design Mach number |
| `Diameter of throat > Diameter of inlet` | `Di` smaller than the throat | Increase `Di` or reduce `De` |
| `Sensible tangential connection to circular throat failed` | The convergent cubic overshoots above the inlet radius | Change `Di`, `L_spline` or `th_min` |
| `BL-corrected wall x is not monotonic` | `delta*` shift larger than the wall point spacing | Check `P0`, `T0` and the nozzle size |
| `Could not hit length exactly` | Cubic connector did not converge to the requested arc length | Change `L_circ` |

---

## Assumptions and limitations

- Calorically perfect gas; edge conditions are isentropic from `P0`, `T0`.
- The BL correlation is for a **turbulent**, **adiabatic** wall boundary layer on a flat plate (Edenfield, 1968). It does not include the effect of pressure gradient, wall cooling, transition or the transverse curvature of small nozzles. For high-Mach, cold-wall or low-Reynolds-number nozzles, check the result against a boundary layer code or viscous CFD.
- `delta*` is applied as a displacement of the wall along its normal. The correction is not iterated with the core flow.
- The convergent section shape (cubic) and the throat arc are geometric choices, not the result of a flow calculation. The upstream radius of curvature of the throat (`R_arc*R*`) should be consistent with the throat curvature used in the MOC run.
- `connect_line_circle_cubic` (a cubic connector with a prescribed arc length) is evaluated but its result is not used in the assembled wall; the convergent wall comes from `conv_cont2`.
- The CNC files are a polyline. The tool path tolerance is set by `dx_cnc` and `cnc_decimals`.

---

## Repository layout

```
|-- verify_and_export_coords_full.m   main script (with local functions)
|-- docs/
|   `-- THEORY.md                     geometry, scaling and BL equations
|-- examples/
|   `-- M3.5_perfect_axi/             MOC_Grid_BDE sample output, Mach 3.5 perfect nozzle
|-- CHANGELOG.md
|-- CITATION.cff                      citation metadata (GitHub "Cite this repository")
`-- LICENSE                           MIT
```

---

## How to cite

If you use this tool, please cite it (GitHub **"Cite this repository"** uses [`CITATION.cff`](CITATION.cff)), together with the MOC program used to compute the contour and the Edenfield boundary layer method.

```bibtex
@software{reddy_moc_nozzle_contour_export_2026,
  author  = {Reddy, MRK},
  title   = {{MOC Nozzle Contour Export}: physically scaled inviscid and
             boundary-layer-corrected nozzle coordinates for meshing and
             CNC machining},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/TheWa1kingDead/MOC-Nozzle-Contour-Export}
}
```

---

## References

- Edenfield, E. E., "Contoured nozzle design and evaluation for hotshot wind tunnels," AIAA Paper 68-369, 1968.
- Eckert, E. R. G., "Engineering relations for friction and heat transfer to surfaces in high velocity flow," *Journal of the Aeronautical Sciences* 22(8), 585-587, 1955.
- Rice, T., *2D and 3D Method of Characteristic Tools for Complex Nozzle Development*, JHU/APL Report RTDC-TPS-481, 2003. [NTRS 20030067852](https://ntrs.nasa.gov/citations/20030067852)
- Reddy, MRK, *Three-Dimensional Nozzle Design Code* (corrected edition of NASA LEW-20180), 2026. [doi:10.5281/zenodo.23269365](https://doi.org/10.5281/zenodo.23269365)

---

## License

The MATLAB code is released under the [MIT License](LICENSE).

The example file in [`examples/M3.5_perfect_axi/`](examples/M3.5_perfect_axi/) is output of the NASA Three-Dimensional Nozzle Design Code sample case and is distributed under the NASA Open Source Agreement v1.3; see the README in that folder.

## Acknowledgements

This tool was developed by MRK Reddy under the guidance of **Prof. S.K. Karthick**, Department of Mechanical and Aerospace Engineering, Indian Institute of Technology Hyderabad, India.

The example nozzle contour was computed with the Three-Dimensional Nozzle Design Code, originally developed by Tharen Rice at the Johns Hopkins University Applied Physics Laboratory with funding from NASA Glenn Research Center.

## Contact

**MRK Reddy** · ORCID [0009-0006-8420-6672](https://orcid.org/0009-0006-8420-6672) · GitHub [@TheWa1kingDead](https://github.com/TheWa1kingDead)

Bug reports and suggestions: [Issues](../../issues).
