# Changelog

All notable changes to this project are listed here. Versions follow [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-10-10

First public release.

- Reads the wall contour directly from the MOC_Grid_BDE `summary.out` (`Calculated wall contour` block), or from a table holding only that block.
- Scales the non-dimensional contour to physical size from the exit diameter, by the isentropic area ratio (`scale_mode = 'isentropic'`) or by the MOC exit radius (`scale_mode = 'contour'`).
- Builds the full wall: cubic convergent section, circular throat arc and MOC divergent section, continuous in slope.
- Edenfield (1968) turbulent boundary layer displacement correction with Eckert reference temperature and Sutherland viscosity; choice of virtual origin (`x0_mode`).
- Exports inviscid and BL-corrected coordinates for meshing (m) and CNC machining (m, or mm with fixed precision). Output files are named after the input file.
- Checks: exit Mach against `Me`, throat larger than inlet, convergent overshoot, monotonic BL-corrected wall.
- Bundled example: MOC_Grid_BDE Mach 3.5 perfect axisymmetric nozzle.
