# Example: Mach 3.5 perfect axisymmetric nozzle

`summary.out` is the output of **MOC_Grid_BDE** for its bundled sample case (Mach 3.5, perfect nozzle, axisymmetric, gamma = 1.4). The script reads the `Calculated wall contour` block from it.

With the default parameters in `verify_and_export_coords_full.m` the nozzle is scaled to a 100 mm exit diameter with a 3 inch (76.2 mm) inlet, and the boundary layer is computed for p0 = 10 bar and T0 = 300 K.

## Source and license of this file

This file comes from the sample case `MOC_Grid_BDE/outputs_M3.5Perf/` of the Three-Dimensional Nozzle Design Code, NASA Government Agency Original Software Designation LEW-20180, original author Tharen Rice (JHU/APL):

- Original: <https://github.com/nasa/Three-Dimensional-Nozzle-Design-Code>
- Corrected edition: <https://github.com/TheWa1kingDead/Three-Dimensional-Nozzle-Design-Code>

> Copyright 2020 United States Government as represented by the Administrator of the National
> Aeronautics and Space Administration. No copyright is claimed in the United States under
> Title 17, U.S. Code. All Other Rights Reserved.

It is distributed unmodified under the NASA Open Source Agreement v1.3, included here as [`NOSA_license.txt`](NOSA_license.txt). This does not imply endorsement by NASA or JHU/APL. The MIT License of this repository does not apply to the files in this folder.
