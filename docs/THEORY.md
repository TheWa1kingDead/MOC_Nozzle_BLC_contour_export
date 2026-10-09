# Theory and method

This note describes every step of `verify_and_export_coords_full.m`: how the non-dimensional MOC wall is scaled, how the convergent section and the throat are built, and how the boundary layer correction is computed and applied.

Notation: $R^*$ throat radius (half-height for planar nozzles), $D_t = 2R^*$, $D_i$ inlet diameter, $D_e$ exit diameter, $M_e$ design exit Mach number, $\gamma$ ratio of specific heats, $R$ gas constant, $p_0$, $T_0$ stagnation pressure and temperature.

---

## 1. MOC wall contour

MOC_Grid_BDE writes the wall in throat radii, $\bar x = x/R^*$, $\bar y = y/R^*$, with the throat at $(\bar x, \bar y) = (0, 1)$, together with the inviscid wall Mach number $M_w(\bar x)$. These three columns are all the script needs. The block is read from the section headed `Calculated wall contour` in `summary.out`.

## 2. Scaling to physical size

**Isentropic mode** (`scale_mode = 'isentropic'`). The one-dimensional area ratio at $M_e$ is

$$
\frac{A_e}{A^*} = \frac{1}{M_e}\left[\frac{2}{\gamma+1}\left(1+\frac{\gamma-1}{2}M_e^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}
$$

(evaluated with `flowisentropic` from the Aerospace Toolbox). The throat diameter follows from the requested exit diameter:

$$
D_t = \frac{D_e}{\sqrt{A_e/A^*}} \quad \text{(AXI)}, \qquad D_t = \frac{D_e}{A_e/A^*} \quad \text{(TWOD)}
$$

Because the MOC exit radius $\bar y_e$ is not exactly $\sqrt{A_e/A^*}$ (the MOC includes the two-dimensional throat flow), the scaled exit diameter $2\bar y_e R^*$ differs slightly from $D_e$. The script prints both.

**Contour mode** (`scale_mode = 'contour'`). $R^*$ is chosen so the last MOC point is exactly at the requested exit:

$$
R^* = \frac{D_e/2}{\bar y_e}
$$

All coordinates are then $x = \bar x R^*$, $y = \bar y R^*$.

## 3. Throat arc

Upstream of the throat the wall is a circular arc of radius $\bar r_a = 1$ (one throat radius) centred at $(0,\, 2)$:

$$
\bar x_c = \bar r_a \cos\theta, \qquad \bar y_c = \bar r_a \sin\theta + 2, \qquad \theta \in [\theta_{min},\, \theta_{max}] = [225^\circ,\, 270^\circ]
$$

At $\theta = 270^\circ$ the arc reaches the throat $(0,1)$ with zero slope and joins the MOC wall. At $\theta_{min}$ the wall slope is $\tan(-(\theta_{max}-\theta_{min})) = -1$ (a 45 deg convergence half-angle).

The arc radius is the upstream radius of curvature of the throat. It should be consistent with the throat curvature used to set up the transonic solution in the MOC run.

## 4. Convergent section

The convergent wall is a cubic of length $\bar L_s$ (default $\bar L_s = 0.5\,D_i/D_t$ throat radii):

$$
\bar y(s) = a s^3 + b s^2 + \bar y_0, \qquad s \in [0, \bar L_s]
$$

with the four conditions

$$
\bar y(0) = \bar y_0 = \frac{D_i}{D_t},\quad \bar y'(0) = 0,\quad \bar y(\bar L_s) = \bar y_c(\theta_{min}),\quad \bar y'(\bar L_s) = m_L = \tan(-(\theta_{max}-\theta_{min}))
$$

which give, with $\Delta = \bar y(\bar L_s) - \bar y_0$,

$$
a = \frac{m_L}{\bar L_s^2} - \frac{2\Delta}{\bar L_s^3}, \qquad b = -\frac{m_L}{\bar L_s} + \frac{3\Delta}{\bar L_s^2}
$$

The cubic is shifted so it ends at the start of the arc, $\bar x = \bar x_c(\theta_{min})$. The wall is therefore continuous in position and slope from the inlet to the nozzle exit.

## 5. Assembly

Convergent cubic, arc and MOC wall are concatenated, sorted in $x$, and points closer than $10^{-10}$ throat radii are merged (the arc end and the first MOC point are the same throat point).

## 6. Boundary layer displacement thickness

Along the divergent wall the boundary layer edge conditions are the isentropic values at the MOC wall Mach number:

$$
T_e = \frac{T_0}{1+\frac{\gamma-1}{2}M_w^2}, \qquad p_e = p_0\left(\frac{T_e}{T_0}\right)^{\frac{\gamma}{\gamma-1}}, \qquad u_e = M_w\sqrt{\gamma R T_e}
$$

**Wall temperature.** Adiabatic wall, turbulent recovery factor $r = Pr_t^{1/3}$ with $Pr_t = 0.9$:

$$
T_w = T_{aw} = T_e + r\,(T_0 - T_e)
$$

**Reference temperature** (Eckert, 1955):

$$
T^* = \tfrac{1}{2}(T_w + T_e) + 0.22\,(T_{aw} - T_e)
$$

**Viscosity** (Sutherland, $\mu_0 = 1.716\times10^{-5}$ Pa s, $T_{ref} = 273.11$ K, $S = 110.56$ K):

$$
\mu^* = \mu_0\left(\frac{T^*}{T_{ref}}\right)^{3/2}\frac{T_{ref}+S}{T^*+S}
$$

**Reference Reynolds number** based on the running length from a virtual origin $x_0$:

$$
\rho^* = \frac{p_e}{R\,T^*}, \qquad Re^*_x = \frac{\rho^* u_e (x - x_0)}{\mu^*}
$$

**Displacement thickness** (Edenfield, 1968):

$$
\delta^* = 0.42\,(x - x_0)\,\left(Re^*_x\right)^{-0.2775}
$$

**Virtual origin.**

- `x0_mode = 'tangent'`: a straight line is fitted to the first 20% (at least 5 points) of the divergent wall and extrapolated to $y = 0$. The tangent at the throat itself cannot be used because $dy/dx \to 0$ there and the extrapolation goes to $-\infty$.
- `x0_mode = 'inlet'`: $x_0$ is the start of the convergent section.

## 7. Applying the correction

With the local wall angle $\theta_w = \operatorname{atan2}(dy, dx)$ (central differences), each divergent wall point is moved outward along the wall normal:

$$
x_{BLC} = x - \delta^*\sin\theta_w, \qquad y_{BLC} = y + \delta^*\cos\theta_w
$$

In the convergent section the shift is radial only and grows linearly from zero at the inlet to $\delta^*$ at the throat:

$$
y_{BLC} = y + \frac{x - x_{inlet}}{x_{throat} - x_{inlet}}\,\delta^*_{throat}
$$

The corrected wall is checked to remain single-valued in $x$.

The suggested exit diameter for the next design iteration is

$$
D_{e,\,next} = D_e + 2\,\delta^*_{exit}
$$

## 8. CNC re-sampling

Both walls are linearly interpolated onto a uniform axial grid $x_k = x_1 + k\,\Delta x_{cnc}$. For the `prt = 3` file the values are converted to millimetres, rounded to `cnc_decimals` places and repeated rows are removed.

---

## References

- Edenfield, E. E., "Contoured nozzle design and evaluation for hotshot wind tunnels," AIAA Paper 68-369, 1968.
- Eckert, E. R. G., "Engineering relations for friction and heat transfer to surfaces in high velocity flow," *Journal of the Aeronautical Sciences* 22(8), 585-587, 1955.
- Sutherland, W., "The viscosity of gases and molecular force," *Philosophical Magazine* 36(223), 507-531, 1893.
- Rice, T., *2D and 3D Method of Characteristic Tools for Complex Nozzle Development*, JHU/APL Report RTDC-TPS-481, 2003.
