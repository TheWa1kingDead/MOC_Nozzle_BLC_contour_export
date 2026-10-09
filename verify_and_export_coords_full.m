%% Full nozzle coordinate generator with BL correction
% Scales the MOC_Grid_BDE wall contour (summary.out) to physical
% dimensions, adds a convergent section and a circular throat arc, then
% applies the Edenfield (1968) turbulent boundary layer displacement
% correction to the divergent wall. Writes inviscid and BL-corrected
% coordinates for meshing and for CNC machining.
% Supports axisymmetric (AXI) and planar (TWOD) nozzles.
%
% Requires: MATLAB R2019b or newer, Aerospace Toolbox (flowisentropic).
% See README.md for a full description of inputs, method and outputs.
%
% MRK Reddy, 2026. MIT License, see LICENSE.

clc; close all; clearvars;

set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');


%% File location and load data

% Folder and name of the MOC_Grid_BDE output. Either the full summary.out
% written by MOC_Grid_BDE, or a plain table holding only the "Calculated
% wall contour" block (columns J, X/R*, R/R*, Mach, ...).
% Leave fname empty ('') to pick the file with a dialog.
pname = fullfile(fileparts(mfilename('fullpath')), 'examples', 'M3.5_perfect_axi');
fname = 'summary.out';

if isempty(fname) || ~isfile(fullfile(pname, fname))
    [fname, pname] = uigetfile({'*.out;*.OUT;*.txt;*.dat', 'MOC wall files'; '*.*', 'All files'}, ...
                               'Select MOC_Grid_BDE summary.out or wall contour file');
    if isequal(fname, 0)
        error('No input file selected.');
    end
end

[x, y, M_wall] = read_moc_wall(fullfile(pname, fname));
%   x      : nondimensional axial  coordinate  x/R*
%   y      : nondimensional radial coordinate  y/R*  (half-height for TWOD)
%   M_wall : Mach number along wall

[~, base_name] = fileparts(fname);       % prefix for all exported files
out_dir        = pname;                  % exported files go next to the input


%% Control parameters

prt   = 0;                               % 0 = no export
                                         % 1 = BLC export (CNC + mesh)
                                         % 2 = inviscid export (CNC + mesh)
                                         % 3 = BLC CNC file in mm, rounded, duplicates removed
                                         % 4 = all of the above
geom  = 'AXI';                           % 'AXI' (axisymmetric) or 'TWOD' (planar)

Me    = 3.5;                             % Design exit Mach number (same as the MOC run)
g     = 1.4;                             % Gamma
Rgas  = 287;                             % Gas constant (J/kg/K)

Di    = 3*25.4e-3;                       % Inlet diameter (m)   (inlet height for TWOD)
De    = 100e-3;                          % Exit diameter (m)    (exit height for TWOD)

T0    = 300;                             % Stagnation temperature (K)
P0    = 10e5;                            % Stagnation pressure (Pa)

dx_cnc       = 10e-5;                    % CNC output resolution (m)
cnc_decimals = 2;                        % Decimals (mm) kept in the prt = 3 CNC file

scale_mode = 'isentropic';               % 'isentropic' : R* from the isentropic Ae/A* at Me (original behaviour)
                                         % 'contour'    : R* so that the MOC exit radius equals De/2 exactly
x0_mode    = 'tangent';                  % BL virtual origin:
                                         % 'tangent' : divergent-wall tangent extrapolated to the axis (original behaviour)
                                         % 'inlet'   : start of the convergent section


%% Check the MOC data against the design Mach number

if abs(M_wall(end) - Me) > 0.01
    warning('Exit Mach in the MOC file (%.4f) differs from Me (%.4f). Check Me.', M_wall(end), Me);
end


%% Throat diameter

[~, ~, ~, ~, a_rat] = flowisentropic(g, Me);

if strcmp(scale_mode, 'contour')
    Dt_r = De / y(end);                  % MOC exit radius y(end)*R* = De/2
elseif strcmp(geom, 'AXI')
    Dt_r = De / sqrt(a_rat);             % AXI  : Ae/A* = (De/Dt)^2
elseif strcmp(geom, 'TWOD')
    Dt_r = De / a_rat;                   % TWOD : Ae/A* = De/Dt  (h_e/h_t)
else
    error('geom must be ''AXI'' or ''TWOD''.');
end

%% Circular throat

th_min = 225;                            % Throat start angle (deg) [Default is 225]
th_max = 270;                            % Throat end angle (deg)
d_th   = 1;                              % Angular step

th_rng = th_min : d_th : th_max;
R_arc  = 1;                              % nondimensional arc radius
xc     = R_arc * cosd(th_rng);           % arc x-coords
yc     = R_arc * sind(th_rng) + R_arc;   % arc y-coords

dth = th_max - th_min;                   % total arc angle


y_line  = 2*(Di/Dt_r);                   % inlet wall height
cx = 0; cy = 2;
R_circ = (Dt_r/Dt_r)/2;                  % circle parameters for cubic connector
xb_c = xc(1); yb_c = yc(1);              % tangent-start point on throat arc
L_circ = 2.5*Di/Dt_r;                    % connector arc length

out = connect_line_circle_cubic(y_line, cx, cy, R_circ, xb_c, yb_c, L_circ);

%% Convergent spline contour

D1_conv  = 2*(Di/Dt_r);                  % inlet diameter
Dt_conv  = 2*(yc(1) + 1);                % throat-exit diameter
L_spline = 0.5*Di/Dt_r;                  % convergent spline length (default = 0.5)
mL       = tand(-dth);                   % slope at convergent exit

wall_points = conv_cont2(D1_conv, Dt_conv, L_spline, mL);

%% Assemble and scale to physical dimensions

xf = [wall_points(:,1) - (L_spline + abs(xc(1))); xc'; x];
yf = [wall_points(:,2); yc'+1; y];

% Sort and drop coincident points. The arc end (270 deg) and the first
% MOC point are both the throat (0, 1); without a tolerance they survive
% as two points 1e-16 apart.
[xf, isrt] = sort(xf);
yf         = yf(isrt);
keep       = [true; diff(xf) > 1e-10];
xb_nd      = xf(keep);
yb_nd      = yf(keep);

xb = xb_nd * (Dt_r/2);                   % scale to physical (m)
yb = yb_nd * (Dt_r/2);

%% Edenfield 1968 Boundary layer correction

% Locate divergent section in assembled coordinates.
% Divergent section = MOC data, starts at the first x(1) in non-dim coords.
x_div_start_nd = x(1);
idx_div        = find(xb_nd >= x_div_start_nd - 1e-10, 1, 'first');

xb_div = xb(idx_div:end);                % divergent wall x [m]
yb_div = yb(idx_div:end);                % divergent wall y [m]

% Mach distribution along divergent wall.
% Interpolate from the MOC (x, M_wall) to the assembled unique points.
x_div_nd = xb_nd(idx_div:end);
M_div    = interp1(x, M_wall, x_div_nd, 'linear', 'extrap');

% Isentropic conditions at BL edge (= wall value of the inviscid solution)
T_e = T0 ./ (1 + (g-1)/2 .* M_div.^2);            % static temperature  [K]
p_e = P0 .* (T_e ./ T0).^(g/(g-1));               % static pressure     [Pa]

% BL virtual origin x0
if strcmp(x0_mode, 'inlet')
    x0_bl = xb(1);                                % BL starts at the convergent inlet [m]
else
    % Back-extrapolate the divergent wall tangent to the axis (y=0).
    % The first 1-2 points cannot be used: dy/dx is about 0 at the throat,
    % which sends x0 to -inf. Instead fit a line over the first 20% of the
    % divergent section and extrapolate that tangent to y=0.
    n_fit  = max(5, round(length(xb_div)*0.20));            % first 20% of divergent section
    p_fit  = polyfit(xb_div(1:n_fit), yb_div(1:n_fit), 1);  % p(1)=slope, p(2)=intercept
    x0_bl  = -p_fit(2) / p_fit(1);                          % x where fitted line hits y=0  [m]
end

% Displacement thickness along divergent wall
delta_star = blc_edenfield_wall(xb_div, x0_bl, p_e, M_div, T_e, T0, g, Rgas);

% Wall-normal outward shift.
% Wall angle at each point via central differences (gradient handles endpoints)
theta_w   = atan2(gradient(yb_div), gradient(xb_div));    % [rad]

x_blc_div = xb_div - delta_star .* sin(theta_w);          % shift upstream  [m]
y_blc_div = yb_div + delta_star .* cos(theta_w);          % shift outward   [m]

% Convergent section: linear y-ramp to match the throat junction.
% Ramps from zero shift at inlet to delta_star(1) at throat, radially only.
x_conv = xb(1:idx_div-1);
y_conv = yb(1:idx_div-1);
frac         = (x_conv - x_conv(1)) / (x_conv(end) - x_conv(1));
y_blc_conv   = y_conv + frac * delta_star(1);
x_blc_conv   = x_conv;

% Assemble full BLC wall and interpolate for CNC
xb_blc = [x_blc_conv; x_blc_div];
yb_blc = [y_blc_conv; y_blc_div];

if any(diff(xb_blc) <= 0)
    error(['BL-corrected wall x is not monotonic (delta* too large for the wall spacing). ', ...
           'Check P0, T0 and the size of the nozzle.']);
end

xi_blc = xb_blc(1) : dx_cnc : xb_blc(end);
yi_blc = interp1(xb_blc, yb_blc, xi_blc);


%% CNC interpolation

% Inviscid CNC interpolation
xi = xb(1) : dx_cnc : xb(end);
yi = interp1(xb, yb, xi);

% BLC CNC interpolation
xv = xi_blc(1) : dx_cnc : xi_blc(end);
yv = interp1(xi_blc, yi_blc, xv);

%% Sanity checks for geometry

err_throat = Dt_r / Di;
if err_throat > 1
    warning('Diameter of throat > Diameter of inlet. Proceeding anyway.');
end

err_spline = yb(1) / (max(wall_points(:,2) * Dt_r * 0.5));
if err_spline < 1
    warning('Sensible tangential connection to circular throat failed, change diameters or circular throat angle.');
end

%% BL correction outputs

delta_exit   = delta_star(end);             % delta* at nozzle exit (m)
De_inv       = 2 * yb(end);                 % exit diameter of the scaled inviscid contour (m)
De_eff       = 2 * y_blc_div(end);          % effective exit diameter with BLC (m)
De_suggested = De + 2*delta_exit;           % suggested De for next iteration (m)

fprintf('\n Boundary Layer Correction Summary (%s, Me = %.2f)\n', geom, Me);
fprintf('  Input file                       :  %s\n', fullfile(pname, fname));
fprintf('  Scaling mode                     :  %s\n', scale_mode);
fprintf('  Throat diameter (inviscid)       :  %.4f mm\n', Dt_r*1e3);
fprintf('  BL virtual origin  x0 (%s)  : %+.4f mm\n', x0_mode, x0_bl*1e3);
fprintf('  delta* at throat                 :  %.4f mm\n', delta_star(1)*1e3);
fprintf('  delta* at exit                   :  %.4f mm\n', delta_exit*1e3);
fprintf('  Requested exit diameter (De)     :  %.4f mm\n', De*1e3);
fprintf('  Inviscid exit diameter (contour) :  %.4f mm\n', De_inv*1e3);
fprintf('  BLC-corrected exit diameter      :  %.4f mm\n', De_eff*1e3);
fprintf('  Suggested De for next MOC run    :  %.4f mm\n', De_suggested*1e3);

fprintf('\n  BLC throat diameter    :  %.4f mm\n', 2*min(yb_blc)*1e3);
fprintf('  Total nozzle length    :  %.4f mm\n', (max(xb) - min(xb))*1e3);


%% All Plots

figure('Color', 'w');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, sprintf('Nozzle wall contour; $M_e = %.2f$; %s', Me, geom), 'Interpreter', 'latex');

% Tile 1: Extrapolated inviscid wall contour
nexttile;
plot(xb*1e3, yb*1e3, '.-', 'MarkerSize', 4);
xlabel('Axial distance [mm]');
ylabel('Radial distance [mm]');
title('Extrapolated inviscid wall contour');
xlim([xb(1)*1e3 - 100 xb(end)*1e3 + 100]);
ylim([yb(1)*1e3 - 30 yb(end)*1e3 + 25]);
axis equal;
grid on;

% Tile 2: Inviscid vs BLC wall contour
nexttile;
plot(xb*1e3, yb*1e3, '.-k', 'DisplayName', 'Inviscid wall', 'MarkerSize', 4);
hold on;
plot(xb_blc*1e3, yb_blc*1e3, '.-r', 'DisplayName', 'BLC wall', 'MarkerSize', 4);
xlabel('Axial distance [mm]');
ylabel('Radial distance [mm]');
title('Inviscid vs BLC wall');
legend('Location', 'southeast');
xlim([xb(1)*1e3 - 100 xb(end)*1e3 + 100]);
ylim([yb_blc(1)*1e3 - 30 yb_blc(end)*1e3 + 25]);
pbaspect([3 1 1]);
grid on;

% Tile 3: Displacement thickness with Mach distribution
nexttile;
yyaxis left;
plot(xb_div*1e3, delta_star*1e3, '-', 'MarkerSize', 4);
xlabel('Axial distance [mm]');
ylabel('$\delta^*$ [mm]');
yyaxis right;
plot(xb_div*1e3, M_div, '--');
ylabel('Mach number [-]');
title('Displacement thickness with Mach distribution');
xlim([xb_div(1)*1e3 - 100 xb_div(end)*1e3 + 100]);
pbaspect([3 1 1]);
grid on;

% Tile 4: CNC coordinate comparison (inviscid vs BLC)
nexttile;
plot(xi*1e3, yi*1e3, '.-k', 'MarkerSize', 3);
hold on;
plot(xv*1e3, yv*1e3, '.-r', 'MarkerSize', 3);
xlabel('Axial distance [mm]');
ylabel('Radial distance [mm]');
title('CNC coords: inviscid vs BLC');
legend('Inviscid', 'BL corrected', 'Location', 'southeast');
xlim([xi(1)*1e3 - 100 xi(end)*1e3 + 100]);
ylim([yv(1)*1e3 - 30 yv(end)*1e3 + 25]);
pbaspect([3 1 1]);
grid on;

%% Export
% All files are tab separated with columns x, y, z (z = 0).
% *_cnc_*  : uniform dx_cnc spacing, metres (prt 1, 2) or mm (prt 3)
% *_mesh_* : the assembled wall points, metres, for a mesh generator

out_file = @(suffix) fullfile(out_dir, sprintf('%s_%s.txt', base_name, suffix));

if prt == 1 || prt == 4
    % BL-corrected coordinates
    writematrix([xi_blc' yi_blc' zeros(size(xi_blc'))],  out_file('cnc_blc'),        'Delimiter', 'tab');
    writematrix([xb_blc  yb_blc  zeros(size(xb_blc))],   out_file('mesh_blc'),       'Delimiter', 'tab');
end
if prt == 2 || prt == 4
    % Inviscid coordinates
    writematrix([xi'  yi'  zeros(size(xi'))],            out_file('cnc_inviscid'),   'Delimiter', 'tab');
    writematrix([xb   yb   zeros(size(xb))],             out_file('mesh_inviscid'),  'Delimiter', 'tab');
end
if prt == 3 || prt == 4
    % BLC CNC coordinates in mm with a fixed number of decimals instead of
    % full double precision, and repeated points removed
    data_xy        = [xi_blc(:), yi_blc(:)] * 1e3;                 % convert to mm
    data_xy_r      = round(data_xy, cnc_decimals);
    data_xy_unique = unique(data_xy_r, 'rows', 'stable');
    data_final     = [data_xy_unique, zeros(size(data_xy_unique, 1), 1)];

    fmt    = sprintf('%%.%df\\t%%.%df\\t%%.%df\\n', cnc_decimals, cnc_decimals, cnc_decimals);
    fileID = fopen(out_file('cnc_blc_mm'), 'w');
    fprintf(fileID, fmt, data_final');
    fclose(fileID);
end
if prt > 0
    fprintf('\n  Coordinates written to %s\n', out_dir);
end

%% Helper functions - input

function [x, y, M] = read_moc_wall(fpath)

% READ_MOC_WALL   Wall contour from a MOC_Grid_BDE output file
%
%   [x, y, M] = read_moc_wall(fpath)
%
%   fpath may be the full summary.out (the "Calculated wall contour" block
%   is located and read), or a plain numeric table that holds only that
%   block: J, X/R*, R/R*, Mach, ... (header lines are allowed).
%
%   Outputs (column vectors):
%     x : X/R*,  y : R/R*,  M : wall Mach number

    txt   = fileread(fpath);
    lines = regexp(txt, '\r?\n', 'split');
    k     = find(contains(lines, 'Calculated wall contour'), 1, 'first');

    if ~isempty(k)
        data = [];
        for i = k+2 : numel(lines)                     % skip the column header line
            v = sscanf(lines{i}, '%f').';
            if isempty(v), break; end                  % blank line ends the block
            data(end+1, 1:numel(v)) = v;               %#ok<AGROW>
        end
    else
        a = importdata(fpath);
        if isstruct(a), a = a.data; end
        data = a;
    end

    if isempty(data) || size(data, 2) < 4
        error('Could not read a wall contour (J, X/R*, R/R*, Mach) from %s', fpath);
    end

    x = data(:,2);
    y = data(:,3);
    M = data(:,4);
end

%% Helper functions - BL correction

function delta_s = blc_edenfield_wall(xx, x0, pp, ma, tt, t0, g, Rgas)

% BLC_EDENFIELD_WALL   Edenfield (1968) turbulent BL displacement thickness
%
%   delta_s = blc_edenfield_wall(xx, x0, pp, ma, tt, t0, g, Rgas)
%
%   Compressible turbulent flat-plate BL with Eckert reference temperature
%   and Sutherland viscosity. Adiabatic wall assumed (Tw = T_recovery).
%
%   Inputs (all column vectors unless noted as scalar):
%     xx      : physical x-coords along wall          (m)
%     x0      : BL virtual origin                     (m)      scalar
%     pp      : static pressure along wall            (Pa)
%     ma      : Mach number along wall                (-)
%     tt      : static temperature along wall         (K)
%     t0      : stagnation temperature                (K)      scalar
%     g       : ratio of specific heats               (-)      scalar
%     Rgas    : gas constant                          (J/kg/K) scalar
%
%   Output:
%     delta_s : displacement thickness at each x      (m)
%
%   Ref: Edenfield E.E., "Contoured nozzle design and evaluation for
%        hotshot wind tunnels", AIAA Paper 68-369, 1968.

    Pr_t    = 0.9;                                     % turbulent Prandtl number
    r_rec   = Pr_t^(1/3);                              % recovery factor (turbulent)

    t_ad    = tt + r_rec .* (t0 - tt);                 % adiabatic recovery temperature (K)
    tw      = t_ad;                                    % adiabatic wall:  Tw = T_recovery

    t_ref   = 0.5*(tw + tt) + 0.22*(t_ad - tt);        % Eckert reference temperature (K)
    mu_ref  = sutherland_visc(t_ref);                  % viscosity at T_ref (Pa.s)
    rho_ref = pp ./ (Rgas .* t_ref);                   % reference density (kg/m^3)
    Ue      = ma .* sqrt(g .* Rgas .* tt);             % BL edge velocity (m/s)

    Re_ref  = rho_ref .* Ue .* (xx - x0) ./ mu_ref;    % reference Reynolds number

    delta_s = 0.42 .* (xx - x0) .* Re_ref.^(-0.2775);  % displacement thickness (m)

    delta_s = max(delta_s, 0);                         % guard against negative values
end


function mu = sutherland_visc(T)

%   mu = sutherland_visc(T)
%   T   : temperature [K]  (scalar or vector)
%   mu  : dynamic viscosity [Pa.s]

    mu0 = 1.716e-5;                                    % reference viscosity    (Pa.s)
    T0s = 273.11;                                      % reference temperature  (K)
    S   = 110.56;                                      % Sutherland constant    (K)

    mu  = mu0 .* (T./T0s).^1.5 .* (T0s + S) ./ (T + S);
end


%% Helper functions - geometry

function out = connect_line_circle_cubic(y_line, cx, cy, R, xb, yb, L, N)

% Build a cubic y(x) = a(x-xa)^3 + b(x-xa)^2 + y_line on [xa, xb] that:
%  - is tangent to horizontal line y=y_line at x=xa  (y'(xa)=0)
%  - hits the circle at (xb,yb) and is tangent to it there (y'(xb)=m_tan)
%  - has arc length exactly L (within tolerance)

    if nargin < 8 || isempty(N), N = 400; end
    if R <= 0 || L <= 0
        error('R and L must both be > 0.');
    end

    r_err = abs(hypot(xb-cx, yb-cy) - R);
    if r_err > 1e-6
        warning('Point (xb,yb) is off the circle by %.3g. Proceeding anyway.', r_err);
    end

    if abs(yb - cy) < 1e-12
        error(['Circle tangent is vertical at (xb,yb) -> infinite slope. ', ...
               'Use a rotated/parametric version instead.']);
    end

    function m = circle_tangent_slope_oriented(xa)
        r = [xb-cx, yb-cy];
        t = [-r(2), r(1)];
        if dot(t, [xa - xb, y_line - yb]) < 0
            t = -t;
        end
        m = t(2)/t(1);
    end

    function [a,b,mt,h] = coeffs_from_xa(xa)
        h = xb - xa;
        if abs(h) < 1e-12
            a = NaN; b = NaN; mt = NaN; return;
        end
        mt = circle_tangent_slope_oriented(xa);
        S1 = (yb - y_line) / (h^2);
        S2 = mt / h;
        a  = (S2 - 2*S1) / h;
        b  = 3*S1 - S2;
    end

    function s = arc_len(xa)
        [a,b,~,~] = coeffs_from_xa(xa);
        if ~isfinite(a) || ~isfinite(b)
            s = Inf; return;
        end
        f = @(x) sqrt(1 + (3*a*(x-xa).^2 + 2*b*(x-xa)).^2);
        s = integral(f, min(xa,xb), max(xa,xb), 'AbsTol',1e-9, 'RelTol',1e-9, 'ArrayValued',true);
    end

    F = @(xa) arc_len(xa) - L;

    B0 = max([2*L, 6*R, 1]);
    xa_sol = NaN; bracket = [];
    for grow = [1, 2, 4, 8]
        xL = xb - B0*grow;
        xR = xb + B0*grow;
        M  = 60;
        xs = linspace(xL, xR, M);
        Fs = arrayfun(F, xs);
        idx = find(sign(Fs(1:end-1)) .* sign(Fs(2:end)) <= 0 & isfinite(Fs(1:end-1)) & isfinite(Fs(2:end)), 1, 'first');
        if ~isempty(idx)
            bracket = [xs(idx), xs(idx+1)];
            break;
        end
    end

    if ~isempty(bracket)
        xa_sol = fzero(F, bracket);
    else
        xa_sol = fminbnd(@(x) F(x).^2, xb - 8*B0, xb + 8*B0);
    end

    [a,b,mt,h] = coeffs_from_xa(xa_sol); %#ok<ASGLU>
    achievedL  = arc_len(xa_sol);
    if abs(achievedL - L) > 1e-6
        warning('Could not hit length exactly. Achieved %.9g vs target %.9g (diff %.3g).', achievedL, L, abs(achievedL-L));
    end

    xs = linspace(min(xa_sol, xb), max(xa_sol, xb), N).';
    u  = xs - xa_sol;
    ys = a*u.^3 + b*u.^2 + y_line;

    out.xy     = [xs, ys];
    out.xa     = xa_sol;
    out.coeffs = [a, b];
    out.m_tan  = mt;
    out.length = achievedL;
end


function wall_points = conv_cont2(D1, Dt, L, mL, N)

% Cubic convergent profile with y'(0) = 0 and y'(L) = mL
% y(0) = D1/2, y(L) = Dt/2

    if nargin < 4 || isempty(mL), mL = 0; end
    if nargin < 5 || isempty(N),   N   = 1000; end

    x = linspace(0, L, N).';

    y0    = D1/2;
    yL    = Dt/2;
    delta = yL - y0;

    a =  (mL / L^2) - (2*delta / L^3);
    b = -(mL / L)   + (3*delta / L^2);

    y = a*x.^3 + b*x.^2 + y0;

    wall_points = [x, y];
end
