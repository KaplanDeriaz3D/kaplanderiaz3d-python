# KaplanDeriaz3D — User Guide

**KaplanDeriaz3D** is an open-source parametric design tool for the 3D blade
geometry of Kaplan (axial-flow) and Deriaz (diagonal-flow) hydraulic
turbines. Given a small set of global hydraulic parameters (net head, flow
rate, rotational speed, geometry bounds), the application derives the blade
angle distribution directly from the turbomachinery Euler equation and
integrates it into a full 3D solid, ready for CAD/CFD or 3D printing.

This guide goes deliberately deeper than a typical tooltip reference: every
formula the application actually evaluates is written out explicitly here,
so you can understand exactly what the tool computes without having to
read the source code. It documents the **base version** of the tool: the
validated hydrodynamic core (infinite-blade, zero-thickness Euler model),
with no industrial post-processing adjustments.

---

## Table of Contents

1. [Interface Overview](#1-interface-overview)
2. [General Tab — Design Parameters](#2-general-tab--design-parameters)
3. [The Interpolation Laws, Explicitly](#3-the-interpolation-laws-explicitly)
4. [The Hydrodynamic Model, Step by Step](#4-the-hydrodynamic-model-step-by-step)
5. [Hydrofoil / Mesh Tab](#5-hydrofoil--mesh-tab)
6. [The Hydrofoil Profile Formulas, Explicitly](#6-the-hydrofoil-profile-formulas-explicitly)
   - [6.1 What "thickness direction" actually is](#61-what-thickness-direction-actually-is)
7. [Export Tab](#7-export-tab)
8. [3D Visualization Panel](#8-3d-visualization-panel)
9. [Results Table, With Every Formula](#9-results-table-with-every-formula)
10. [Menu Bar](#10-menu-bar)
11. [Known Limitations and Scope](#11-known-limitations-and-scope)

---

## 1. Interface Overview

The window is split into two halves:

- **Left panel** — three tabs (`General`, `Hydrofoil / mesh`, `Export`)
  containing every input parameter, organized by purpose.
- **Right panel** — the 3D viewer (two side-by-side viewports: the mean
  surface and the solid blade) plus the results table underneath.

A menu bar (`File`, `View`, `Help`) at the top gives access to
configuration save/load, view controls, and theme switching.

Every input field's **label** (not the field itself) carries a tooltip —
hover over any label to see a short explanation. This guide gives the full
mathematical picture behind those tooltips.

The workflow is always the same, regardless of turbine type:

1. Set your parameters in the `General` tab.
2. Set your hydrofoil profile in the `Hydrofoil / mesh` tab.
3. Press **COMPUTE SOLID 3D DESIGN**.
4. Inspect the 3D render and the results table.
5. Export via the `Export` tab or the `File` menu.

---

## 2. General Tab — Design Parameters

### Turbine Type

A dropdown to choose between **Kaplan (Axial)** and **Deriaz (Diagonal)**.
Switching this resets the hydraulic parameters below to sensible defaults
for that turbine type, and swaps the geometry group (Kaplan Geometry vs.
Deriaz Geometry) accordingly.

### Design Parameters group

| Field | Symbol | Default | Range | Meaning |
|---|---|---|---|---|
| Rotational Speed (RPM) | `RPM` | 450 (Kaplan) / 150 (Deriaz) | 1 – 10 000 | Shaft speed of the runner. |
| Rotation Direction | — | Counter-Clockwise | CCW / CW | Sign convention for the blade twist integration (`rot_sign` = +1 or −1). |
| Flow Rate Q0 (m³/s) | `Q0` | 12 / 110 | 0.01 – 10 000 | Nominal design discharge. |
| Net Head Hn (m) | `Hn` | 15 / 60 | 0.1 – 2000 | Available net head across the machine. |
| Gravity g (m/s²) | `g` | 9.81 | 1 – 20 | Local gravitational acceleration. |
| **Water density rho (kg/m³)** | `rho` | **1000** | 1 – 10 000 | Fluid density, used **only** for the power/torque figures in the results table (Section 9) — it plays no role anywhere in the blade-angle or geometry derivation itself, which is why it doesn't appear in the equations of Section 4. |
| Hydraulic Eff. (eta_h) | `eta_h` | 0.90 / 0.91 | 0.1 – 1.0 | Hydraulic efficiency; scales the head into the "infinite head" the Euler equation actually uses (see Section 4). |
| Volumetric Eff. (eta_v) | `eta_v` | 0.96 / 0.97 | 0.1 – 1.0 | Volumetric efficiency; scales the flow rate into the "real flow" used in the velocity triangles. |
| Mechanical Eff. (eta_m) | `eta_o` | 0.98 | 0.1 – 1.0 | Mechanical/organic efficiency; used only for the net shaft power output (Section 9), not for the blade geometry. |
| Solidity (sigma) | `sigma_target` | 1.25 / 1.40 | 0.5 – 3.0 | Target chord-to-pitch ratio, used to determine the optimal number of blades (see Section 4). |
| Interpolation Scheme | — | Cubic (Standard) | 5 options | The law used to blend the blade angle between its inlet and outlet values along the chord — see Section 3 for every formula. |

> **Why isn't water density used in the physics?** The Euler-equation
> derivation that produces the blade geometry (Section 4) works entirely
> in terms of velocities and angles — density cancels out of it
> completely, the same way it does in the classical derivation of the
> Euler turbomachinery equation. It only re-enters the picture when the
> tool converts hydraulic power into a physical force/torque figure for
> the results table (`P_hyd = rho * g * Q0 * Hn`), which is why `rho` sits
> in the General tab (it's a genuine design input) but never appears
> inside `compute_kaplan_backbone()` or `compute_deriaz_backbone()`.

### Kaplan Geometry (visible only when Kaplan is selected)

| Field | Symbol | Default | Range | Meaning |
|---|---|---|---|---|
| Hub Radius R_hub (m) | `R_hub` | 0.30 | 0.01 – 50 | Radius of the cylindrical hub surface. |
| Tip Radius R_tip (m) | `R_tip` | 0.65 | 0.01 – 50 | Radius of the cylindrical shroud/casing surface. |
| Axial Length L_z (m) | `L_z` | 0.35 | 0.01 – 20 | Axial extent of the blade, from leading edge (`z=0`) to trailing edge (`z=-L_z`). |

Geometrically, every streamline in a Kaplan runner lives on its own
**cylinder** of constant radius `r` between `R_hub` and `R_tip` — this is
the literal meaning of "axial flow": the radius never changes as the flow
moves through the runner, only the axial position `z` and the
circumferential angle `theta` do.

### Deriaz Geometry (visible only when Deriaz is selected)

| Field | Symbol | Default | Range | Meaning |
|---|---|---|---|---|
| Inner Spherical Radius Re_int (m) | `Re_int` | 2.00 | 0.01 – 100 | Radius of the spherical hub surface. |
| Outer Spherical Radius Re_ext (m) | `Re_ext` | 2.60 | 0.01 – 100 | Radius of the spherical shroud/casing surface. |
| Inlet Cone Angle gamma1 (deg) | `gamma1` | 30.0 | 0 – 89 | Spherical angle at the leading edge. |
| Outlet Cone Angle gamma2 (deg) | `gamma2` | 60.0 | 0 – 89 | Spherical angle at the trailing edge. |

**What gamma actually is, geometrically.** Every Deriaz streamline lives
on its own **sphere** of constant radius `Re` between `Re_int` and
`Re_ext`, centred at the origin. A point's position on that sphere is
parametrized by the angle `gamma`, measured from the **equatorial plane**
(`z = 0`) toward the downstream pole, exactly like a latitude angle:

```
rc(gamma) = Re * cos(gamma)      <- projected cylindrical radius
z(gamma)  = -Re * sin(gamma)     <- axial position (negative = downstream)
```

At `gamma = 0` the point sits on the equator (`z = 0`, maximum
cylindrical radius `rc = Re`); as `gamma` increases toward 90°, the point
moves toward the pole (`z` more negative — further downstream — and `rc`
shrinking toward 0). Because `gamma1 < gamma2` by construction (the
leading edge sits closer to the equator, the trailing edge closer to the
pole), water genuinely travels diagonally — partly radial, partly axial —
which is exactly what gives the Deriaz turbine its name and its
double-regulation capability. `rc(gamma)` is the direct Deriaz analogue of
the constant radius `r` in the Kaplan case; it's what actually enters the
peripheral-speed and blade-angle formulas in Section 4.

---

## 3. The Interpolation Laws, Explicitly

The blade angle `beta` is only computed directly at the leading edge
(`beta1`) and the trailing edge (`beta2`) — see Section 4 for those
formulas. Everywhere in between, `beta` is obtained by **blending**
`beta1` and `beta2` through one of five interpolation laws, using a
normalized chordwise coordinate `t` that runs from `t=0` at the leading
edge to `t=1` at the trailing edge:

```
beta(t) = beta1 + (beta2 - beta1) * f(t)
```

where `f(t)` is one of the five functions below, all satisfying
`f(0) = 0` and `f(1) = 1` exactly — so `beta` always equals `beta1` at the
inlet and `beta2` at the outlet regardless of which law you pick; only the
*shape* of the transition between them changes, which changes the assumed
blade loading (pressure) distribution.

| Interpolation Scheme | Formula `f(t)` | Loading behaviour |
|---|---|---|
| **Cubic (Standard)** | `f(t) = 3t² − 2t³` | Smooth S-curve, zero slope at both ends (a Hermite smoothstep). Peak loading around 30–40% of the chord. **Default.** |
| **Linear (Uniform)** | `f(t) = t` | Constant angle gradient; uniform loading along the whole profile. |
| **Cosine (Smooth)** | `f(t) = 0.5 · (1 − cos(π·t))` | Ultra-smooth transition near both edges, minimizing localized cavitation spikes. |
| **Inlet Loaded (Attack)** | `f(t) = 1 − (1 − t)²` | Steeper deflection near the leading edge; maximum loading at the inlet. |
| **Outlet Loaded (Discharge)** | `f(t) = t²` | Steeper deflection near the trailing edge; maximum loading at the outlet. |

For **Kaplan**, `t` is the normalized axial position: `t = |z| / L_z`.
For **Deriaz**, `t` is the normalized spherical angle:
`t = (gamma − gamma1) / (gamma2 − gamma1)`. Either way, `t` sweeps from 0
to 1 across the chord, and the same five formulas above apply unchanged.

---

## 4. The Hydrodynamic Model, Step by Step

Both turbine types are derived from the **Euler turbomachinery equation**,
under the classical idealization of an **infinite number of
infinitely thin blades** (no slip, no blockage — see Section 11 for what
this leaves out).

### 4.1 Quantities shared by both turbine types

```
eta_t  = eta_h * eta_v * eta_o                  (overall efficiency)
omega  = 2*pi*RPM / 60                          (angular speed, rad/s)
Q_real = Q0 * eta_v                             (effective flow rate)
H_inf  = Hn * eta_h                             (effective/"infinite" head)
nq     = RPM * sqrt(Q0) / (g*Hn)^0.75           (specific speed — note g in the denominator)
```

`nq` is the dimensionless specific speed of the design point — this is
the one figure in the whole model where gravitational acceleration
appears in the **denominator**, raised to the 0.75 power alongside the
head; it characterizes the "shape" of turbine the design point calls for
(low `nq` → Pelton-like; high `nq` → propeller-like), independent of the
machine's absolute size.

### 4.2 Kaplan (axial) — `compute_kaplan_backbone()`

Geometric setup:

```
R_m       = sqrt((R_tip² + R_hub²) / 2)         (mean hydraulic radius)
z_vec     = linspace(0, -L_z, N_cuerda)         (chordwise axial stations)
r_vec     = linspace(R_hub, R_tip, N_radios)    (spanwise radial stations)
Area_paso = pi * (R_tip² - R_hub²)              (annular flow area)
V_z       = Q_real / Area_paso                  (uniform meridional/axial velocity)
```

At each radial station `r`:

```
U        = omega * r                            (peripheral/blade speed)
V_theta1 = g * H_inf / (omega * r)               (required inlet tangential velocity)
```

This comes directly from the Euler equation with the standard
zero-exit-swirl design condition (`V_theta2 = 0`, i.e. the flow leaves the
runner with no residual rotation) and the infinite-blade assumption. If
`V_theta1 >= U` at any radius, the velocity triangle is physically
impossible (the blade would have to push the flow backward relative to
its own motion) — the tool reports this as an error rather than producing
invalid geometry.

The inlet and outlet blade angles (measured from the tangential
direction) follow directly from the velocity triangle:

```
beta1 = atan2(V_z, U - V_theta1)
beta2 = atan2(V_z, U)                           (V_theta2 = 0 by design)
```

The chordwise `beta(z)` distribution is then obtained by blending `beta1`
and `beta2` through the selected interpolation law (Section 3), and the
backbone is built by integrating the local rate of circumferential twist
station by station along `z`:

```
d(theta)/dz = -rot_sign * cot(beta(z)) / r
```

Each streamline lives on its own cylinder of constant radius `r` — only
`theta` (not `r`) evolves along it, which is the geometric meaning of
"axial flow."

Finally:

```
L_chord_reference = arc length of the streamline at r = R_m
Z_optimo = clamp( round( (2*pi*R_m) / (L_chord_reference / sigma_target) ), 3, 12 )
```

`Z_optimo`, the blade count, comes from a **target solidity** (chord-to-
pitch ratio) evaluated at the mean radius: the pitch (circumferential
spacing between blades) is `2*pi*R_m / Z`, so requiring
`chord / pitch = sigma_target` gives the formula above directly, clamped
to a structurally/hydraulically sensible range of 3–12 blades.

### 4.3 Deriaz (diagonal) — `compute_deriaz_backbone()`

Same underlying physics, adapted to the spherical geometry described in
Section 2. The mean radius used for the reference chord and blade count
is the **log-mean** spherical radius (the natural mean for a domain whose
effective circumference scales with `cos(gamma)*Re`, not a constant `r`):

```
Re_medio = (Re_ext - Re_int) / ln(Re_ext / Re_int)
```

At each spherical radius `Re`, the inlet and outlet stations are
evaluated separately (since the projected radius differs between them):

```
rc1 = Re * cos(gamma1);  U1 = omega * rc1
Vm1 = Q_real / (2*pi * Re*cos(gamma1) * (Re_ext - Re_int))     (inlet meridional velocity)
Vtheta1 = g * H_inf / U1

rc2 = Re * cos(gamma2);  U2 = omega * rc2
Vm2 = Q_real / (2*pi * Re*cos(gamma2) * (Re_ext - Re_int))     (outlet meridional velocity)
```

Same physical-feasibility check as Kaplan (`Vtheta1 >= U1` → error), then:

```
beta1 = atan2(Vm1, U1 - Vtheta1)
beta2 = atan2(Vm2, U2)                          (V_theta2 = 0 by design)
```

The chordwise `beta(gamma)` distribution is obtained the same way as
Kaplan (Section 3, using `t = (gamma-gamma1)/(gamma2-gamma1)`), and the
backbone is integrated station by station along `gamma` instead of `z`:

```
d(theta)/d(gamma) = -rot_sign * cot(beta(gamma)) / cos(gamma)
```

— the direct spherical analogue of the Kaplan cylindrical integration.
Finally, `L_chord_reference` and `Z_optimo` are computed the same way as
Kaplan (Section 4.2), evaluated at the `Re_medio` station, with the blade
count clamped to 4–12 blades instead of 3–12 (Deriaz runners structurally
require a minimum of 4 blades for the pivot/closing mechanism).

---

## 5. Hydrofoil / Mesh Tab

### Hydrofoil Profile Design group

**Hydrofoil Profile** — five available families (formulas in Section 6):

- **NACA 00XX (Standard Symmetric)** — classic symmetric NACA thickness
  profile.
- **Customized (4-Digit NACA Series)** — lets you specify camber (`m`) and
  its chordwise position (`p`) independently at the hub and the tip; a new
  sub-panel appears below the dropdown when this is selected.
- **Reversible Hydrofoil (Pump-Turbine)** — symmetric elliptical thickness,
  no camber, suited to bidirectional (pump/turbine) flow.
- **Anti-Cavitation (Flat Pressure)** — forward-loaded thickness profile
  aimed at flattening the pressure-side distribution.
- **Low-Torque S-Camber** — S-shaped camber line intended to reduce the net
  actuation torque on the blade pitch mechanism.

**Max Rel. Thickness (t/c)** (`self.spin_thickness`, default **0.085**,
range 0.01–0.40) — the reference relative thickness, measured at the hub
chord.

**Hub to Tip Thickness Ratio** (`self.spin_hub_to_tip`, default **0.65**,
range 0.10–1.00) — how much thinner the blade gets at the tip relative to
the hub (0.65 means the tip retains 65% of the hub's relative thickness).

**Customized sub-panel** (only when "Customized" is selected): `m_hub`,
`p_hub`, `m_tip`, `p_tip` — see Section 6 for exactly how these four
numbers become the camber line.

### Advanced Options (Mesh Resolution)

| Field | Default | Range | Meaning |
|---|---|---|---|
| Streamlines (N_radii) | 100 | 5 – 200 | Number of radial evaluation stations from hub to tip. More = smoother spanwise geometry, slower computation. |
| Chord Stations (dz/dgamma/ds) | 200 | 15 – 500 | Number of chordwise points per streamline. More = smoother chordwise geometry, slower computation. |

---

## 6. The Hydrofoil Profile Formulas, Explicitly

Before the profile formulas themselves, one conversion step matters and
is easy to miss: **the "Max Rel. Thickness (t/c)" you set in the GUI is
not used directly, station by station.** `build_solid_blade()` first
fixes an *absolute* reference thickness from the HUB chord alone, then
re-derives a *local* relative thickness for every other station from
that fixed absolute value and that station's own chord length, and only
*then* applies the hub-to-tip taper you actually control:

```
t_abs_ref      = t_rel_input * c_hub          (fixed once, from the hub chord)
t_rel_local(i) = t_abs_ref / c_i              (re-relativized to station i's own chord)
```

The effect: if two stations happened to have the *same* chord length,
`t_rel_local` would be identical between them regardless of any chord
differences elsewhere on the blade — the absolute (metric) thickness
tracks the chord automatically, station by station, rather than the
input `t/c` being applied as a flat percentage everywhere. `t_rel_local`
is what actually gets passed into `generate_hydro_profile()` as `t_c_base`
below — **not** the raw GUI input.

`generate_hydro_profile()` then returns a camber line `yc(x)` and a
half-thickness distribution `yt(x)`, as functions of the normalized
chordwise position `x` in [0, 1], for one radial station. This is where
the hub-to-tip taper you set in the GUI is actually applied, on top of
the re-relativization above:

```
r_star   = i_radio / (N_radios - 1)                           (0 at hub, 1 at tip)
t_c      = t_rel_local * (1 - (1 - hub_to_tip_ratio)*r_star)   (your taper, applied here)
m_factor = (1 - r_star)^1.5                                    (camber taper toward the tip)
```

The standard 4-digit NACA thickness formula, used by three of the five
families below:

```
yt_naca(x) = (t_c/0.2) * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x² + 0.2843*x³ - 0.1015*x⁴)
```

| Profile | Camber `yc(x)` | Thickness `yt(x)` |
|---|---|---|
| **NACA 00XX** | `0` | `yt_naca(x)` |
| **Reversible Hydrofoil** | `0` | `t_c * sqrt(x*(1-x))` |
| **Anti-Cavitation** | `yt(x) / 2 * m_factor` | `2.6896 * t_c * x * (1-x)^1.5` |
| **Low-Torque S-Camber** | `0.8 * x*(1-x)*(0.5-x) * m_factor` | `yt_naca(x)` |
| **Customized** | see below | `yt_naca(x)` |

**Customized camber line** — the classic 4-digit NACA camber formula,
with `m` (max camber, %) and `p` (its chordwise position, tenths of
chord) linearly blended across the span from the tip values to the hub
values using `m_factor`:

```
m = m_tip/100 + (m_hub/100 - m_tip/100) * m_factor
p = max( p_tip/10 + (p_hub/10 - p_tip/10) * m_factor,  0.05 )

for x < p:  yc(x) = (m/p²)   * (2*p*x - x²)
for x >= p: yc(x) = (m/(1-p)²) * ((1-2p) + 2*p*x - x²)
```

Once `yc` and `yt` are known at a station, they still need to be placed
in 3D. `build_solid_blade()` does this by offsetting the extrados
(pressure side) and intrados (suction side) surfaces away from the mean
surface **along a local thickness direction** — the vector math behind
that direction is worth spelling out explicitly too, because it is the
one part of the geometry pipeline most likely to look "wrong" if you
change it without understanding why it's built this way.

### 6.1 What "thickness direction" actually is

The naive approach would be to offset along the mean surface's own full
normal vector, `N_naive = normalize(Tu × Tv)`, where `Tu` and `Tv` are the
chordwise and spanwise tangent directions (estimated numerically with
`np.gradient()`). **This is deliberately NOT what the tool does**,
because that full normal mixes the chordwise curvature of the profile
with the blade's spanwise twist rate — offsetting along it makes a thick
blade visibly "lean sideways" as thickness increases, since the twist
component grows right along with the thickness offset.

Instead, `compute_thickness_direction()` builds a direction confined to
the local **chord/circumferential plane**, the same way a 2D airfoil
section is classically defined before being stacked along the span:

```
Tu_hat  = normalize(Tu)                              (unit chordwise tangent)

rc      = sqrt(X² + Y²)                              (projected cylindrical radius)
e_theta = (-Y/rc, X/rc, 0)                            (local circumferential unit vector)
```

`e_theta` is the tangent to the circle of radius `rc` centred on the Z
axis at that point — valid for both Kaplan (cylinders) and Deriaz
(spheres), since both are surfaces of revolution about the Z axis. The
thickness direction is then `e_theta`, **projected to be exactly
perpendicular to the chordwise tangent** (a single step of Gram-Schmidt
orthogonalization):

```
N = e_theta - (e_theta · Tu_hat) * Tu_hat
N_thickness = normalize(N)
```

**Fallback for the degenerate case.** If `e_theta` happens to be nearly
parallel to `Tu_hat` at some point (so the projection above nearly
cancels out, `|N| ≈ 0`), the classic full surface normal
`normalize(Tu × Tv)` is used at that point instead — this only matters at
isolated, unusual points in the geometry, not across the blade generally.

`build_solid_blade()` then places the two surfaces as:

```
extrados = mean_surface + (yc·c_i + yt·c_i) * N_thickness
intrados = mean_surface + (yc·c_i - yt·c_i) * N_thickness
```

where `c_i` is that station's own chord length (camber shifts both
surfaces together along `N_thickness`; thickness pushes them apart along
the same direction).

> **A note on transparency:** the source file also contains a
> `compute_surface_normals()` function implementing exactly the naive
> `normalize(Tu × Tv)` approach described above — it is **not called
> anywhere** in the current pipeline; it's kept in the file purely as a
> documented reference for why the fix in `compute_thickness_direction()`
> was needed, and is mentioned here so nothing in the source file is left
> unexplained.

---

## 7. Export Tab

### Export Manager group

- **Export Geometry** — `Solid Blade (Full)` (extrados + intrados + end
  caps, a watertight closed mesh) or `Mean Surface Only` (a single
  zero-thickness sheet).
- **CAD Format** — `STL` (triangulated mesh, universally compatible) or
  `IGES` (a point-cloud fallback writer — **not** a true B-Rep solid; for a
  real IGES/STEP solid, pass the exported STL through a CAD kernel such as
  FreeCAD or OpenCASCADE).
- **Export CAD Geometry** button — opens a save dialog and writes the file.
  The "Solid Blade (Full)" STL is also the geometry to use for any CFD
  meshing work — it's watertight and ready to mesh with no further
  preprocessing.
- **Export Blade Angles (Excel)** button — writes a two-column `.xlsx`
  table of inlet/outlet blade angle (`beta1`, `beta2`, in degrees) vs.
  radius, for feeding into external tools.

Both export buttons are disabled until a geometry has been computed at
least once.

---

## 8. 3D Visualization Panel

Two side-by-side PyVista viewports:

- **Left — Mean Surface Design.** The zero-thickness backbone, colour-mapped
  by radius (Kaplan) or projected cylindrical radius (Deriaz).
- **Right — Solid Runner.** The full thickened blade(s), same colour
  mapping.

**Reset 3D View** button (top-right) restores the default isometric camera
angle on both viewports. Camera position/zoom is otherwise preserved
between recomputes of the *same* turbine type, so you don't lose your
viewing angle every time you tweak a parameter and hit Compute again — it
only resets automatically when you switch between Kaplan and Deriaz (since
the two geometries live at very different scales).

---

## 9. Results Table, With Every Formula

Populated after every successful computation:

| Row | Formula |
|---|---|
| Airfoil Selection | (the selected hydrofoil profile family, text only) |
| Total Efficiency (eta_t) | `eta_t = eta_h * eta_v * eta_o` |
| Specific Speed (nq) | `nq = RPM * sqrt(Q0) / (g*Hn)^0.75` |
| Integrated 3D Chord (L_chord_ref) | Arc length of the mean-radius streamline (Section 4) |
| Number of Blades (Z) | `Z_optimo` (Section 4) |
| Mean Hydraulic/Spherical Radius | `R_m = sqrt((R_tip²+R_hub²)/2)` (Kaplan) or `Re_medio = (Re_ext-Re_int)/ln(Re_ext/Re_int)` (Deriaz) |
| Gross Hydraulic Power (P_hyd) | `P_hyd = rho * g * Q0 * Hn` |
| Net Mechanical Power (P_mech) | `P_mech = P_hyd * eta_t` |
| Shaft Torque (T) | `T = P_mech / omega` |
| Angular Speed (omega) | `omega = 2*pi*RPM / 60` |

Note that `P_hyd`, `P_mech`, and `T` are the only three figures in the
entire results table that use the *gross* design values (`Q0`, `Hn`)
together with density `rho` — every other row, and the blade geometry
itself, is derived from the *effective* values (`Q_real`, `H_inf`) defined
in Section 4.1.

---

## 10. Menu Bar

### File

- **Save Configuration...** — writes every input parameter to a JSON file.
- **Load Configuration...** — restores a previously saved parameter set.
- **Reset to Defaults** — restores the pristine Kaplan-default state (as if
  the app had just started), and clears any previous computation results
  so the reset is visually unambiguous.
- **Export CAD Geometry...** / **Export Blade Angles (Excel)...** —
  shortcuts to the same exports available in the Export tab.
- **Close**

### View

- **Reset 3D View** — same as the button in the visualization panel.
- **Copy Results Table** — copies the results table to the clipboard as
  tab-separated text.
- **Light Mode** / **Dark Mode** — theme switch (dark is the default).

### Help

- **About KaplanDeriaz3D...** — version and licensing information.

---

## 11. Known Limitations and Scope

Documented here explicitly, in the same spirit of transparency as the
tooltips throughout the app (see also the project's `CONTRIBUTING_ROADMAP`
document, which turns every one of these into a concrete, actionable
contribution opportunity):

- **Infinite-blade, zero-thickness assumption.** The core hydrodynamic
  model (Section 4) does not account for finite-blade slip or
  blade-thickness flow blockage. This is the standard textbook
  simplification.
- **No centrifugal contribution to the meridional velocity.** The
  meridional/axial velocity (`V_z` for Kaplan, `Vm1`/`Vm2` for Deriaz) is
  derived purely from mass-flow continuity (`Q_real / Area`), assumed
  uniform across the annulus. The model does **not** account for the
  contribution of centrifugal force to the meridional velocity
  distribution — a real effect (particularly relevant in Deriaz pump
  mode, where the diagonal flow path means the fluid genuinely
  accelerates radially, not just axially) that a radial-equilibrium
  treatment would capture but this streamline-by-streamline model does
  not.
- **No CFD validation.** Geometry produced by this tool is a preliminary,
  parametric design output — it has not been validated against
  Navier-Stokes CFD simulation, and does not predict efficiency, pressure
  distribution, or cavitation risk on its own.
- **No structural analysis.** Blade thickness is set by the hydrofoil
  profile selection alone; no stress, fatigue, or mechanical-actuation
  analysis is performed.
- **No industrial post-processing.** This version outputs the geometry
  exactly as derived from the inverse hydrodynamic design — it does not
  include any finishing operations (edge sweep, trailing-edge alignment,
  corner rounding, wake trimming) that a manufactured runner would
  typically go through. Applying such operations, if needed, is left to
  downstream CAD work outside this tool.
- **Francis turbines are not supported.** Only Kaplan and Deriaz geometry
  generation is implemented.

---

*This guide documents the KaplanDeriaz3D open-source project. For the
full mathematical derivations and the specific-speed/solidity theory
behind the formulas above, see the accompanying thesis document. For a
detailed roadmap of exactly what's missing and how to contribute, see
`CONTRIBUTING_ROADMAP.md` in the same repository.*
