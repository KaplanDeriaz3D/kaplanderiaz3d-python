# KaplanDeriaz3D - Hydraulic Turbine Blade Designer (Professional Edition)
# Copyright (C) 2026 Juan Fernandez Lozano
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


"""
KaplanDeriaz3D - Hydraulic Turbine Blade Designer (Professional Edition)
=========================================================================
Python/PySide6 + PyVista port of KaplanDeriaz3D_Airfoil_Version_Professional.m (MATLAB App Designer).

Requirements:
    pip install PySide6 pyvista pyvistaqt numpy openpyxl

Compatibility note: the code is written against PySide6, but the API used
(QWidget, QGridLayout, QDoubleSpinBox, Signal/Slot, etc.) is nearly
identical in PyQt6. To port to PyQt6:
    - change 'from PySide6...' imports to 'from PyQt6...'
    - change 'Signal' to 'pyqtSignal' (no custom Signals are used here)
    - 'app.exec()' stays the same in both (PyQt6 also uses exec())

Run:
    python kaplan_deriaz_python_app_Professional.py
"""

import sys
import json
import datetime
import numpy as np

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFrame, QVBoxLayout, QHBoxLayout,
    QGridLayout, QGroupBox, QLabel, QComboBox, QDoubleSpinBox,
    QSpinBox, QPushButton, QTableWidget, QTableWidgetItem, QCheckBox,
    QMessageBox, QFileDialog, QHeaderView, QScrollArea, QTabWidget,
    QStatusBar,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QAction, QKeySequence

import pyvista as pv
from pyvistaqt import QtInteractor

import openpyxl

# =============================================================================
# MATHEMATICAL CORE (no Qt dependencies)
# =============================================================================

def evaluate_interpolation(interp_type, t):
    """Evaluates the blade-angle interpolation LAW: given a normalized
    chordwise (or gamma-wise, for Deriaz) coordinate t in [0, 1], returns a
    monotonically increasing weight in [0, 1] used to blend the local blade
    angle beta between its inlet value (beta1, at t=0) and its outlet value
    (beta2, at t=1): beta(t) = beta1 + (beta2 - beta1) * evaluate_interpolation(t).

    Five interpolation families are offered, each corresponding to a
    different assumed blade loading (pressure) distribution along the
    chord:
      - Cubic (Standard):  3t^2 - 2t^3.   Smooth S-curve, zero slope at both
        ends (Hermite smoothstep). Peak loading around 30-40% chord.
      - Linear (Uniform):  t.             Constant angle gradient; uniform
        loading along the whole profile.
      - Cosine (Smooth):   0.5*(1-cos(pi*t)). Ultra-smooth transition near
        both edges, minimizing localized cavitation spikes.
      - Inlet Loaded:      1-(1-t)^2.     Steeper deflection near the inlet;
        maximum loading at the leading edge.
      - Outlet Loaded:     t^2.           Steeper deflection near the
        outlet; maximum loading at the trailing edge.
    All five satisfy the two boundary conditions f(0)=0, f(1)=1 exactly, so
    beta always equals beta1 at the leading edge and beta2 at the trailing
    edge regardless of which law is selected; the choice of law only
    changes how beta transitions between those two fixed endpoints."""
    t = np.asarray(t, dtype=float)
    if interp_type == 'Cubic (Standard)':
        return 3 * (t ** 2) - 2 * (t ** 3)
    elif interp_type == 'Linear (Uniform)':
        return t
    elif interp_type == 'Cosine (Smooth)':
        return 0.5 * (1 - np.cos(np.pi * t))
    elif interp_type == 'Inlet Loaded (Attack)':
        return 1 - (1 - t) ** 2
    elif interp_type == 'Outlet Loaded (Discharge)':
        return t ** 2
    else:
        return 3 * (t ** 2) - 2 * (t ** 3)


def compute_kaplan_backbone(RPM, Q0, Hn, g, eta_h, eta_v, eta_o, sigma_target,
                             interp_type, rot_sign, N_radios, N_cuerda,
                             R_hub, R_tip, L_z):
    """SIMPLIFIED axial (Kaplan) hydrodynamic model: derives and integrates
    the blade mean-surface backbone directly from the turbomachinery Euler
    equation, assuming an infinite number of zero-thickness blades (the
    classic textbook idealization; see the thesis Section 6.4.2 discussion
    for the finite-blade correction this deliberately omits).

    PHYSICS, step by step, at each radial station r in [R_hub, R_tip]:
      1. Peripheral (blade) speed: U = omega * r.
      2. Required inlet tangential velocity, from the Euler equation with
         zero exit swirl (V_theta2 = 0) and infinite blades:
             V_theta1 = g * H_inf / (omega * r)
         (H_inf = Hn * eta_h is the hydraulic-efficiency-corrected head).
         If V_theta1 would reach or exceed U, the velocity triangle is
         physically impossible (the blade would need to turn the flow
         backward relative to its own motion) - this is reported as an
         error rather than silently producing an invalid geometry.
      3. Meridional (through-flow) velocity from continuity, assumed
         uniform across the annulus: V_z = Q_real / Area_paso.
      4. Inlet/outlet blade angles beta1, beta2 (measured from the
         tangential direction) from the velocity triangle:
             beta1 = atan2(V_z, U - V_theta1)
             beta2 = atan2(V_z, U)             (V_theta2 = 0 by design)
      5. The chordwise beta(z) distribution is obtained by blending beta1
         and beta2 through evaluate_interpolation() using the user-selected
         law (cubic, linear, cosine, inlet- or outlet-loaded).
      6. The backbone itself is built by integrating the blade differential
         equation station by station along z, using the geometric relation
         between the blade angle and the local rate of circumferential
         twist: d(theta)/dz = -cot(beta(z)) / r. Each streamline lives on
         its own cylinder of constant radius r (this is what "axial flow"
         means geometrically), so only theta (not r) evolves along it.

    Returns a dict with the full backbone (X_mid, Y_mid, Z_mid on a
    (N_cuerda x N_radios) grid), the inlet/outlet angle vectors, the
    solidity-derived optimal blade count Z_optimo, and the summary
    hydraulic quantities used elsewhere in the app. On any physical
    infeasibility, returns {'error': <message>} instead."""
    if R_hub >= R_tip:
        return {'error': 'Geometric bounds error: R_hub must be smaller than R_tip.'}

    eta_t = eta_h * eta_v * eta_o
    omega = (2 * np.pi * RPM) / 60.0
    Q_real = Q0 * eta_v
    H_inf = Hn * eta_h
    nq = RPM * np.sqrt(Q0) / ((g * Hn) ** 0.75)

    R_m = np.sqrt((R_tip ** 2 + R_hub ** 2) / 2.0)
    z_vec = np.linspace(0.0, -L_z, N_cuerda)
    r_vec = np.linspace(R_hub, R_tip, N_radios)
    Area_paso = np.pi * (R_tip ** 2 - R_hub ** 2)
    V_z = Q_real / Area_paso

    X_mid = np.zeros((N_cuerda, N_radios))
    Y_mid = np.zeros((N_cuerda, N_radios))
    Z_mid = np.zeros((N_cuerda, N_radios))
    b1_vec = np.zeros(N_radios)
    b2_vec = np.zeros(N_radios)

    for i in range(N_radios):
        r = r_vec[i]
        U = omega * r
        V_theta1 = (g * H_inf) / (omega * r)

        if V_theta1 >= U:
            return {'error': (f"Physical boundary instability at r={r:.2f} m: "
                               f"V_theta1 ({V_theta1:.2f} m/s) >= U ({U:.2f} m/s). "
                               f"Increase RPM or decrease head.")}

        beta1 = np.arctan2(V_z, (U - V_theta1))
        beta2 = np.arctan2(V_z, U)
        b1_vec[i] = np.degrees(beta1)
        b2_vec[i] = np.degrees(beta2)

        t_dim = np.abs(z_vec) / L_z
        polinomio = evaluate_interpolation(interp_type, t_dim)
        beta_z = beta1 + (beta2 - beta1) * polinomio

        theta_rel = 0.0
        for j in range(N_cuerda):
            if j > 0:
                dz = abs(z_vec[j] - z_vec[j - 1])
                theta_rel -= rot_sign * (1.0 / np.tan(beta_z[j]) / r) * dz
            X_mid[j, i] = r * np.cos(theta_rel)
            Y_mid[j, i] = r * np.sin(theta_rel)
            Z_mid[j, i] = z_vec[j]

    idx_mid = int(np.argmin(np.abs(r_vec - R_m)))
    dx = np.diff(X_mid[:, idx_mid])
    dy = np.diff(Y_mid[:, idx_mid])
    dz = np.diff(Z_mid[:, idx_mid])
    L_chord_reference = float(np.sum(np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)))

    Z_optimo = int(max(3, min(round((2 * np.pi * R_m) / (L_chord_reference / sigma_target)), 12)))

    R_color = np.tile(r_vec, (N_cuerda, 1))

    return {
        'error': None,
        'X_mid': X_mid, 'Y_mid': Y_mid, 'Z_mid': Z_mid,
        'b1_vec': b1_vec, 'b2_vec': b2_vec,
        'R_color': R_color, 'radius_vec': r_vec,
        'Z_optimo': Z_optimo, 'L_chord_reference': L_chord_reference,
        'rm_label': 'Mean Hydraulic Radius (R_m)', 'rm_val': f'{R_m:.3f} m',
        'eta_t': eta_t, 'omega': omega, 'Q_real': Q_real, 'H_inf': H_inf, 'nq': nq,
        'hub_radius': R_hub, 'tip_radius': R_tip, 'L_z': L_z,
    }


def compute_deriaz_backbone(RPM, Q0, Hn, g, eta_h, eta_v, eta_o, sigma_target,
                             interp_type, rot_sign, N_radios, N_cuerda,
                             Re_int, Re_ext, gamma1_deg, gamma2_deg):
    """SIMPLIFIED diagonal (Deriaz) hydrodynamic model: the same
    Euler-equation-based approach as compute_kaplan_backbone() above,
    adapted to the diagonal (spherical) flow geometry that gives the Deriaz
    turbine its double-regulation capability.

    KEY GEOMETRIC DIFFERENCE FROM KAPLAN: instead of streamlines living on
    cylinders of constant radius, Deriaz streamlines live on SPHERES of
    constant radius Re (the mean spherical hub/casing radii). Each
    streamline is parametrized by the cone angle gamma, sweeping from
    gamma1 (leading edge) to gamma2 (trailing edge); the projected
    cylindrical radius at any point along a streamline is
    rc(gamma) = Re * cos(gamma), which is what actually enters the
    peripheral-speed and blade-angle formulas below (this is the direct
    analogue of "r" in the Kaplan case, but it now varies along the
    streamline instead of staying constant).

    PHYSICS, at each spherical radius Re in [Re_int, Re_ext]:
      1. Inlet/outlet projected radii and peripheral speeds:
             rc1 = Re*cos(gamma1),  U1 = omega*rc1
             rc2 = Re*cos(gamma2),  U2 = omega*rc2
      2. Meridional velocity from continuity through the spherical annulus
         cross-section, evaluated at inlet and outlet respectively.
      3. Required inlet tangential velocity from the Euler equation with
         zero exit swirl: V_theta1 = g*H_inf / U1 (same physical
         infeasibility check as Kaplan if this would reach or exceed U1).
      4. Blade angles from the velocity triangle at each end, and the
         chordwise beta(gamma) distribution obtained the same way as
         Kaplan, via evaluate_interpolation().
      5. The backbone is integrated station by station along gamma
         (instead of along z), using
         d(theta)/d(gamma) = -cot(beta(gamma)) / cos(gamma), the spherical
         analogue of the Kaplan cylindrical integration.

    Returns the same kind of result dict as compute_kaplan_backbone(),
    with Deriaz-specific fields (Re_int, Re_ext, gamma1, gamma2) added and
    'R_color'/'radius_vec' populated with the projected cylindrical radius
    field, used for the radius-based colour mapping in the 3D viewer."""
    if Re_int >= Re_ext or gamma1_deg >= gamma2_deg:
        return {'error': 'Geometric bounds error: check spherical radii or slope angle bounds.'}

    eta_t = eta_h * eta_v * eta_o
    omega = (2 * np.pi * RPM) / 60.0
    Q_real = Q0 * eta_v
    H_inf = Hn * eta_h
    nq = RPM * np.sqrt(Q0) / ((g * Hn) ** 0.75)

    gamma1 = np.radians(gamma1_deg)
    gamma2 = np.radians(gamma2_deg)
    Re_medio = (Re_ext - Re_int) / np.log(Re_ext / Re_int)

    Re_vec = np.linspace(Re_int, Re_ext, N_radios)
    gamma_vec = np.linspace(gamma1, gamma2, N_cuerda)

    X_mid = np.zeros((N_cuerda, N_radios))
    Y_mid = np.zeros((N_cuerda, N_radios))
    Z_mid = np.zeros((N_cuerda, N_radios))
    RC_surf = np.zeros((N_cuerda, N_radios))
    b1_vec = np.zeros(N_radios)
    b2_vec = np.zeros(N_radios)

    for i in range(N_radios):
        Re = Re_vec[i]
        rc1 = Re * np.cos(gamma1)
        U1 = omega * rc1
        Vm1 = Q_real / (2 * np.pi * Re * np.cos(gamma1) * (Re_ext - Re_int))
        Vtheta1 = (g * H_inf) / U1

        if Vtheta1 >= U1:
            return {'error': (f"Physical boundary instability at Re={Re:.2f} m: "
                               f"V_theta1 ({Vtheta1:.2f} m/s) >= U1 ({U1:.2f} m/s). "
                               f"Increase RPM or decrease head.")}

        beta1 = np.arctan2(Vm1, (U1 - Vtheta1))
        rc2 = Re * np.cos(gamma2)
        U2 = omega * rc2
        Vm2 = Q_real / (2 * np.pi * Re * np.cos(gamma2) * (Re_ext - Re_int))
        beta2 = np.arctan2(Vm2, U2)
        b1_vec[i] = np.degrees(beta1)
        b2_vec[i] = np.degrees(beta2)

        t = (gamma_vec - gamma1) / (gamma2 - gamma1)
        polinomio = evaluate_interpolation(interp_type, t)
        beta_gamma = beta1 + (beta2 - beta1) * polinomio

        theta_rel = 0.0
        for j in range(N_cuerda):
            gamma_actual = gamma_vec[j]
            if j > 0:
                d_gamma = gamma_vec[j] - gamma_vec[j - 1]
                theta_rel -= rot_sign * (1.0 / np.tan(beta_gamma[j]) / np.cos(gamma_actual)) * d_gamma
            rc_local = Re * np.cos(gamma_actual)
            X_mid[j, i] = rc_local * np.cos(theta_rel)
            Y_mid[j, i] = rc_local * np.sin(theta_rel)
            Z_mid[j, i] = -Re * np.sin(gamma_actual)
            RC_surf[j, i] = rc_local

    idx_mid = int(np.argmin(np.abs(Re_vec - Re_medio)))
    dx = np.diff(X_mid[:, idx_mid])
    dy = np.diff(Y_mid[:, idx_mid])
    dz = np.diff(Z_mid[:, idx_mid])
    L_chord_reference = float(np.sum(np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)))

    rc_promedio_medio = float(np.mean(RC_surf[:, idx_mid]))
    Z_optimo = int(max(4, min(round((2 * np.pi * rc_promedio_medio) / (L_chord_reference / sigma_target)), 12)))

    return {
        'error': None,
        'X_mid': X_mid, 'Y_mid': Y_mid, 'Z_mid': Z_mid,
        'b1_vec': b1_vec, 'b2_vec': b2_vec,
        'R_color': RC_surf, 'radius_vec': Re_vec,
        'Z_optimo': Z_optimo,
        'L_chord_reference': L_chord_reference,
        'rm_label': 'Mean Spherical Radius (Re_mean)', 'rm_val': f'{Re_medio:.3f} m',
        'eta_t': eta_t, 'omega': omega, 'Q_real': Q_real, 'H_inf': H_inf, 'nq': nq,
        'Re_int': Re_int, 'Re_ext': Re_ext, 'gamma1': gamma1, 'gamma2': gamma2,
    }

def generate_hydro_profile(profile_type, n_points, t_c_base, i_radio, N_radios,
                            hub_to_tip_ratio, m_hub=2.0, p_hub=4.0, m_tip=0.0, p_tip=4.0):
    """Generates the 2D hydrofoil camber (yc) and thickness (yt) distributions,
    as functions of normalized chordwise position x_norm in [0, 1], for a
    single radial station i_radio (out of N_radios total stations).

    The relative thickness t/c is tapered linearly from t_c_base at the hub
    (r_star=0) down to t_c_base*hub_to_tip_ratio at the tip (r_star=1), and
    an analogous (1-r_star)^1.5 factor tapers camber-based profiles toward
    the tip. Five profile families are supported, matching the options
    exposed in the GUI:
      - NACA 00XX (Standard Symmetric): classic 4-digit NACA thickness
        formula, zero camber (yc=0).
      - Reversible Hydrofoil (Pump-Turbine): symmetric elliptical thickness
        (yt = t_c*sqrt(x*(1-x))), zero camber - suited to bidirectional
        (pump/turbine) flow.
      - Anti-Cavitation (Flat Pressure): forward-loaded thickness profile
        (yt = 2.6896*t_c*x*(1-x)^1.5) with camber derived from half the
        thickness distribution, aimed at flattening the pressure-side
        pressure distribution to mitigate cavitation.
      - Low-Torque S-Camber: an S-shaped camber line (cubic in x, changing
        sign) combined with standard NACA thickness, intended to reduce
        the net pitching torque on the blade actuation mechanism.
      - Customized (4-Digit NACA Series): user-specified NACA camber
        parameters (max camber m, its chordwise position p) at hub and
        tip, linearly blended across the span via the same tip-taper
        factor, with standard NACA thickness.

    Returns (yc, yt): two n_points-long arrays giving the camber line and
    half-thickness distribution at this radial station, both already
    scaled to physical units (not normalized by chord)."""
    x_norm = np.linspace(0.0, 1.0, n_points)
    r_star = (i_radio / (N_radios - 1)) if N_radios > 1 else 0.0

    t_c = t_c_base * (1 - (1 - hub_to_tip_ratio) * r_star)
    m_factor = (1 - r_star) ** 1.5

    naca_thickness = (t_c / 0.2) * (0.2969 * np.sqrt(x_norm) - 0.1260 * x_norm
                                     - 0.3516 * x_norm ** 2 + 0.2843 * x_norm ** 3
                                     - 0.1015 * x_norm ** 4)

    if 'Reversible Hydrofoil' in profile_type:
        yc = np.zeros_like(x_norm)
        yt = t_c * np.sqrt(np.clip(x_norm * (1 - x_norm), 0, None))

    elif 'Anti-Cavitation' in profile_type:
        yt = (2.6896 * t_c) * x_norm * (1 - x_norm) ** 1.5
        yc = (yt / 2) * m_factor

    elif 'Low-Torque S-Camber' in profile_type:
        yc = 0.8 * x_norm * (1 - x_norm) * (0.5 - x_norm) * m_factor
        yt = naca_thickness

    elif 'Customized' in profile_type:
        m_hub_val = m_hub / 100.0
        p_hub_val = p_hub / 10.0
        m_tip_val = m_tip / 100.0
        p_tip_val = p_tip / 10.0

        m = m_tip_val + (m_hub_val - m_tip_val) * m_factor
        p = p_tip_val + (p_hub_val - p_tip_val) * m_factor
        p = max(p, 0.05)

        yc = np.zeros_like(x_norm)
        if abs(m) > 0:
            front = x_norm < p
            yc[front] = (m / p ** 2) * (2 * p * x_norm[front] - x_norm[front] ** 2)
            back = ~front
            yc[back] = (m / (1 - p) ** 2) * ((1 - 2 * p) + 2 * p * x_norm[back] - x_norm[back] ** 2)
        yt = naca_thickness

    else:  # 'NACA 00XX (Standard Symmetric)' / default fallback
        yc = np.zeros_like(x_norm)
        yt = naca_thickness

    return yc, yt


def compute_surface_normals(X, Y, Z):
    """Classic full surface-normal computation: cross(Tu, Tv) of the two
    parametric tangent directions, normalized to unit length.

    NOT USED ELSEWHERE IN THIS FILE - kept for reference only. This was the
    ORIGINAL method used to offset the blade thickness away from the mean
    surface, and it has a known geometric flaw: it mixes the chordwise
    curvature of the profile with the blade's spanwise twist rate between
    radial stations, which makes a solid blade generated this way visibly
    "lean sideways" as its thickness increases (the twist component grows
    proportionally with thickness). It was superseded by
    compute_thickness_direction() below, which instead builds a direction
    confined to the local chord/circumferential plane and does not exhibit
    this artifact. See that function's own docstring for the full
    explanation and the fix."""
    dXdu, dXdv = np.gradient(X)
    dYdu, dYdv = np.gradient(Y)
    dZdu, dZdv = np.gradient(Z)

    Tu = np.stack([dXdu, dYdu, dZdu], axis=-1)
    Tv = np.stack([dXdv, dYdv, dZdv], axis=-1)

    N = np.cross(Tu, Tv)
    norm = np.linalg.norm(N, axis=-1, keepdims=True)
    norm[norm == 0] = 1.0
    N = N / norm

    return N[..., 0], N[..., 1], N[..., 2]


def compute_thickness_direction(X, Y, Z):
    """
    Thickness direction used to offset the extrados/intrados away from the
    mean surface.

    Instead of using the FULL normal of the mean surface (cross(Tu, Tv),
    which mixes chordwise curvature with the blade's twist rate between
    radial stations), a direction is built that lives within the local
    section plane (chord + circumferential), matching how 2D airfoil
    sections are classically defined when stacked along radial stations.
    This prevents the thickness from "dragging" a twist component that
    grows proportionally with thickness and makes the blade appear to lean
    sideways (no longer parallel to the inner surface) as it gets thicker -
    the unwanted effect observed with thick profiles.

    Only affects how the SOLID (extrados/intrados) is generated from the
    mean surface; the mean surface itself (backbone) is not modified.
    """
    dXdu, dXdv = np.gradient(X)
    dYdu, dYdv = np.gradient(Y)
    dZdu, dZdv = np.gradient(Z)

    Tu = np.stack([dXdu, dYdu, dZdu], axis=-1)   # tangente a lo largo de la cuerda
    Tv = np.stack([dXdv, dYdv, dZdv], axis=-1)   # tangente a lo largo del radio/envergadura

    def _normalize(V):
        n = np.linalg.norm(V, axis=-1, keepdims=True)
        n = np.where(n == 0, 1.0, n)
        return V / n

    Tu_hat = _normalize(Tu)

    # Local circumferential direction (tangent to the circle of radius
    # sqrt(X^2+Y^2) centered on the Z axis) - valid for both axial
    # (Kaplan) and conical (Deriaz) geometry, since both are surfaces
    # of revolution about the Z axis.
    r_cyl = np.sqrt(X ** 2 + Y ** 2)
    r_safe = np.where(r_cyl < 1e-9, 1.0, r_cyl)
    e_theta = np.stack([-Y / r_safe, X / r_safe, np.zeros_like(X)], axis=-1)

    dot = np.sum(e_theta * Tu_hat, axis=-1, keepdims=True)
    N = e_theta - dot * Tu_hat
    N_norm = np.linalg.norm(N, axis=-1, keepdims=True)

    # Fall back to the classic surface normal in the degenerate case
    # where e_theta ends up nearly parallel to the chordwise tangent.
    fallback = _normalize(np.cross(Tu, Tv))
    degenerate = N_norm < 1e-6
    N_safe = N / np.where(N_norm == 0, 1.0, N_norm)
    N_final = np.where(degenerate, fallback, N_safe)

    return N_final[..., 0], N_final[..., 1], N_final[..., 2]


# =============================================================================
# INDUSTRIAL POST-PROCESSING ADJUSTMENTS (optional, off by default)
# =============================================================================
# Several geometric operations, all applied to the mean-surface backbone
# BEFORE the hydrofoil thickness is generated, so build_solid_blade() always
# regenerates the profile fresh on the already-adjusted mean line and tapers
# correctly at any new edge created by these operations. They approximate
# real CAD/manufacturing finishing steps applied to a computed "ideal"
# hydrodynamic runner before it can be built.
#
# NONE of these operations is simply "the correct one" with the others being
# lesser approximations - every one of them trades some property of the pure
# hydrodynamic design for another, and which trade-off is preferable is a
# genuine engineering decision, not a bug to be fixed. Each function's
# docstring below states plainly what is gained and what is given up.
#
# Only the leading-edge sweep (Kaplan only - see its docstring for why) is
# free of any cost at all.
# =============================================================================

def trim_backbone_by_arclength_variable(X_mid, Y_mid, Z_mid, chord_fraction_per_streamline):
    """Resamples every streamline of the mean-surface backbone onto a new
    parametrization, keeping only a PER-STREAMLINE fraction of each
    streamline's own arc length STARTING FROM THE LEADING EDGE (i.e. this
    trims from the trailing-edge side only), re-sampled with the same
    number of points evenly spaced in arc length.

    This is the low-level building block behind trim_backbone_by_arclength()
    (uniform trailing-edge wake trim): a thin wrapper that calls this with
    the same fraction for every streamline. For an operation that can trim
    from BOTH ends of a streamline (needed to round a corner symmetrically
    at both the leading- and trailing-edge side), see
    trim_backbone_by_arclength_window() below instead."""
    N_cuerda, N_radios = X_mid.shape
    X_new = np.zeros_like(X_mid)
    Y_new = np.zeros_like(Y_mid)
    Z_new = np.zeros_like(Z_mid)

    for i in range(N_radios):
        dx = np.diff(X_mid[:, i])
        dy = np.diff(Y_mid[:, i])
        dz = np.diff(Z_mid[:, i])
        ds = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        s_cum = np.concatenate([[0.0], np.cumsum(ds)])
        total_length = s_cum[-1]
        target_length = chord_fraction_per_streamline[i] * total_length

        s_new = np.linspace(0.0, target_length, N_cuerda)
        X_new[:, i] = np.interp(s_new, s_cum, X_mid[:, i])
        Y_new[:, i] = np.interp(s_new, s_cum, Y_mid[:, i])
        Z_new[:, i] = np.interp(s_new, s_cum, Z_mid[:, i])

    return X_new, Y_new, Z_new


def trim_backbone_by_arclength_window(X_mid, Y_mid, Z_mid, start_fraction_per_streamline,
                                       end_fraction_per_streamline, min_window_frac=0.02):
    """Generalizes trim_backbone_by_arclength_variable() to trim from BOTH
    ends of every streamline at once: keeps only the arc-length WINDOW
    [start_fraction, end_fraction] of each streamline's own total length
    (both fractions given per streamline), re-sampled with the same number
    of points evenly spaced across the kept window.

    Used by apply_corner_rounding() below to round a corner symmetrically
    at both its leading-edge side and its trailing-edge side. A safety
    floor (min_window_frac) guarantees the kept window can never invert or
    collapse to (near) zero length, however aggressively start/end fractions
    are requested - if they would, the window is recentred on their
    midpoint and widened back out to the floor width."""
    N_cuerda, N_radios = X_mid.shape
    X_new = np.zeros_like(X_mid)
    Y_new = np.zeros_like(Y_mid)
    Z_new = np.zeros_like(Z_mid)

    for i in range(N_radios):
        dx = np.diff(X_mid[:, i])
        dy = np.diff(Y_mid[:, i])
        dz = np.diff(Z_mid[:, i])
        ds = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        s_cum = np.concatenate([[0.0], np.cumsum(ds)])
        total_length = s_cum[-1]

        f0 = start_fraction_per_streamline[i]
        f1 = end_fraction_per_streamline[i]
        if f1 - f0 < min_window_frac:
            mid = 0.5 * (f0 + f1)
            f0, f1 = mid - min_window_frac / 2.0, mid + min_window_frac / 2.0
            f0, f1 = max(f0, 0.0), min(f1, 1.0)

        s0 = f0 * total_length
        s1 = f1 * total_length
        s_new = np.linspace(s0, s1, N_cuerda)
        X_new[:, i] = np.interp(s_new, s_cum, X_mid[:, i])
        Y_new[:, i] = np.interp(s_new, s_cum, Y_mid[:, i])
        Z_new[:, i] = np.interp(s_new, s_cum, Z_mid[:, i])

    return X_new, Y_new, Z_new


def trim_backbone_by_arclength(X_mid, Y_mid, Z_mid, chord_fraction_kept):
    """Trims the trailing-edge wake off the mean-surface backbone: removes
    the last (1 - chord_fraction_kept) of every streamline's OWN arc
    length, uniformly across the whole span.

    No commercial Deriaz or Kaplan runner is left in service with a bare,
    untrimmed wake: the trailing edge is always cut back, and the airfoil
    thickness is re-tapered on the shorter chord that results. Because this
    function trims by a common ARC-LENGTH FRACTION (see
    trim_backbone_by_arclength_variable() above for why that matters) and
    the hydrofoil profile is always regenerated AFTER this trim runs, on
    the new shorter chord, the thickness naturally tapers back to (near)
    zero at the new trailing edge with no special-casing required."""
    N_radios = X_mid.shape[1]
    fractions = np.full(N_radios, chord_fraction_kept)
    return trim_backbone_by_arclength_variable(X_mid, Y_mid, Z_mid, fractions)


def apply_leading_edge_sweep(X_mid, Y_mid, Z_mid, radius_vec, sweep_angle_deg):
    """Rigid axial (Z) translation of each streamline, growing linearly with
    radial distance from the hub, at a rate set by a true sweep angle in
    the standard aerospace/turbomachinery sense:
        tan(sweep_angle) = axial offset / spanwise distance from the hub.

    KAPLAN ONLY - DO NOT USE ON A DERIAZ BACKBONE. This was verified
    directly: for a Kaplan blade, every streamline lives on a cylinder of
    constant radius r, and a cylinder is invariant under translation along
    its own axis (Z) - a point translated in Z stays on the exact same
    cylinder, U and the whole velocity triangle at that point are governed
    only by r and are completely unaffected. For a DERIAZ blade, streamlines
    instead live on a SPHERE of constant radius Re about the hub/casing
    common centre; a per-streamline constant Z-offset does NOT preserve
    that sphere (only X and Y stay fixed, so the point's true distance to
    the centre, sqrt(X^2+Y^2+Z^2), changes and drifts away from Re for
    every streamline except the one left with zero offset). This was
    confirmed numerically: applying this function to a Deriaz backbone
    pulls every non-hub streamline measurably off its design sphere,
    which breaks the physical mechanism the Deriaz double-regulation
    depends on (blades pivoting on, and clearing, a truly spherical
    hub/casing surface). The GUI only exposes this control for Kaplan;
    there is currently no equivalent sweep operation offered for Deriaz.

    HYDRAULIC COST (Kaplan): none. A per-streamline CONSTANT offset changes
    neither the local chord (arc length is invariant under a rigid
    translation), nor the local beta angle (a derivative-based, purely
    local property), nor the total swept angle theta of the streamline -
    it only repositions the whole streamline in space."""
    Z_new = Z_mid.copy()
    N_radios = X_mid.shape[1]
    for i in range(N_radios):
        offset = -(radius_vec[i] - radius_vec[0]) * np.tan(np.radians(sweep_angle_deg))
        Z_new[:, i] = Z_mid[:, i] + offset
    return X_mid.copy(), Y_mid.copy(), Z_new


def apply_corner_rounding(X_mid, Y_mid, Z_mid, radius_vec, chord_ref,
                           round_le_tip, le_tip_extent_frac,
                           round_te_tip, te_tip_extent_frac,
                           round_le_root, le_root_extent_frac,
                           round_te_root, te_root_extent_frac,
                           min_chord_floor=0.15):
    """Rounds any combination of the blade planform's four corners
    independently - leading-edge/tip, trailing-edge/tip, leading-edge/
    root, trailing-edge/root - each with its own on/off switch and its
    own fillet extent, using the exact same geometric construction a CAD
    fillet uses to round a rectangular corner, applied via
    trim_backbone_by_arclength_window() above.

    GEOMETRY: think of a single streamline's own local coordinate frame as
    (dr, ds), where dr is the distance from the span boundary (0 exactly at
    the tip or root, growing toward midspan) and ds is the distance from
    one of that streamline's own edges (leading OR trailing). A fillet of
    radius R replaces the sharp (dr=0, ds=0) corner with the arc of a
    circle of radius R centered at (R, R); solving that circle for ds
    gives, for dr in [0, R]:
        ds_boundary(dr) = R - sqrt(R^2 - (R - dr)^2)
    which is exactly how far back from that edge the streamline must be
    trimmed. At dr = 0 (right at the tip/root) the trim equals the full
    radius R; at dr = R (one fillet-radius away) the trim is zero.

    Each of the four corners has its own independent (on/off, extent)
    pair and is blended in separately - so, for example, only the
    trailing-edge/tip corner can be rounded while the other three are
    left sharp, matching a real anti-cavitation lip that may only be
    applied on one side, or any other combination.

    SCALING NOTE: the fillet radius R is defined as a fraction of the
    blade's REFERENCE CHORD length, not of the spanwise extent - chord and
    span are frequently very different lengths, so scaling R by span made
    an earlier version of this function's rounding effectively invisible
    for any reasonable "percent" setting.

    A hard floor (min_chord_floor) caps how much of the local chord any
    SINGLE corner's trim can consume, and trim_backbone_by_arclength_
    window()'s own safety net (its min_window_frac) additionally
    guarantees the final kept window can never collapse to zero or invert
    even if several corners' trims overlap and compound on the same
    streamline (e.g. a short-span blade with both tip corners and both
    root corners all rounded aggressively at once).

    ENGINEERING BASIS (thesis-relevant honesty note): the TIP rounding
    approximates the anti-cavitation lip that many commercial Kaplan
    runners carry on the suction side near the periphery, a documented
    feature used to mitigate tip-clearance/tip-vortex cavitation erosion.
    The ROOT rounding is offered for the same geometric flexibility and is
    mechanically plausible (reducing stress concentration at the blade
    root, analogous to a standard mechanical fillet), but no
    Kaplan/Deriaz-specific published source was found to confirm it as a
    normalized industrial practice - its default extent should be treated
    as illustrative only, unlike the tip option."""
    N_radios = X_mid.shape[1]
    chord_vec = np.array([
        np.sum(np.sqrt(np.diff(X_mid[:, i]) ** 2 + np.diff(Y_mid[:, i]) ** 2 + np.diff(Z_mid[:, i]) ** 2))
        for i in range(N_radios)
    ])
    start_frac = np.zeros(N_radios)  # how far in from the LEADING edge to trim
    end_frac = np.ones(N_radios)     # how far in from the TRAILING edge to trim (as a fraction kept)

    def _corner_trim_frac(i, extent_frac, at_tip):
        R = extent_frac * chord_ref
        dr = (radius_vec[-1] - radius_vec[i]) if at_tip else (radius_vec[i] - radius_vec[0])
        if dr >= R:
            return 0.0
        u = R - dr
        ds_boundary = R - np.sqrt(max(R ** 2 - u ** 2, 0.0))
        local_chord = chord_vec[i]
        # Cap this single corner's own trim so it alone can never consume
        # more than (1 - min_chord_floor) of the local chord; the window
        # function's own floor additionally guards against several
        # corners compounding past that on the same streamline.
        return min(ds_boundary / local_chord, 1.0 - min_chord_floor)

    def _blend_le(round_this, extent_frac, at_tip):
        if not round_this or extent_frac <= 0:
            return
        for i in range(N_radios):
            trim_frac = _corner_trim_frac(i, extent_frac, at_tip)
            start_frac[i] = max(start_frac[i], trim_frac)

    def _blend_te(round_this, extent_frac, at_tip):
        if not round_this or extent_frac <= 0:
            return
        for i in range(N_radios):
            trim_frac = _corner_trim_frac(i, extent_frac, at_tip)
            end_frac[i] = min(end_frac[i], 1.0 - trim_frac)

    _blend_le(round_le_tip, le_tip_extent_frac, at_tip=True)
    _blend_te(round_te_tip, te_tip_extent_frac, at_tip=True)
    _blend_le(round_le_root, le_root_extent_frac, at_tip=False)
    _blend_te(round_te_root, te_root_extent_frac, at_tip=False)

    return trim_backbone_by_arclength_window(X_mid, Y_mid, Z_mid, start_frac, end_frac)


def align_trailing_edge_radially(X_mid, Y_mid, Z_mid, reference_X=None, reference_Y=None, reference_Z=None):
    """Cuts the mean-surface backbone with a single flat plane and keeps
    only the portion of the blade behind it, on either a Kaplan or a
    Deriaz backbone.

    THE CUTTING PLANE, GEOMETRICALLY: it is built from two lines, both
    passing through the rotation axis:
      1. The rotation axis itself (the Z axis).
      2. The line from the origin to the ACTUAL, natural trailing-edge
         point of the HUB's own streamline (for Kaplan: the point at
         r=R_hub, z=-L_z, at whatever circumferential angle theta that
         streamline's own twist integration naturally produced; for
         Deriaz: the equivalent point at Re=Re_int, gamma=gamma2).
    Two lines through a common point define a plane; because both lines
    here pass through the axis, that plane is exactly a constant-theta
    half-plane pair (theta_hub and theta_hub+180 degrees) - i.e. a flat
    plane containing the shaft axis, at the SAME circumferential angle as
    the hub's own trailing edge. "Behind" the cut (the direction the flow
    has already travelled through) is what remains after this function
    runs; "in front of" the cut (any leftover chord past that theta) is
    the trimmed-away wake.

    Because this plane always passes through the axis, cutting a
    streamline with it is equivalent to trimming that streamline's arc
    length until its own theta first reaches theta_hub - which is exactly
    what this function computes for every streamline, using the SAME
    theta_hub for the whole span. The hub itself is therefore always the
    anchor, by construction, not "whichever streamline happens to sweep
    the least" (an earlier version of this function picked the anchor
    that way; anchoring explicitly at the hub is simpler, more
    predictable, and matches how a real double-regulated runner's blade
    root is normally used as the fixed reference for this kind of finish
    cut).

    SAME CONSTRUCTION FOR KAPLAN AND DERIAZ, VERIFIED SPHERE-PRESERVING.
    The plane-cut idea does not depend on whether the underlying
    streamlines live on cylinders (Kaplan) or spheres (Deriaz): cutting at
    a constant theta only ever trims arc length along each streamline's
    OWN already-computed curve, it never displaces a point off of that
    curve, so a Deriaz streamline trimmed this way stays exactly on its
    own design sphere Re (confirmed numerically: maximum deviation on the
    order of 1e-5 m, pure floating-point/interpolation noise, across every
    streamline on a representative case).

    reference_X/Y/Z: OPTIONAL pristine backbone (before ANY other
    post-processing step has touched it) used ONLY to compute theta_hub.
    THIS MATTERS: if this function runs after apply_corner_rounding() has
    already shortened the hub streamline, computing theta_hub from the
    ALREADY-TRIMMED geometry would silently anchor the plane at a point
    BEFORE the hub's true design trailing edge (z=-L_z for Kaplan, gamma2
    for Deriaz) rather than at it - confirmed numerically in an earlier
    version of this function (a ~26 degree difference between the buggy,
    already-trimmed-derived anchor and the correct one on a representative
    case). Always pass the untouched, pristine backbone as
    reference_X/Y/Z when this runs alongside any other trimming operation
    in the post-processing chain; if omitted, theta_hub is computed from
    X_mid/Y_mid/Z_mid themselves (only safe when nothing has trimmed them
    yet). This is also why, in the GUI's post-processing pipeline, corner
    rounding is applied AFTER this function runs, never before: rounding
    the hub corner first and only then computing theta_hub from the
    already-rounded hub would anchor the cutting plane at the wrong angle,
    and rounding the tip/root corners of a blade that this function has
    NOT yet trimmed means the rounded region can later be cut away
    entirely by this function, silently erasing the rounding - applying
    this cut first and rounding the resulting, already-cut edges
    afterward avoids both problems.

    ABSOLUTE THETA, NOT "THETA MINUS THETA[0]" - A SECOND BUG THIS AVOIDS.
    An earlier version of this function measured each streamline's swept
    angle as theta - theta[0], using whatever the FIRST POINT of the
    CURRENT array happens to be as the zero reference. That is only safe
    if index 0 still corresponds to the true design leading edge. If
    apply_corner_rounding() has ALREADY resampled this same streamline
    (trimming its leading-edge side too, to round that corner), the
    current array's index 0 no longer sits at the true design leading
    edge at all - it sits wherever the rounding happened to cut to, which
    is a point where theta is no longer 0. Using "theta - theta[0]" there
    silently measures the sweep from the WRONG zero point, while
    theta_hub above was computed from the reference geometry's TRUE zero
    (confirmed numerically: a rounded tip streamline's own theta[0] was
    -19.7 degrees instead of 0). Mixing those two different zero points
    is what produced visibly broken, "bow-tie" shaped blades when corner
    rounding and this alignment were both enabled together in an earlier
    pipeline ordering. The fix is to use theta directly (no offset): by
    construction, every streamline in this application's coordinate
    convention starts its integration at theta = 0 exactly at its own true
    design leading edge (both backbone functions initialize theta_rel =
    0.0 before integrating), so raw theta already measures "how far this
    point has swept from the true design leading edge" correctly,
    regardless of whether an earlier step has since resampled the array
    and moved what index 0 means.

    NEITHER OPTION HERE IS SIMPLY "THE CORRECT ONE" - THIS IS A GENUINE
    TWO-SIDED TRADE-OFF, NOT A FIX FOR A "WRONG" GEOMETRY:

    Leaving this control OFF keeps the natural, hydrodynamically pure
    geometry: every streamline integrated exactly to its own Euler-
    equation design point, with zero artificial energy loss (all Kaplan
    streamlines already share the same axial extent by construction; all
    Deriaz streamlines already share the same spherical angle gamma2).
    Its cost is geometric/structural instead: theta at the discharge edge
    varies substantially across the span (often tens of degrees between
    hub and tip), which can concentrate stress at that edge, promote
    vibration, and leave an irregular wake shape downstream - real
    problems, not merely cosmetic ones.

    Applying THIS function removes that geometric problem by cutting with
    the plane described above, at the cost of trimming every streamline
    other than the hub short of its own design point: that streamline
    then delivers less specific energy than intended and discharges with
    a residual tangential velocity component, itself a possible seed for
    downstream vortices and draft-tube pressure pulsations."""
    N_cuerda, N_radios = X_mid.shape
    ref_X = X_mid if reference_X is None else reference_X
    ref_Y = Y_mid if reference_Y is None else reference_Y

    # theta_hub: the circumferential angle of the HUB's own natural
    # trailing edge (last point of streamline index 0 in the pristine
    # reference geometry) - this single angle defines the cutting plane,
    # together with the rotation axis, for every streamline on the blade.
    theta_hub_curve = np.unwrap(np.arctan2(ref_Y[:, 0], ref_X[:, 0]))
    target = theta_hub_curve[-1]  # ref theta[0] is 0 by construction

    X_new = np.zeros_like(X_mid)
    Y_new = np.zeros_like(Y_mid)
    Z_new = np.zeros_like(Z_mid)
    for i in range(N_radios):
        theta = np.unwrap(np.arctan2(Y_mid[:, i], X_mid[:, i]))  # NOT offset by theta[0]
        dx = np.diff(X_mid[:, i])
        dy = np.diff(Y_mid[:, i])
        dz = np.diff(Z_mid[:, i])
        ds = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        s_cum = np.concatenate([[0.0], np.cumsum(ds)])

        # np.interp() requires an increasing x-array; theta may be
        # monotonically decreasing depending on the rotation-direction sign
        # convention, so flip both arrays together in that case.
        if theta[-1] < theta[0]:
            theta_for_interp = theta[::-1]
            s_for_interp = s_cum[::-1]
        else:
            theta_for_interp = theta
            s_for_interp = s_cum

        target_len = np.interp(target, theta_for_interp, s_for_interp)
        target_len = min(max(target_len, 0.0), s_cum[-1])  # safety clamp

        s_new = np.linspace(0.0, target_len, N_cuerda)
        X_new[:, i] = np.interp(s_new, s_cum, X_mid[:, i])
        Y_new[:, i] = np.interp(s_new, s_cum, Y_mid[:, i])
        Z_new[:, i] = np.interp(s_new, s_cum, Z_mid[:, i])

    return X_new, Y_new, Z_new, target



def build_solid_blade(X_mid, Y_mid, Z_mid, profile_type, t_rel_input,
                       N_cuerda, N_radios, hub_to_tip_ratio,
                       m_hub=2.0, p_hub=4.0, m_tip=0.0, p_tip=4.0):
    """Thickens the 2D mean-surface backbone (X_mid, Y_mid, Z_mid) into a
    3D solid blade by offsetting an extrados (pressure side) and an
    intrados (suction side) surface away from it, along the local
    thickness direction returned by compute_thickness_direction().

    For each radial station i:
      1. The local chord length c_i is measured as the cumulative arc
         length of that streamline (s_backbone), NOT assumed constant -
         this is what lets the thickness taper correctly on a backbone
         that has already been trimmed/reshaped by any of the industrial
         post-processing operations above, since those change c_i directly.
      2. The absolute reference thickness (t_abs_ref), fixed once from the
         user-specified relative thickness at the HUB chord, is converted
         back into a LOCAL relative thickness t_rel_local = t_abs_ref / c_i
         - so the absolute (metric) thickness stays roughly constant across
         the span even though the relative thickness the profile generator
         actually uses varies with local chord.
      3. generate_hydro_profile() returns the normalized camber (yc) and
         half-thickness (yt) distributions for this station and profile
         family; camber_offset and thick_offset scale them back to
         physical units using c_i.
      4. The extrados and intrados points are placed at
         mean_surface +/- (camber_offset + thick_offset) along the local
         thickness direction (camber shifts both sides together, thickness
         pushes them apart).
      5. The hub column (i=0) of both surfaces is then re-projected back
         onto the mean surface's own local radius, exactly - see the
         "HUB FLUSH-ATTACHMENT CORRECTION" comment below for why this is
         needed even though the thickness direction has no radial
         component at the hub.

    Returns the six (N_cuerda x N_radios) coordinate arrays for the
    extrados and intrados surfaces."""
    s_backbone = np.zeros((N_cuerda, N_radios))
    dx_m = np.diff(X_mid, axis=0)
    dy_m = np.diff(Y_mid, axis=0)
    dz_m = np.diff(Z_mid, axis=0)
    ds_m = np.sqrt(dx_m ** 2 + dy_m ** 2 + dz_m ** 2)
    s_backbone[1:, :] = np.cumsum(ds_m, axis=0)
    chord_lengths = s_backbone[-1, :]

    c_hub = chord_lengths[0]
    t_abs_ref = t_rel_input * c_hub

    Nx, Ny, Nz = compute_thickness_direction(X_mid, Y_mid, Z_mid)

    X_ext = np.zeros_like(X_mid); Y_ext = np.zeros_like(X_mid); Z_ext = np.zeros_like(X_mid)
    X_int = np.zeros_like(X_mid); Y_int = np.zeros_like(X_mid); Z_int = np.zeros_like(X_mid)

    for i in range(N_radios):
        c_i = chord_lengths[i]
        t_rel_local = t_abs_ref / c_i

        yc_base, yt_base = generate_hydro_profile(profile_type, N_cuerda, t_rel_local, i, N_radios,
                                                    hub_to_tip_ratio, m_hub, p_hub, m_tip, p_tip)

        camber_offset = yc_base * c_i
        thick_offset = yt_base * c_i

        offset_ext = camber_offset + thick_offset
        offset_int = camber_offset - thick_offset

        nx_i, ny_i, nz_i = Nx[:, i], Ny[:, i], Nz[:, i]

        X_ext[:, i] = X_mid[:, i] + offset_ext * nx_i
        Y_ext[:, i] = Y_mid[:, i] + offset_ext * ny_i
        Z_ext[:, i] = Z_mid[:, i] + offset_ext * nz_i

        X_int[:, i] = X_mid[:, i] + offset_int * nx_i
        Y_int[:, i] = Y_mid[:, i] + offset_int * ny_i
        Z_int[:, i] = Z_mid[:, i] + offset_int * nz_i

    # --- HUB FLUSH-ATTACHMENT CORRECTION (re-projection at i=0) ---
    # Even though the thickness direction is purely tangential at the hub
    # (zero radial component, by construction: every hub streamline point
    # sits at a fixed cylindrical/spherical radius, so both its chordwise
    # tangent and the circumferential direction e_theta are automatically
    # perpendicular to the radial direction there), offsetting a point by
    # a FINITE amount in a purely tangential direction still increases its
    # distance from the axis very slightly - the same chord-vs-arc effect
    # that makes a straight chord of a circle always fall strictly inside
    # the circle itself. For a point at radius r offset by a tangential
    # distance d, the new radius is sqrt(r^2+d^2) = r + d^2/(2r) + ...,
    # i.e. a SECOND-ORDER (quadratic-in-thickness) growing gap - confirmed
    # numerically (gap approximately quadruples when thickness doubles).
    # This is invisible at realistic thickness values but becomes
    # measurable at unusually large relative thickness.
    #
    # Since a real blade root must sit exactly flush against the hub
    # surface with no air gap, however small, this correction re-projects
    # ONLY the hub column (i=0) of the extrados/intrados surfaces back
    # onto the mean surface's own true hub geometry at each chordwise
    # station, preserving only the circumferential angle (theta) the
    # tangential offset produced.
    #
    # THIS MUST CORRECT THE FULL 3D POSITION, NOT JUST THE CYLINDRICAL
    # (X,Y) RADIUS - an earlier version of this correction only rescaled
    # (X,Y) to match the mean surface's cylindrical radius sqrt(X^2+Y^2),
    # leaving Z untouched. That is exactly right for a Kaplan hub, which
    # is a CYLINDER (defined purely by a fixed cylindrical radius, with Z
    # free) - but wrong for a Deriaz hub, which is a SPHERE (defined by a
    # fixed distance sqrt(X^2+Y^2+Z^2) from the ORIGIN, not just from the
    # axis): correcting only the cylindrical radius there left Z
    # untouched, so the corrected point's true 3D distance from the
    # origin still drifted off Re_int - confirmed numerically, the
    # intrados dipped as far as ~5% under Re_int and the extrados as far
    # as ~2.5% over it at extreme thickness, i.e. exactly the "sinks into
    # the sphere" symptom this fixes. Forcing BOTH the cylindrical radius
    # AND Z to match the mean surface's own values at the same chordwise
    # station is correct for either surface of revolution at once,
    # without needing to know which one it is: for Kaplan, Z_mid is
    # already identical across every streamline regardless of radius, so
    # forcing it here is a no-op; for Deriaz, forcing rc and Z together
    # reproduces the mean surface's own (rc, z) pair exactly, whose
    # rc^2+z^2 already equals Re_int^2 by construction of the backbone
    # itself - so the corrected point's full 3D radius lands on Re_int
    # exactly, not merely its cylindrical projection.
    rc_hub_target = np.sqrt(X_mid[:, 0] ** 2 + Y_mid[:, 0] ** 2)
    z_hub_target = Z_mid[:, 0]

    for X_surf, Y_surf, Z_surf in ((X_ext, Y_ext, Z_ext), (X_int, Y_int, Z_int)):
        theta_offset = np.arctan2(Y_surf[:, 0], X_surf[:, 0])
        X_surf[:, 0] = rc_hub_target * np.cos(theta_offset)
        Y_surf[:, 0] = rc_hub_target * np.sin(theta_offset)
        Z_surf[:, 0] = z_hub_target

    return X_ext, Y_ext, Z_ext, X_int, Y_int, Z_int


def surf2patch_triangles(X, Y, Z):
    """Converts a structured (R x C) point grid into a triangulated mesh:
    every quad cell of the grid is split into two triangles, sharing the
    diagonal from the cell's (row, col) corner to its (row+1, col+1)
    corner. Returns (F, V): F is an (n_triangles x 3) integer array of
    vertex indices into V, and V is an (R*C x 3) array of the flattened
    vertex coordinates (row-major / 'C' order, matching how F's indices
    are computed: index = row * C + col)."""
    R, C = X.shape
    V = np.column_stack([X.ravel(order='C'), Y.ravel(order='C'), Z.ravel(order='C')])
    r_idx, c_idx = np.meshgrid(np.arange(R - 1), np.arange(C - 1), indexing='ij')
    i1 = (r_idx * C + c_idx).ravel()
    i2 = (r_idx * C + c_idx + 1).ravel()
    i3 = ((r_idx + 1) * C + c_idx + 1).ravel()
    i4 = ((r_idx + 1) * C + c_idx).ravel()
    F1 = np.column_stack([i1, i2, i3])
    F2 = np.column_stack([i1, i3, i4])
    F = np.vstack([F1, F2]).astype(np.int64)
    return F, V


def build_export_mesh(X_ext, Y_ext, Z_ext, X_int, Y_int, Z_int, X_mid, Y_mid, Z_mid, is_solid):
    """Assembles the final, exportable triangulated mesh for ONE blade.

    If is_solid is False, simply triangulates the mean surface alone (a
    single, zero-thickness sheet) and returns it - used for the
    "Mean Surface Only" export option.

    If is_solid is True, builds a WATERTIGHT closed solid from six
    triangulated patches:
      - the extrados (pressure side) surface,
      - the intrados (suction side) surface, with its triangle winding
        order flipped ([:, [0,2,1]]) so its normals point outward
        (opposite sense to the extrados, since it faces the other way),
      - four thin "cap" strips closing the remaining open edges: leading
        edge (row 0), trailing edge (row -1, winding flipped), hub end
        (column 0), and tip end (column -1, winding flipped) - each cap is
        itself triangulated by treating the two matching extrados/intrados
        edge curves as a 2-row structured grid via surf2patch_triangles().

    All six patches are concatenated with independently-offset vertex
    index blocks, then np.unique(..., return_inverse=True) welds any
    duplicate vertices shared between adjacent patches (e.g. the leading
    edge cap and the extrados patch legitimately share their first row of
    points) into single, shared mesh vertices - this is what makes the
    result a single watertight manifold instead of six disconnected
    pieces, which matters for any downstream CAD import or 3D printing."""
    if not is_solid:
        return surf2patch_triangles(X_mid, Y_mid, Z_mid)

    F_ext, V_ext = surf2patch_triangles(X_ext, Y_ext, Z_ext)

    F_int, V_int = surf2patch_triangles(X_int, Y_int, Z_int)
    F_int = F_int[:, [0, 2, 1]]

    F_le, V_le = surf2patch_triangles(
        np.vstack([X_ext[0, :], X_int[0, :]]),
        np.vstack([Y_ext[0, :], Y_int[0, :]]),
        np.vstack([Z_ext[0, :], Z_int[0, :]]))

    F_te, V_te = surf2patch_triangles(
        np.vstack([X_ext[-1, :], X_int[-1, :]]),
        np.vstack([Y_ext[-1, :], Y_int[-1, :]]),
        np.vstack([Z_ext[-1, :], Z_int[-1, :]]))
    F_te = F_te[:, [0, 2, 1]]

    F_hub, V_hub = surf2patch_triangles(
        np.column_stack([X_ext[:, 0], X_int[:, 0]]),
        np.column_stack([Y_ext[:, 0], Y_int[:, 0]]),
        np.column_stack([Z_ext[:, 0], Z_int[:, 0]]))

    F_tip, V_tip = surf2patch_triangles(
        np.column_stack([X_ext[:, -1], X_int[:, -1]]),
        np.column_stack([Y_ext[:, -1], Y_int[:, -1]]),
        np.column_stack([Z_ext[:, -1], Z_int[:, -1]]))
    F_tip = F_tip[:, [0, 2, 1]]

    V_all = np.vstack([V_ext, V_int, V_le, V_te, V_hub, V_tip])

    off_ext = 0
    off_int = off_ext + len(V_ext)
    off_le = off_int + len(V_int)
    off_te = off_le + len(V_le)
    off_hub = off_te + len(V_te)
    off_tip = off_hub + len(V_hub)

    F_all = np.vstack([
        F_ext + off_ext, F_int + off_int, F_le + off_le,
        F_te + off_te, F_hub + off_hub, F_tip + off_tip,
    ])

    V_final, inverse = np.unique(V_all, axis=0, return_inverse=True)
    F_final = inverse[F_all]
    return F_final, V_final


def write_stl_ascii(filename, F, V):
    """Writes a mesh (F: triangle vertex-index array, V: vertex coordinate
    array) to an ASCII STL file. Each triangle's outward normal is
    computed directly from its own vertices via cross((p2-p1), (p3-p1)),
    normalized to unit length (degenerate zero-area triangles are written
    with a zero normal rather than raising a division error). ASCII STL
    (as opposed to the more compact binary STL variant) is used for
    maximum compatibility and human-readability, at the cost of a larger
    file size - acceptable given these are single-blade meshes, not
    full assemblies."""
    with open(filename, 'w') as f:
        f.write('solid HydroTurbineBlade\n')
        for tri in F:
            p1, p2, p3 = V[tri[0]], V[tri[1]], V[tri[2]]
            normal = np.cross(p2 - p1, p3 - p1)
            n_len = np.linalg.norm(normal)
            normal = normal / n_len if n_len > 0 else np.zeros(3)
            f.write(f'  facet normal {normal[0]:e} {normal[1]:e} {normal[2]:e}\n')
            f.write('    outer loop\n')
            f.write(f'      vertex {p1[0]:e} {p1[1]:e} {p1[2]:e}\n')
            f.write(f'      vertex {p2[0]:e} {p2[1]:e} {p2[2]:e}\n')
            f.write(f'      vertex {p3[0]:e} {p3[1]:e} {p3[2]:e}\n')
            f.write('    endloop\n')
            f.write('  endfacet\n')
        f.write('endsolid HydroTurbineBlade\n')


def write_iges_fallback(filename, F, V):
    """IGES fallback writer, matching the original MATLAB behaviour: a point
    cloud (type 116 entities), NOT a valid B-Rep solid. For a real
    IGES/STEP solid, run the exported STL through a CAD kernel
    (e.g. FreeCAD/OpenCASCADE)."""
    with open(filename, 'w') as f:
        f.write('S      1\n')
        f.write('G240101.000000;1H;1H;1H;1H;38;15;1;1.0;1;1H;1.0;1H;1H;;;G      1\n')
        p_line = 1
        for tri in F:
            for idx in tri:
                p = V[idx]
                f.write(f'116,{p[0]:f},{p[1]:f},{p[2]:f},0,0,0,1.0,0,0;P{p_line:7d}\n')
                p_line += 1
        f.write(f'T{p_line - 1:6d}\n')


# =============================================================================
# WIDGETS AUXILIARES
# =============================================================================

def make_double_spin(value, minimum, maximum, decimals=2, step=0.1):
    """Creates a QDoubleSpinBox. Deliberately has NO tooltip of its own:
    tooltips live only on the field's label (see add_row), never on the
    input widget itself."""
    sp = QDoubleSpinBox()
    sp.setDecimals(decimals)
    sp.setRange(minimum, maximum)
    sp.setSingleStep(step)
    sp.setValue(value)
    return sp


def make_int_spin(value, minimum, maximum):
    """Creates a QSpinBox. See make_double_spin() docstring re: tooltips."""
    sp = QSpinBox()
    sp.setRange(minimum, maximum)
    sp.setValue(value)
    return sp


def add_row(grid, row, label_text, widget, tooltip=None, bold=False):
    """Adds a label + widget pair to a QGridLayout. The tooltip (if any) is
    applied ONLY to the label, so hovering the parameter name shows the
    explanation while the input field itself stays tooltip-free."""
    lbl = QLabel(label_text)
    if bold:
        f = lbl.font(); f.setBold(True); lbl.setFont(f)
    if tooltip:
        lbl.setToolTip(tooltip)
    grid.addWidget(lbl, row, 0)
    grid.addWidget(widget, row, 1)
    return lbl


# =============================================================================
# MAIN WINDOW
# =============================================================================

# =============================================================================
# STYLESHEETS (QSS) - light mode / dark mode
# =============================================================================
APP_STYLESHEET_LIGHT = """
QMainWindow, QWidget {
    background-color: #F1F3F7;
    color: #1E2430;
    font-family: "Segoe UI", "Cantarell", sans-serif;
    font-size: 10pt;
}

QWidget#HeaderBanner {
    background-color: #F0F4F8; /* Fondo claro (azul/gris muy suave) */
    border-bottom: 1px solid #D1D9E6; /* Optional subtle border to define the boundary */
}

QLabel#HeaderTitle {
    color: #1A2536; /* Azul oscuro/gris para texto principal */
    font-size: 15pt;
    font-weight: 600;
}

QLabel#HeaderSubtitle {
    color: #5A6B82; /* Medium gray for the subtitle */
    font-size: 9pt;
}

QTabWidget::pane {
    border: 1px solid #D4D9E2;
    border-radius: 0px;
    background-color: #FFFFFF;
    top: -1px;
}
QTabBar::tab {
    background: #E4E8F0;
    color: #384154;
    padding: 7px 16px;
    border-top-left-radius: 0px;
    border-top-right-radius: 0px;
    margin-right: 2px;
    font-weight: 600;
}
QTabBar::tab:selected {
    background: #FFFFFF;
    color: #14406E;
    border: 1px solid #D4D9E2;
    border-bottom: none;
}
QTabBar::tab:hover { background: #EFF3FA; }

QGroupBox {
    background-color: #FFFFFF;
    border: 1px solid #DDE2EA;
    border-radius: 0px;
    margin-top: 10px;
    padding: 7px 6px 6px 6px;
    font-weight: 600;
    color: #23324A;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 6px;
    color: #14406E;
    background-color: #FFFFFF;
}

QLabel { color: #2B3242; }

QComboBox, QDoubleSpinBox, QSpinBox {
    background-color: #FFFFFF;
    border: 1px solid #C7CEDA;
    border-radius: 0px;
    padding: 3px 6px;
    min-height: 20px;
    selection-background-color: #1F87E6;
}
QComboBox:hover, QDoubleSpinBox:hover, QSpinBox:hover { border: 1px solid #1F87E6; }
QComboBox:focus, QDoubleSpinBox:focus, QSpinBox:focus { border: 1px solid #14406E; }
QComboBox::drop-down { border: none; width: 18px; }

QPushButton {
    border-radius: 0px;
    padding: 7px 10px;
    font-weight: 600;
    border: none;
}
QPushButton:disabled { background-color: #C6CCD6 !important; color: #8A8F99 !important; }
QPushButton#ComputeBtn { background-color: #1F87E6; color: white; font-size: 10.5pt; }
QPushButton#ComputeBtn:hover { background-color: #1B76C9; }
QPushButton#ComputeBtn:pressed { background-color: #155C9E; }
QPushButton#ExportCadBtn { background-color: #2E9E4C; color: white; }
QPushButton#ExportCadBtn:hover { background-color: #268A41; }
QPushButton#ExportExcelBtn { background-color: #D97318; color: white; }
QPushButton#ExportExcelBtn:hover { background-color: #C1650F; }
QPushButton#SecondaryBtn { background-color: #E4E8F0; color: #23324A; }
QPushButton#SecondaryBtn:hover { background-color: #D6DCE7; }

QGroupBox#ExportGroup { background-color: #F3FAF4; border-color: #CDE9D2; }
QGroupBox#HydrofoilGroup { background-color: #F7F8FF; border-color: #DADEF5; }

QTableWidget {
    background-color: #FFFFFF;
    gridline-color: #E7EAF0;
    border: 1px solid #DDE2EA;
    border-radius: 0px;
    selection-background-color: #DCEBFB;
    selection-color: #14406E;
}
QHeaderView::section {
    background-color: #EEF1F6;
    color: #23324A;
    padding: 5px;
    border: none;
    border-bottom: 1px solid #DDE2EA;
    font-weight: 600;
}

QMenuBar { background-color: #16233D; color: #E7ECF5; padding: 2px; }
QMenuBar::item { background: transparent; padding: 5px 10px; border-radius: 0px; }
QMenuBar::item:selected { background: #24365C; }
QMenu { background-color: #FFFFFF; border: 1px solid #D4D9E2; }
QMenu::item { padding: 6px 22px 6px 14px; }
QMenu::item:selected { background-color: #1F87E6; color: white; }
QMenu::separator { height: 1px; background: #E4E8F0; margin: 4px 6px; }

QStatusBar { background-color: #EEF1F6; border-top: 1px solid #DDE2EA; }
QStatusBar QLabel { color: #4A5468; }

QScrollBar:vertical { background: #F1F3F7; width: 11px; margin: 0; }
QScrollBar::handle:vertical { background: #C7CEDA; border-radius: 0px; min-height: 24px; }
QScrollBar::handle:vertical:hover { background: #ABB4C4; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }

QSplitter::handle { background-color: #E4E8F0; }

QToolTip {
    background-color: #FFFFE9;
    color: #1E2430;
    border: 1px solid #C9C9A8;
    padding: 4px 7px;
}
"""

APP_STYLESHEET_DARK = """
QMainWindow, QWidget {
    background-color: #1B1F2A;
    color: #F0F4F8;
    font-family: "Segoe UI", "Cantarell", sans-serif;
    font-size: 10pt;
}

QLabel#HeaderTitle {
    color: #FFFFFF;
    font-size: 15pt;
    font-weight: 600;
}
QLabel#HeaderSubtitle {
    color: #C2D6EE;
    font-size: 9pt;
}
QWidget#HeaderBanner {
    background-color: #0C1220;
}

QTabWidget::pane {
    border: 1px solid #333B4D;
    border-radius: 0px;
    background-color: #232837;
    top: -1px;
}
QTabBar::tab {
    background: #262C3B;
    color: #E6F0FA;
    padding: 7px 16px;
    border-top-left-radius: 0px;
    border-top-right-radius: 0px;
    margin-right: 2px;
    font-weight: 600;
}
QTabBar::tab:selected {
    background: #232837;
    color: #8BBFEE;
    border: 1px solid #333B4D;
    border-bottom: none;
}
QTabBar::tab:hover { background: #2D3446; }

QGroupBox {
    background-color: #232837;
    border: 1px solid #333B4D;
    border-radius: 0px;
    margin-top: 10px;
    padding: 7px 6px 6px 6px;
    font-weight: 500; /* Reducido de 600 a 500 para suavizar el peso */
    color: #F0F4F8;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 6px;
    color: #8BBFEE; /* Azul mate menos estridente */
    font-weight: 500; /* Softer bold weight */
    background-color: #232837;
}

QLabel { color: #F0F4F8; }

QComboBox, QDoubleSpinBox, QSpinBox {
    background-color: #1B1F2A;
    color: #F0F4F8;
    border: 1px solid #3B4457;
    border-radius: 0px;
    padding: 3px 6px;
    min-height: 20px;
    selection-background-color: #1F87E6;
}
QComboBox:hover, QDoubleSpinBox:hover, QSpinBox:hover { border: 1px solid #1F87E6; }
QComboBox:focus, QDoubleSpinBox:focus, QSpinBox:focus { border: 1px solid #8BBFEE; }
QComboBox::drop-down { border: none; width: 18px; }
QComboBox QAbstractItemView {
    background-color: #232837; color: #F0F4F8; selection-background-color: #1F87E6;
}

QPushButton {
    border-radius: 0px;
    padding: 7px 10px;
    font-weight: 600;
    border: none;
}
QPushButton:disabled { background-color: #3A4152 !important; color: #8C9BB0 !important; }
QPushButton#ComputeBtn { background-color: #1F87E6; color: #FFFFFF; font-size: 10.5pt; }
QPushButton#ComputeBtn:hover { background-color: #2E96F5; }
QPushButton#ComputeBtn:pressed { background-color: #155C9E; }
QPushButton#ExportCadBtn { background-color: #2E9E4C; color: #FFFFFF; }
QPushButton#ExportCadBtn:hover { background-color: #37B058; }
QPushButton#ExportExcelBtn { background-color: #D97318; color: #FFFFFF; }
QPushButton#ExportExcelBtn:hover { background-color: #EA8226; }
QPushButton#SecondaryBtn { background-color: #2D3446; color: #F0F4F8; }
QPushButton#SecondaryBtn:hover { background-color: #3A4359; }

QGroupBox#ExportGroup { background-color: #1A2333; border-color: #2B3D5B; }
QGroupBox#HydrofoilGroup { background-color: #1D2130; border-color: #333B58; }

QTableWidget {
    background-color: #232837;
    color: #F0F4F8;
    gridline-color: #333B4D;
    border: 1px solid #333B4D;
    border-radius: 0px;
    selection-background-color: #234873;
    selection-color: #FFFFFF;
}
QHeaderView::section {
    background-color: #262C3B;
    color: #F0F4F8;
    padding: 5px;
    border: none;
    border-bottom: 1px solid #333B4D;
    font-weight: 600;
}

QMenuBar { background-color: #0C1220; color: #F0F4F8; padding: 2px; }
QMenuBar::item { background: transparent; padding: 5px 10px; border-radius: 0px; }
QMenuBar::item:selected { background: #232E47; }
QMenu { background-color: #232837; color: #F0F4F8; border: 1px solid #333B4D; }
QMenu::item { padding: 6px 22px 6px 14px; }
QMenu::item:selected { background-color: #1F87E6; color: #FFFFFF; }
QMenu::separator { height: 1px; background: #333B4D; margin: 4px 6px; }

QStatusBar { background-color: #1B1F2A; border-top: 1px solid #333B4D; }
QStatusBar QLabel { color: #DCE8F5; }

QScrollBar:vertical { background: #1B1F2A; width: 11px; margin: 0; }
QScrollBar::handle:vertical { background: #3B4457; border-radius: 0px; min-height: 24px; }
QScrollBar::handle:vertical:hover { background: #4C5773; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }

QSplitter::handle { background-color: #333B4D; }

QToolTip {
    background-color: #2A3040;
    color: #F0F4F8;
    border: 1px solid #4A5468;
    padding: 4px 7px;
}
"""


class KaplanDeriaz3DApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("KaplanDeriaz3D - Hydraulic Turbine Blade Designer (Professional Edition)")
        self.resize(1720, 1020)

        # --- State / computed geometry ---
        self.is_computed = False
        self.X_mid = self.Y_mid = self.Z_mid = None
        self.X_sol_ext = self.Y_sol_ext = self.Z_sol_ext = None
        self.X_sol_int = self.Y_sol_int = self.Z_sol_int = None
        self.radius_vec = self.beta1_deg = self.beta2_deg = None
        self.Z_optimo_current = None
        self.R_color = None
        self.last_is_kaplan = True
        self.last_backbone = None
        self._has_rendered_once = False
        self.dark_mode = False

        self._build_menu_bar()
        self._build_ui()
        self._build_status_bar()
        self.apply_theme(dark=True)
        self.on_turbine_type_change()
        self.on_profile_type_change()

    # -------------------------------------------------------------------
    # MENU BAR (File / View / Help)
    # -------------------------------------------------------------------
    def _build_menu_bar(self):
        menubar = self.menuBar()

        menu_file = menubar.addMenu("&File")
        act_save_cfg = QAction("Save Configuration...", self)
        act_save_cfg.triggered.connect(self.save_configuration)
        menu_file.addAction(act_save_cfg)
        act_load_cfg = QAction("Load Configuration...", self)
        act_load_cfg.triggered.connect(self.load_configuration)
        menu_file.addAction(act_load_cfg)
        menu_file.addSeparator()
        act_reset_defaults = QAction("Reset to Defaults", self)
        act_reset_defaults.triggered.connect(self.reset_to_defaults)
        menu_file.addAction(act_reset_defaults)
        menu_file.addSeparator()
        act_export_cad = QAction("Export CAD Geometry...", self)
        act_export_cad.triggered.connect(self.export_cad)
        menu_file.addAction(act_export_cad)
        act_export_xlsx = QAction("Export Blade Angles (Excel)...", self)
        act_export_xlsx.triggered.connect(self.export_excel)
        menu_file.addAction(act_export_xlsx)
        menu_file.addSeparator()
        act_close = QAction("Close", self)
        act_close.triggered.connect(self.close)
        menu_file.addAction(act_close)

        menu_view = menubar.addMenu("&View")
        act_reset_view = QAction("Reset 3D View", self)
        act_reset_view.triggered.connect(self.reset_3d_view)
        menu_view.addAction(act_reset_view)
        act_copy_table = QAction("Copy Results Table", self)
        act_copy_table.triggered.connect(self.copy_results_to_clipboard)
        menu_view.addAction(act_copy_table)
        menu_view.addSeparator()

        self.act_light_mode = QAction("Light Mode", self, checkable=True)
        self.act_light_mode.triggered.connect(lambda: self.apply_theme(dark=False))
        menu_view.addAction(self.act_light_mode)
        self.act_dark_mode = QAction("Dark Mode", self, checkable=True)
        self.act_dark_mode.setChecked(True)
        self.act_dark_mode.triggered.connect(lambda: self.apply_theme(dark=True))
        menu_view.addAction(self.act_dark_mode)

        menu_help = menubar.addMenu("&Help")
        act_about = QAction("About KaplanDeriaz3D...", self)
        act_about.triggered.connect(self.show_about)
        menu_help.addAction(act_about)

    # -------------------------------------------------------------------
    # TEMA (Light / Dark)
    # -------------------------------------------------------------------
    def apply_theme(self, dark):
        self.dark_mode = dark
        self.setStyleSheet(APP_STYLESHEET_DARK if dark else APP_STYLESHEET_LIGHT)
        self.act_dark_mode.setChecked(dark)
        self.act_light_mode.setChecked(not dark)

        panel_color = "#232837" if dark else "#FFFFFF"
        border_color = "#333B4D" if dark else "#DDE2EA"
        
        # Update the outer render container
        self.renders_container.setStyleSheet(
            f"background-color: {panel_color}; border: 1px solid {border_color}; border-radius: 0px;")

        # --- Always refresh the 3D viewport backgrounds when the theme changes ---
        if hasattr(self, 'plotter_mid') and hasattr(self, 'plotter_solid'):
            self.plotter_mid.set_background(self._plot_bg_color())
            self.plotter_solid.set_background(self._plot_bg_color())

        # If geometry has already been computed, re-render so the 3D
        # viewers (background, text, colorbar) pick up the new theme.
        if self.is_computed and self.last_backbone is not None:
            self._render(self.last_backbone, self.last_is_kaplan)
        else:
            # Forces a visual refresh even if the scene is still empty
            self.plotter_mid.render()
            self.plotter_solid.render()

    def _plot_bg_color(self):
        return '#232837' if self.dark_mode else 'white'

    def _plot_text_color(self):
        return 'white' if self.dark_mode else 'black'

    # -------------------------------------------------------------------
    # STATUS BAR
    # -------------------------------------------------------------------
    def _build_status_bar(self):
        status = QStatusBar()
        self.setStatusBar(status)
        self.status_label = QLabel("Ready. Set your parameters and press COMPUTE.")
        self.last_computed_label = QLabel("")
        status.addWidget(self.status_label, stretch=1)
        status.addPermanentWidget(self.last_computed_label)

    def _set_status(self, text, kind='info'):
        colors = {'info': '#4A5468', 'busy': '#1B76C9', 'ok': '#1C7C33', 'error': '#B3261E'}
        if self.dark_mode:
            colors = {'info': '#ABB2C4', 'busy': '#6EB3FF', 'ok': '#59D67C', 'error': '#FF8A80'}
        self.status_label.setStyleSheet(f"color: {colors.get(kind, colors['info'])}; font-weight: 600;")
        self.status_label.setText(text)

    # -------------------------------------------------------------------
    # UI CONSTRUCTION
    # -------------------------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        outer_layout = QVBoxLayout(central)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        # ---- Cabecera (banner) ----
        # QFrame is used instead of QWidget: QFrame reliably paints its
        # background from the stylesheet across platforms/Windows styles,
        # whereas a plain QWidget can sometimes fail to paint the QSS
        # "background-color" (even with WA_StyledBackground enabled)
        # depending on the active system style.
        header = QFrame()
        header.setObjectName("HeaderBanner")
        header.setAttribute(Qt.WA_StyledBackground, True)
        header.setFrameShape(QFrame.NoFrame)
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(18, 8, 18, 8)
        header_layout.setSpacing(1)
        title_lbl = QLabel("KaplanDeriaz3D - Hydraulic Turbine Blade Designer")
        title_lbl.setObjectName("HeaderTitle")
        subtitle_lbl = QLabel(
            "Solid Turbine Designer with CAD / Excel Export - Kaplan (Axial) & Deriaz (Diagonal) Runners")
        subtitle_lbl.setObjectName("HeaderSubtitle")
        header_layout.addWidget(title_lbl)
        header_layout.addWidget(subtitle_lbl)
        outer_layout.addWidget(header)

        content = QWidget()
        main_layout = QHBoxLayout(content)
        main_layout.setContentsMargins(8, 8, 8, 6)
        outer_layout.addWidget(content, stretch=1)

        # ---- Left column: tabs (fixed width, not draggable) ----
        # A generous fixed width is used instead of a QSplitter so the
        # parameters panel can never be accidentally narrowed and truncate
        # texto (p.ej. "Chord Stations (dz/dgamma/ds):").
        self.tabs = QTabWidget()
        self.tabs.setFixedWidth(520)
        main_layout.addWidget(self.tabs)

        self._build_tab_general()
        self._build_tab_hydro()
        self._build_tab_export()

        # ---- Right column: PyVista renders + results table ----
        right = QWidget()
        right_layout = QVBoxLayout(right)
        main_layout.addWidget(right, stretch=1)

        # Small toolbar above the renders (title + Reset View button)
        renders_toolbar = QHBoxLayout()
        renders_title = QLabel("3D Visualization")
        f = renders_title.font(); f.setBold(True); f.setPointSize(f.pointSize() + 1)
        renders_title.setFont(f)
        renders_toolbar.addWidget(renders_title)
        renders_toolbar.addStretch()
        self.btn_reset_view = QPushButton(" Reset 3D View")
        self.btn_reset_view.setObjectName("SecondaryBtn")
        self.btn_reset_view.setToolTip(
            "Restores the default camera position on both 3D viewports\n"
            "without recomputing any geometry.")
        self.btn_reset_view.clicked.connect(self.reset_3d_view)
        renders_toolbar.addWidget(self.btn_reset_view)
        right_layout.addLayout(renders_toolbar)

        # Horizontal container for the two PyVista 3D viewers
        self.renders_container = QFrame()
        self.renders_container.setAttribute(Qt.WA_StyledBackground, True)
        self.renders_container.setFrameShape(QFrame.NoFrame)
        renders_layout = QHBoxLayout(self.renders_container)
        renders_layout.setContentsMargins(0, 0, 0, 0)

        # Viewer 1: Mean surface
        self.plotter_mid = QtInteractor(self.renders_container)
        self.plotter_mid.set_background(self._plot_bg_color())
        renders_layout.addWidget(self.plotter_mid)

        # Viewer 2: Solid blade profile
        self.plotter_solid = QtInteractor(self.renders_container)
        self.plotter_solid.set_background(self._plot_bg_color())
        renders_layout.addWidget(self.plotter_solid)

        # Constrain mouse rotation to azimuth (left/right) and elevation
        # (up/down) only, with a fixed "up" vector (world Z axis) - this
        # prevents the camera from ever tilting/rolling the horizon.
        # Note: the exact keyword arguments accepted by enable_terrain_style()
        # differ across pyvista versions, so call it defensively.
        for plotter in (self.plotter_mid, self.plotter_solid):
            try:
                plotter.enable_terrain_style(mouse_wheel_zooming=True)
            except TypeError:
                plotter.enable_terrain_style()
            try:
                self._set_default_camera_view(plotter)
            except Exception:
                pass  # Empty scene; camera will be set properly on first compute.

        # Cap the render area height so it doesn't crowd out the results table
        self.renders_container.setMaximumHeight(500)
        self.renders_container.setStyleSheet(
            f"background-color: {self._plot_bg_color()}; border: 1px solid #333B4D; border-radius: 0px;")
        right_layout.addWidget(self.renders_container, stretch=3)


        # Bottom results table
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(['Analyzed Parameter', 'Computed Value'])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.verticalHeader().setDefaultSectionSize(24)
        self.table.verticalHeader().setVisible(False)
        
        # --- HEIGHT AND SCROLL ADJUSTMENT ---
        self.table.setMinimumHeight(190)
        # Allow the table to grow a bit more if the window is tall enough
        self.table.setMaximumHeight(350) 
        
        # Forces the vertical scrollbar to appear once rows exceed the available space
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        right_layout.addWidget(self.table, stretch=1)

    # -------------------------------------------------------------------
    # TAB 1: GENERAL PARAMETERS
    # -------------------------------------------------------------------
    def _build_tab_general(self):
        tab = QWidget()
        self.tabs.addTab(tab, " General ")
        outer = QVBoxLayout(tab)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        container = QWidget()
        scroll.setWidget(container)
        layout = QVBoxLayout(container)
        layout.setSpacing(6)
        self._general_container = container  # needed to force a relayout when
        # switching turbine type (Kaplan/Deriaz geometry panels have a
        # different number of rows), so the Compute button and everything
        # below it stays reachable via the scrollbar instead of getting
        # visually stuck off the bottom of the visible area.

        # --- Group: Design Parameters ---
        group = QGroupBox(" Design Parameters ")
        grid = QGridLayout(group)
        grid.setVerticalSpacing(4)

        r = 0
        self.combo_type = QComboBox()
        self.combo_type.addItems(['Kaplan (Axial)', 'Deriaz (Diagonal)'])
        self.combo_type.currentIndexChanged.connect(self.on_turbine_type_change)
        add_row(grid, r, 'Turbine Type:', self.combo_type, bold=False, tooltip=(
            "Hydraulic Turbine Type:\n"
            "- Kaplan (Axial): Purely axial flow. Optimal for low heads (Hn < 40 m) and large flow rates.\n"
            "- Deriaz (Diagonal): Diagonal/mixed flow. Ideal for medium heads (40 m < Hn < 200 m) "
            "and reversible pumped storage."))

        r += 1
        self.spin_rpm = make_double_spin(450, 1, 10000, decimals=1, step=10)
        add_row(grid, r, 'Rotational Speed (RPM):', self.spin_rpm, tooltip=(
            "Runner Rotational Speed (N in RPM):\n"
            "- Sets the angular speed of the shaft (\u03C9 = 2\u03C0N/60).\n"
            "- Directly impacts the inlet and outlet velocity triangles by altering the peripheral "
            "speed (u = \u03C9\u00B7r).\n"
            "- Must synchronize with the grid frequency according to the generator pole pairs."))

        r += 1
        self.combo_rot = QComboBox()
        self.combo_rot.addItems(['Counter-Clockwise (Standard)', 'Clockwise'])
        add_row(grid, r, 'Rotation Direction:', self.combo_rot, tooltip=(
            "Runner Direction of Rotation:\n"
            "- Counter-Clockwise (Standard): Standard for most turbines when viewed from above.\n"
            "- Clockwise: Inverts blade kinematic orientation and tangential flow components."))

        r += 1
        self.spin_q0 = make_double_spin(12, 0.01, 10000, decimals=2, step=0.5)
        add_row(grid, r, 'Flow Rate Q0 (m3/s):', self.spin_q0, tooltip=(
            "Nominal Volumetric Flow Rate (Q0):\n"
            "- Total volume of water passing through the runner per second.\n"
            "- Along with the flow area, it determines the meridian flow velocity "
            "(Vm = Q_effective / A_flow), fixing the streamline slopes."))

        r += 1
        self.spin_hn = make_double_spin(15, 0.1, 2000, decimals=2, step=1)
        add_row(grid, r, 'Net Head Hn (m):', self.spin_hn, tooltip=(
            "Net Head / Effective Head (Hn):\n"
            "- Specific energy available per unit weight of fluid between turbine inlet and outlet.\n"
            "- Governs total hydraulic power available (Ph = \u03C1\u00B7g\u00B7Q0\u00B7Hn) and the "
            "tangential velocity V\u03B8 for all streamlines."))

        r += 1
        self.spin_g = make_double_spin(9.81, 1.0, 20.0, decimals=2, step=0.01)
        add_row(grid, r, 'Gravity g (m/s2):', self.spin_g,
                tooltip="Local gravitational acceleration (typically 9.81 m/s^2).")
        
        r += 1
        self.spin_rho = make_double_spin(1000, 1.0, 10000.0, decimals=2, step=1)
        add_row(grid, r, 'Water density rho (kg/m^3):', self.spin_rho,
                tooltip="Water liquid density (typically 1000 kg/m^3).")

        r += 1
        self.spin_eta_h = make_double_spin(0.90, 0.1, 1.0, decimals=2, step=0.01)
        add_row(grid, r, 'Hydraulic Eff. (eta_h):', self.spin_eta_h, tooltip=(
            "Hydraulic Efficiency (\u03B7_h):\n"
            "- Evaluates head losses due to friction, boundary layer separation, and inlet shock.\n"
            "- Modifies effective Eulerian head transferred to the shaft: H_euler = Hn \u00B7 \u03B7_h."))

        r += 1
        self.spin_eta_v = make_double_spin(0.96, 0.1, 1.0, decimals=2, step=0.01)
        add_row(grid, r, 'Volumetric Eff. (eta_v):', self.spin_eta_v, tooltip=(
            "Volumetric Efficiency (\u03B7_v):\n"
            "- Quantifies leakage losses through peripheral tip/hub clearances.\n"
            "- Effective flow rate performing useful work: Q_effective = Q0 \u00B7 \u03B7_v."))

        r += 1
        self.spin_eta_o = make_double_spin(0.98, 0.1, 1.0, decimals=2, step=0.01)
        add_row(grid, r, 'Mechanical Eff. (eta_m):', self.spin_eta_o, tooltip=(
            "Mechanical Efficiency (\u03B7_m):\n"
            "- Accounts for mechanical friction losses in bearings, seals, and shaft components.\n"
            "- Total Turbine Efficiency: \u03B7_total = \u03B7_h \u00B7 \u03B7_v \u00B7 \u03B7_m."))

        r += 1
        self.spin_sigma = make_double_spin(1.25, 0.5, 3.0, decimals=2, step=0.05)
        add_row(grid, r, 'Solidity (sigma):', self.spin_sigma, tooltip=(
            "Runner Solidity (\u03C3 = Lc / t):\n"
            "- Ratio between mean blade chord (Lc) and tangential pitch (t = 2\u03C0r / Z).\n"
            "- DIRECT BLADE IMPACT: Higher solidity (\u03C3) requires more blades (Z) or longer chords "
            "to distribute hydrodynamic loading.\n"
            "- High values (\u03C3 > 1.3) mitigate cavitation risk but increase viscous friction losses."))

        r += 1
        self.combo_interp = QComboBox()
        self.combo_interp.addItems([
            'Cubic (Standard)', 'Linear (Uniform)', 'Cosine (Smooth)',
            'Inlet Loaded (Attack)', 'Outlet Loaded (Discharge)'])
        add_row(grid, r, 'Interpolation Scheme:', self.combo_interp, tooltip=(
            "Blade-Particle Beta Angle Distribution Law (Pressure distribution along chord):\n"
            "- Cubic (Standard): Smooth distribution. Peak hydrodynamic pressure located around "
            "30-40% of chord length.\n"
            "- Linear (Uniform): Constant angle gradient. Pressure is distributed uniformly along "
            "the profile.\n"
            "- Cosine (Smooth): Ultra-smooth transition near leading/trailing edges. Minimizes "
            "localized cavitation spikes.\n"
            "- Inlet Loaded: Steeper deflection at the entry. Maximum pressure zone at the leading "
            "edge (ideal for clean, shock-free flow entry).\n"
            "- Outlet Loaded: Steeper curvature towards the exit. Maximum pressure zone at the "
            "trailing edge (maximizes energy transfer before discharge)."))

        layout.addWidget(group)

        # --- Group: Kaplan Geometry ---
        self.group_kaplan = QGroupBox("Kaplan Geometry")
        kg = QGridLayout(self.group_kaplan)
        self.spin_rhub = make_double_spin(0.30, 0.01, 50, decimals=3, step=0.01)
        self.spin_rtip = make_double_spin(0.65, 0.01, 50, decimals=3, step=0.01)
        self.spin_lz = make_double_spin(0.35, 0.01, 20, decimals=3, step=0.01)
        add_row(kg, 0, 'Hub Radius R_hub (m):', self.spin_rhub,
                tooltip="Hub / Central core radius of the Kaplan turbine in meters.")
        add_row(kg, 1, 'Tip Radius R_tip (m):', self.spin_rtip,
                tooltip="Outer peripheral radius (Blade Tip) in meters.")
        add_row(kg, 2, 'Axial Length L_z (m):', self.spin_lz,
                tooltip="Projected longitudinal/axial length of the runner along the Z-axis.")
        layout.addWidget(self.group_kaplan)

        # --- Group: Deriaz Geometry ---
        self.group_deriaz = QGroupBox("Deriaz Geometry")
        dg = QGridLayout(self.group_deriaz)
        self.spin_re_int = make_double_spin(2.00, 0.01, 100, decimals=3, step=0.05)
        self.spin_re_ext = make_double_spin(2.60, 0.01, 100, decimals=3, step=0.05)
        self.spin_gamma1 = make_double_spin(30.0, 0.0, 89.0, decimals=1, step=1)
        self.spin_gamma2 = make_double_spin(60.0, 0.0, 89.0, decimals=1, step=1)
        add_row(dg, 0, 'Inner Sph. Rad. Re_int:', self.spin_re_int,
                tooltip="Inner spherical radius for the Deriaz hub dome.")
        add_row(dg, 1, 'Outer Sph. Rad. Re_ext:', self.spin_re_ext,
                tooltip="Outer spherical radius for the diagonal outer casing.")
        add_row(dg, 2, 'Inlet Angle gamma1 (deg):', self.spin_gamma1, tooltip=(
            "Inlet Cone Angle (\u03B31 in degrees):\n"
            "- Angle of the surface of revolution at the leading edge relative to the axis of rotation.\n"
            "- Alters 3D spatial orientation and meridian velocity vector component at inlet."))
        add_row(dg, 3, 'Outlet Angle gamma2 (deg):', self.spin_gamma2, tooltip=(
            "Outlet Cone Angle (\u03B32 in degrees):\n"
            "- Slope angle of the flow surface at the trailing edge relative to the axis.\n"
            "- Controls outlet flow divergence entering the draft tube."))
        layout.addWidget(self.group_deriaz)
        layout.addStretch()

        # --- Compute button ---
        # Deliberately placed OUTSIDE the QScrollArea (added to `outer`, the
        # tab's own layout, not to the scrollable `layout`) so it is always
        # fully visible regardless of how tall the parameter groups above it
        # get (Kaplan has 3 geometry rows, Deriaz has 4) - this avoids any
        # dependency on scroll-area relayout timing.
        self.btn_compute = QPushButton("COMPUTE SOLID 3D BLADE DESIGN")
        self.btn_compute.setObjectName("ComputeBtn")
        self.btn_compute.setFont(QFont("Segoe UI", 10, QFont.Bold))
        self.btn_compute.setToolTip("Runs the full hydraulic + geometric solver and regenerates\nthe 3D blade preview with the current parameters (Ctrl+Enter).")
        self.btn_compute.clicked.connect(self.compute_turbine)
        self.btn_compute.setShortcut(QKeySequence("Ctrl+Return"))
        outer.addWidget(self.btn_compute)

    # -------------------------------------------------------------------
    # TAB 2: HYDROFOIL PROFILE AND MESH
    # -------------------------------------------------------------------
    def _build_tab_hydro(self):
        tab = QWidget()
        self.tabs.addTab(tab, " Hydrofoil / mesh ")
        outer = QVBoxLayout(tab)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        container = QWidget()
        scroll.setWidget(container)
        layout = QVBoxLayout(container)

        # --- Group: Hydrofoil Profile Design ---
        group = QGroupBox("Hydrofoil Profile Design")
        group.setObjectName("HydrofoilGroup")
        ag = QGridLayout(group)

        self.combo_naca = QComboBox()
        self.combo_naca.addItems([
            'NACA 00XX (Standard Symmetric)',
            'Customized (4-Digit NACA Series)',
            'Reversible Hydrofoil (Pump-Turbine)',
            'Anti-Cavitation (Flat Pressure)',
            'Low-Torque S-Camber'])
        self.combo_naca.setCurrentText('Reversible Hydrofoil (Pump-Turbine)')
        self.combo_naca.currentIndexChanged.connect(self.on_profile_type_change)
        add_row(ag, 0, 'Hydrofoil Profile:', self.combo_naca, bold=False, tooltip=(
            "Hydrofoil Formulations:\n"
            "1. NACA 00XX: Uncambered (yc=0) | Standard NACA 4-digit thickness yt(x).\n"
            "2. Customized: Parametric m & p interpolation from hub to tip "
            "(m_factor = (1-r*)^1.5). Piecewise parabolic yc(x) with standard yt(x).\n"
            "3. Reversible: Pure symmetric ellipse (yc=0) | yt = t_c\u00B7\u221A(x\u00B7(1-x)).\n"
            "4. Anti-Cavitation: High forward camber | yt = 2.6896\u00B7t_c\u00B7x\u00B7(1-x)^1.5, "
            "yc = (yt/2)\u00B7m_factor.\n"
            "5. Low-Torque S-Camber: Reflexed camber | yc = 0.8\u00B7x\u00B7(1-x)\u00B7(0.5-x)\u00B7m_factor.\n\n"
            "Note: Relative thickness (t_c) scales with the Hub-to-Tip Ratio."))

        self.spin_thickness = make_double_spin(0.085, 0.01, 0.40, decimals=4, step=0.005)
        add_row(ag, 1, 'Max Rel. Thickness (t/c):', self.spin_thickness, tooltip=(
            "Maximum Relative Airfoil Thickness (t/c):\n"
            "- Non-dimensional ratio between max thickness (t) and chord (c).\n"
            "- Higher thickness increases structural stiffness at the hub, but raises drag and "
            "flow separation risk."))

        self.spin_hub_to_tip = make_double_spin(0.65, 0.10, 1.00, decimals=2, step=0.05)
        add_row(ag, 2, 'Hub to Tip Thickness Ratio:', self.spin_hub_to_tip, tooltip=(
            "Hub to Tip Relative Thickness Ratio (t_tip / t_hub):\n"
            "- Ratio between the profile normalized relative thickness at the blade tip and at "
            "the hub/root.\n"
            "- A value of 0.65 means the tip retains 65% of the hub relative thickness "
            "(35% taper loss).\n"
            "- Reduces structural stress from centrifugal forces and improves hydrodynamic "
            "performance near the tip."))

        # Sub-panel: Customized
        self.group_custom = QGroupBox("Customized (4-Digit NACA Series)")
        self.group_custom.setVisible(False)
        cg = QGridLayout(self.group_custom)
        self.spin_m_hub = make_double_spin(2.0, -9.0, 9.0, decimals=1, step=0.5)
        self.spin_p_hub = make_double_spin(4.0, 1.0, 9.0, decimals=1, step=0.5)
        self.spin_m_tip = make_double_spin(0.0, -9.0, 9.0, decimals=1, step=0.5)
        self.spin_p_tip = make_double_spin(4.0, 1.0, 9.0, decimals=1, step=0.5)
        lbl_m_hub = QLabel('m_hub (%):'); lbl_m_hub.setToolTip("Max camber at hub as % of chord (e.g. 2 for 2%).")
        lbl_p_hub = QLabel('p_hub (x10%):'); lbl_p_hub.setToolTip("Location of max camber at hub (e.g. 4 for 40% chord).")
        lbl_m_tip = QLabel('m_tip (%):'); lbl_m_tip.setToolTip("Max camber at tip as % of chord (0 = symmetric).")
        lbl_p_tip = QLabel('p_tip (x10%):'); lbl_p_tip.setToolTip("Location of max camber at tip (e.g. 4 for 40% chord).")
        cg.addWidget(lbl_m_hub, 0, 0); cg.addWidget(self.spin_m_hub, 0, 1)
        cg.addWidget(lbl_p_hub, 0, 2); cg.addWidget(self.spin_p_hub, 0, 3)
        cg.addWidget(lbl_m_tip, 1, 0); cg.addWidget(self.spin_m_tip, 1, 1)
        cg.addWidget(lbl_p_tip, 1, 2); cg.addWidget(self.spin_p_tip, 1, 3)
        ag.addWidget(self.group_custom, 3, 0, 1, 2)

        layout.addWidget(group)

        # --- Group: Advanced Options ---
        group_adv = QGroupBox("Advanced Options (Mesh Resolution)")
        advg = QGridLayout(group_adv)
        self.spin_nradios = make_int_spin(100, 5, 200)
        self.spin_ncuerda = make_int_spin(200, 15, 500)
        add_row(advg, 0, 'Streamlines (N_radii):', self.spin_nradios, tooltip=(
            "Spanwise Mesh Resolution (N_radii):\n"
            "- Number of evaluated radial streamlines from Hub to Tip."))
        add_row(advg, 1, 'Chord Stations (dz/dgamma/ds):', self.spin_ncuerda, tooltip=(
            "Chordwise Profile Resolution (dz/dgamma/ds):\n"
            "- Discretization points along the chord (LE to TE)."))

        # --- Industrial post-processing adjustments (optional, off by default) ---
        self.check_industrial_postproc = QCheckBox("Enable Industrial Post-Processing")
        self.check_industrial_postproc.setChecked(False)
        self.check_industrial_postproc.setToolTip(
            "Applies a set of geometric adjustments to the mean-surface backbone that "
            "approximate real CAD/manufacturing finishing operations performed on a computed "
            "'ideal' hydrodynamic runner before it can be built (leading-edge sweep, tip/root "
            "corner rounding, and - for Kaplan only - trailing-edge radial alignment). All are "
            "applied to the mean line BEFORE the hydrofoil thickness is generated, so the profile "
            "always tapers correctly on the adjusted geometry.\n\n"
            "IMPORTANT: enabling this will inevitably move the geometry away from the "
            "hydraulically optimal shape computed by the inverse design core, in exchange for a "
            "constructible, manufacturing-ready runner. Each individual adjustment below states "
            "its own hydraulic cost (or confirms it has none) in its own tooltip.")
        self.check_industrial_postproc.toggled.connect(self.on_industrial_postproc_toggled)
        advg.addWidget(self.check_industrial_postproc, 2, 0, 1, 2)

        self.group_industrial = QGroupBox("Industrial Post-Processing Adjustments")
        self.group_industrial.setVisible(False)
        ig = QGridLayout(self.group_industrial)

        # a) Trailing-edge radial alignment. Available for BOTH turbine
        # types; see the tooltip below and the align_trailing_edge_radially()
        # docstring for the "no free lunch" framing - neither leaving this
        # off nor turning it on is simply "the correct" choice, they trade
        # different kinds of cost against each other.
        # a) Leading-edge sweep - KAPLAN ONLY. Verified numerically that a
        # uniform per-streamline Z-offset pulls every Deriaz streamline
        # (other than the hub) measurably off its design sphere, which
        # breaks the physical spherical hub/casing compatibility the
        # Deriaz double-regulation mechanism depends on. See
        # apply_leading_edge_sweep()'s docstring for the full derivation.
        # Placed first because it runs first in the post-processing
        # pipeline (see compute_turbine()): a pure translation, so it
        # commutes cleanly with everything that follows.
        self.check_le_sweep = QCheckBox("Apply Leading-Edge Sweep (Axial Translation)")
        self.check_le_sweep.setChecked(False)
        self.check_le_sweep.setToolTip(
            "KAPLAN ONLY (disabled for Deriaz - see below). Rigid axial (Z) translation of each "
            "streamline, growing linearly with radial distance from the hub at a rate set by the "
            "sweep angle below (the standard aerospace/turbomachinery sweep-angle definition). "
            "This is the ONLY adjustment in this panel with NO hydraulic cost on a Kaplan "
            "backbone: since the whole streamline (leading and trailing edge together) is "
            "shifted by the same constant offset, neither its chord length, nor its local beta "
            "angle, nor its total swept angle change at all - only its position in space.\n\n"
            "WHY THIS IS DISABLED FOR DERIAZ: a Kaplan streamline lives on a cylinder of constant "
            "radius, which is invariant under a Z-translation along its own axis, so the trick "
            "above is exact. A Deriaz streamline instead lives on a SPHERE of constant radius Re; "
            "a uniform Z-offset does NOT preserve that sphere (verified numerically - every "
            "streamline except the one left with zero offset drifts measurably off its design "
            "radius), which would break the physical mechanism the Deriaz double-regulation "
            "depends on. There is currently no equivalent sweep operation offered for Deriaz.")
        ig.addWidget(self.check_le_sweep, 0, 0, 1, 2)
        self.spin_sweep_angle = make_double_spin(10.0, 0.0, 45.0, decimals=1, step=1.0)
        add_row(ig, 1, 'Sweep Angle (deg):', self.spin_sweep_angle, tooltip=(
            "Sweep angle, in the standard aerospace/turbomachinery sense: "
            "tan(angle) = axial offset / spanwise distance from the hub. "
            "10 degrees is an illustrative starting value, not a verified industry-normalized "
            "figure for Kaplan runners specifically - adjust freely."))

        # b) Trailing-edge alignment (wake trim), via an explicit plane
        # cut - see align_trailing_edge_radially()'s docstring for the
        # full geometric construction (two lines through the axis: the
        # axis itself, and the line to the hub's own natural trailing-
        # edge point). Runs SECOND in the pipeline, after the sweep above
        # (order-independent, since sweep never changes theta) and BEFORE
        # corner rounding (order that DOES matter - see the docstring's
        # explanation of why rounding must come after this, not before).
        self.check_align_trailing_edge = QCheckBox("Align Trailing Edge (Plane Cut)")
        self.check_align_trailing_edge.setChecked(False)
        self.check_align_trailing_edge.setToolTip(
            "Cuts the blade with a single flat plane and keeps only the portion behind it. The "
            "plane is built from two lines that both pass through the rotation axis: the axis "
            "itself, and the line from the origin to the HUB's own natural trailing-edge point "
            "(for Kaplan: at r=R_hub, z=-L_z; for Deriaz: at Re=Re_int, gamma=gamma2, i.e. the "
            "same anchor point used elsewhere for hub-referenced geometry). Two lines through a "
            "common point define a plane; because both lines pass through the axis, this plane is "
            "a flat cut at constant circumferential angle theta, at exactly the hub's own natural "
            "discharge angle. No free parameter: the hub is always the anchor, by construction.\n\n"
            "NEITHER OPTION IS FREE - THIS IS A GENUINE TRADE-OFF BETWEEN TWO DIFFERENT COSTS, "
            "NOT A 'CORRECT' CHOICE VS. A NECESSARY EVIL:\n\n"
            "Leaving this UNCHECKED keeps the natural, hydrodynamically pure geometry (every "
            "streamline integrated exactly to its Euler-equation design point, with zero "
            "artificial energy loss). Its cost is geometric/structural: the trailing edge lands "
            "at a different circumferential angle at each radius (theta varies significantly "
            "across the span - often tens of degrees between hub and tip), which can concentrate "
            "stress at the discharge edge, promote vibration, and leave an irregular wake shape "
            "downstream.\n\n"
            "CHECKING this trades that geometric cost for an energetic one: every streamline other "
            "than the hub is cut off before its integration reaches the point where the exit "
            "tangential velocity was meant to reach zero - so it delivers less specific energy "
            "than designed and discharges with residual swirl, a possible seed for downstream "
            "vortices and pressure pulsations. Verified to stay exactly on each streamline's own "
            "design surface (cylinder for Kaplan, sphere for Deriaz) - this cut only ever trims "
            "arc length along the already-computed curve, it never displaces a point off of it.\n\n"
            "Which trade-off is preferable depends on your priorities for this specific design "
            "(peak theoretical efficiency vs. discharge-edge mechanical cleanliness) - there is "
            "no universally 'correct' default here.")
        ig.addWidget(self.check_align_trailing_edge, 2, 0, 1, 2)

        # c) Corner rounding - FOUR independent corners (leading-edge/tip,
        # trailing-edge/tip, leading-edge/root, trailing-edge/root), each
        # with its own on/off switch and its own fillet extent. An
        # earlier version of this only offered "tip" and "root" as
        # combined pairs (both corners at that end rounded together with
        # the same extent) - see apply_corner_rounding()'s docstring for
        # the geometry, now generalized to four fully independent corners.
        self.check_round_le_tip = QCheckBox("Round Tip Corner - Leading Edge Side")
        self.check_round_le_tip.setChecked(True)
        self.check_round_le_tip.setToolTip(
            "Rounds the corner of the blade planform where the LEADING edge meets the tip, "
            "independently of the other three corners. Part of the anti-cavitation lip that many "
            "commercial Kaplan runners carry on the suction side near the periphery, to mitigate "
            "tip-clearance/tip-vortex cavitation erosion. Applied to the mean line before the "
            "hydrofoil thickness is generated, using a true 2D corner-fillet construction scaled "
            "to the blade's own reference chord (see apply_corner_rounding() in the source for "
            "the exact geometry).")
        ig.addWidget(self.check_round_le_tip, 3, 0, 1, 2)
        self.spin_le_tip_round_extent = make_double_spin(0.15, 0.0, 1.00, decimals=2, step=0.01)
        add_row(ig, 4, 'LE/Tip Rounding Extent (% of ref. chord):', self.spin_le_tip_round_extent, tooltip=(
            "Fillet radius for this corner, as a fraction of the blade's reference chord length "
            "(not of the spanwise extent). A safety floor prevents the chord from ever collapsing "
            "to zero however large this is set."))

        self.check_round_te_tip = QCheckBox("Round Tip Corner - Trailing Edge Side")
        self.check_round_te_tip.setChecked(True)
        self.check_round_te_tip.setToolTip(
            "Rounds the corner of the blade planform where the TRAILING edge meets the tip, "
            "independently of the other three corners. The other half of the anti-cavitation lip "
            "described in the leading-edge/tip tooltip above - see that tooltip for the full "
            "engineering context.")
        ig.addWidget(self.check_round_te_tip, 5, 0, 1, 2)
        self.spin_te_tip_round_extent = make_double_spin(0.15, 0.0, 1.00, decimals=2, step=0.01)
        add_row(ig, 6, 'TE/Tip Rounding Extent (% of ref. chord):', self.spin_te_tip_round_extent, tooltip=(
            "Fillet radius for this corner, as a fraction of the blade's reference chord length "
            "(not of the spanwise extent). A safety floor prevents the chord from ever collapsing "
            "to zero however large this is set."))

        self.check_round_le_root = QCheckBox("Round Root Corner - Leading Edge Side")
        self.check_round_le_root.setChecked(True)
        self.check_round_le_root.setToolTip(
            "Rounds the corner of the blade planform where the LEADING edge meets the root/hub, "
            "independently of the other three corners. Offered for the same geometric flexibility "
            "as the tip corners above, using the same true 2D corner-fillet construction. "
            "Mechanically plausible (reducing stress concentration at the blade root, analogous "
            "to a standard fillet), but unlike the tip corners, no Kaplan/Deriaz-specific "
            "published source was found to confirm this as a normalized industrial practice - "
            "treat the default extent as illustrative only.")
        ig.addWidget(self.check_round_le_root, 7, 0, 1, 2)
        self.spin_le_root_round_extent = make_double_spin(0.15, 0.0, 1.00, decimals=2, step=0.01)
        add_row(ig, 8, 'LE/Root Rounding Extent (% of ref. chord):', self.spin_le_root_round_extent, tooltip=(
            "Fillet radius for this corner, as a fraction of the blade's reference chord length "
            "(not of the spanwise extent). A safety floor prevents the chord from ever collapsing "
            "to zero however large this is set."))

        self.check_round_te_root = QCheckBox("Round Root Corner - Trailing Edge Side")
        self.check_round_te_root.setChecked(True)
        self.check_round_te_root.setToolTip(
            "Rounds the corner of the blade planform where the TRAILING edge meets the root/hub, "
            "independently of the other three corners. The other half of the root fillet "
            "described in the leading-edge/root tooltip above - see that tooltip for the full "
            "engineering context.")
        ig.addWidget(self.check_round_te_root, 9, 0, 1, 2)
        self.spin_te_root_round_extent = make_double_spin(0.15, 0.0, 1.00, decimals=2, step=0.01)
        add_row(ig, 10, 'TE/Root Rounding Extent (% of ref. chord):', self.spin_te_root_round_extent, tooltip=(
            "Fillet radius for this corner, as a fraction of the blade's reference chord length "
            "(not of the spanwise extent). A safety floor prevents the chord from ever collapsing "
            "to zero however large this is set."))

        advg.addWidget(self.group_industrial, 3, 0, 1, 2)
        layout.addWidget(group_adv)
        layout.addStretch()

        # Compute button placed OUTSIDE the scroll area (see the same
        # pattern/rationale in _build_tab_general) so it is always fully
        # visible regardless of the content height above it.
        self.btn_compute2 = QPushButton("COMPUTE SOLID 3D BLADE DESIGN")
        self.btn_compute2.setObjectName("ComputeBtn")
        self.btn_compute2.setFont(QFont("Segoe UI", 10, QFont.Bold))
        self.btn_compute2.setToolTip("Runs the full hydraulic + geometric solver and regenerates\nthe 3D blade preview with the current parameters.")
        self.btn_compute2.clicked.connect(self.compute_turbine)
        outer.addWidget(self.btn_compute2)

    # -------------------------------------------------------------------
    # TAB 3: EXPORT MANAGER
    # -------------------------------------------------------------------
    def _build_tab_export(self):
        tab = QWidget()
        self.tabs.addTab(tab, " Export ")
        outer = QVBoxLayout(tab)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        container = QWidget()
        scroll.setWidget(container)
        layout = QVBoxLayout(container)

        self.group_export = QGroupBox(" Export Manager ")
        self.group_export.setObjectName("ExportGroup")
        exg = QGridLayout(self.group_export)
        self.combo_export_type = QComboBox()
        self.combo_export_type.addItems(['Solid Blade (Full)', 'Mean Surface Only'])
        self.combo_export_format = QComboBox()
        self.combo_export_format.addItems(['STL (*.stl)', 'IGES (*.igs / *.iges)'])
        add_row(exg, 0, 'Export Geometry:', self.combo_export_type, tooltip=(
            "Export Geometry Scope:\n"
            "- Solid Blade (Full): watertight extrados/intrados shell with end caps, ready for CAD/3D printing.\n"
            "- Mean Surface Only: single-layer camber surface (thin shell, no thickness)."))
        add_row(exg, 1, 'CAD Format:', self.combo_export_format, tooltip=(
            "CAD File Format:\n"
            "- STL: standard triangulated mesh, widely supported by CAD/CAM/3D-printing software.\n"
            "- IGES: point-cloud fallback writer (not a true B-Rep solid); for a real IGES/STEP solid, "
            "pass the exported STL through a CAD kernel (e.g. FreeCAD/OpenCASCADE)."))

        self.btn_export_cad = QPushButton("Export CAD Geometry")
        self.btn_export_cad.setObjectName("ExportCadBtn")
        self.btn_export_cad.setEnabled(False)
        self.btn_export_cad.setFont(QFont("Segoe UI", 9, QFont.Bold))
        self.btn_export_cad.setToolTip("Exports the currently computed geometry to STL or IGES,\naccording to the options selected above.")
        self.btn_export_cad.clicked.connect(self.export_cad)
        exg.addWidget(self.btn_export_cad, 2, 0, 1, 2)

        self.btn_export_excel = QPushButton("Export Blade Angles (Excel)")
        self.btn_export_excel.setObjectName("ExportExcelBtn")
        self.btn_export_excel.setEnabled(False)
        self.btn_export_excel.setFont(QFont("Segoe UI", 9, QFont.Bold))
        self.btn_export_excel.setToolTip("Exports the inlet/outlet blade angle table (Beta1, Beta2)\nper radial station to an .xlsx workbook.")
        self.btn_export_excel.clicked.connect(self.export_excel)
        exg.addWidget(self.btn_export_excel, 3, 0, 1, 2)

        layout.addWidget(self.group_export)

        note = QLabel(
            "Note: geometry must be computed first (General or Hydrofoil / mesh tab) "
            "before these export options become available.")
        note.setWordWrap(True)
        note.setStyleSheet("color: #8A8F99; font-style: italic;")
        layout.addWidget(note)

        layout.addStretch()

    # -------------------------------------------------------------------
    # UI CALLBACKS
    # -------------------------------------------------------------------
    def on_turbine_type_change(self):
        is_kaplan = (self.combo_type.currentText() == 'Kaplan (Axial)')
        self.group_kaplan.setVisible(is_kaplan)
        self.group_deriaz.setVisible(not is_kaplan)
        if is_kaplan:
            self.spin_rpm.setValue(450); self.spin_q0.setValue(12); self.spin_hn.setValue(15)
            self.spin_eta_h.setValue(0.90); self.spin_eta_v.setValue(0.96); self.spin_eta_o.setValue(0.98)
            self.spin_sigma.setValue(1.25)
            self.spin_thickness.setValue(0.085)
            self.combo_naca.setCurrentText('NACA 00XX (Standard Symmetric)')
        else:
            self.spin_rpm.setValue(150); self.spin_q0.setValue(110); self.spin_hn.setValue(60)
            self.spin_eta_h.setValue(0.91); self.spin_eta_v.setValue(0.97); self.spin_eta_o.setValue(0.98)
            self.spin_sigma.setValue(1.40)
            self.spin_thickness.setValue(0.105)
            self.combo_naca.setCurrentText('Reversible Hydrofoil (Pump-Turbine)')

        # The leading-edge sweep is proven geometrically invalid on a Deriaz
        # backbone (it pulls streamlines off their design sphere - see the
        # checkbox's own tooltip and apply_leading_edge_sweep()'s docstring
        # for the full derivation), so it is disabled and force-unchecked
        # outside Kaplan.
        if hasattr(self, 'check_le_sweep'):
            if not is_kaplan:
                self.check_le_sweep.setChecked(False)
            self.check_le_sweep.setEnabled(is_kaplan)
            self.spin_sweep_angle.setEnabled(is_kaplan)

        # Force the scroll area to recompute its scrollable content height
        # right away (Kaplan has 3 geometry rows, Deriaz has 4), so the
        # Compute button and everything below it stays reachable via the
        # scrollbar instead of getting clipped until the next manual resize.
        # Explicitly re-activating the layout (rather than relying on
        # adjustSize(), which fights with QScrollArea's own size management
        # when setWidgetResizable(True) is set) is the reliable way to make
        # this take effect immediately.
        if hasattr(self, '_general_container'):
            container_layout = self._general_container.layout()
            if container_layout is not None:
                container_layout.activate()
            self._general_container.updateGeometry()

    def on_profile_type_change(self):
        self.group_custom.setVisible(self.combo_naca.currentText() == 'Customized (4-Digit NACA Series)')

    def on_industrial_postproc_toggled(self, checked):
        self.group_industrial.setVisible(checked)

    # -------------------------------------------------------------------
    # MAIN COMPUTATION
    # -------------------------------------------------------------------
    def compute_turbine(self):
        """Top-level orchestration for the COMPUTE button (both tabs). Reads
        every input widget, calls the appropriate simplified backbone
        function (compute_kaplan_backbone or compute_deriaz_backbone),
        optionally applies the industrial post-processing chain (corner
        rounding -> sweep -> alignment -> wake trim) to the resulting mean
        surface, generates the solid via build_solid_blade(), stores the
        results on self for later export/rendering, fills the results
        table (including the Euler-equation breakdown and the preliminary
        distributor estimates), and triggers the 3D render. Any physical
        infeasibility detected by the backbone functions (or by the
        industrial post-processing functions) is surfaced as a
        QMessageBox rather than silently producing invalid geometry."""
        self._set_status("Computing turbine geometry...", kind='busy')
        QApplication.processEvents()

        RPM = self.spin_rpm.value(); Q0 = self.spin_q0.value(); Hn = self.spin_hn.value()
        g = self.spin_g.value()
        eta_h = self.spin_eta_h.value(); eta_v = self.spin_eta_v.value(); eta_o = self.spin_eta_o.value()
        sigma_target = self.spin_sigma.value()
        interp_type = self.combo_interp.currentText()
        rot_sign = 1.0 if 'Counter-Clockwise' in self.combo_rot.currentText() else -1.0

        N_radios = self.spin_nradios.value()
        N_cuerda = self.spin_ncuerda.value()

        is_kaplan = (self.combo_type.currentText() == 'Kaplan (Axial)')

        if is_kaplan:
            res = compute_kaplan_backbone(
                RPM, Q0, Hn, g, eta_h, eta_v, eta_o, sigma_target, interp_type, rot_sign,
                N_radios, N_cuerda, self.spin_rhub.value(), self.spin_rtip.value(), self.spin_lz.value())
        else:
            res = compute_deriaz_backbone(
                RPM, Q0, Hn, g, eta_h, eta_v, eta_o, sigma_target, interp_type, rot_sign,
                N_radios, N_cuerda, self.spin_re_int.value(), self.spin_re_ext.value(),
                self.spin_gamma1.value(), self.spin_gamma2.value())

        if res['error'] is not None:
            title = 'Physical Limit Error' if 'instability' in res['error'] else 'Geometric Error'
            self._set_status(f"Error: {res['error']}", kind='error')
            QMessageBox.critical(self, title, res['error'])
            return

        # --- Hydrofoil profile and solid blade generation ---
        profile_type = self.combo_naca.currentText()
        t_rel_input = self.spin_thickness.value()
        hub_to_tip_ratio = self.spin_hub_to_tip.value()
        m_hub, p_hub = self.spin_m_hub.value(), self.spin_p_hub.value()
        m_tip, p_tip = self.spin_m_tip.value(), self.spin_p_tip.value()

        X_mid, Y_mid, Z_mid = res['X_mid'], res['Y_mid'], res['Z_mid']
        radius_vec = res['radius_vec']
        ref_chord = res['L_chord_reference']
        # Pristine copy, captured before ANY post-processing operation
        # touches the backbone - used exclusively as the REFERENCE geometry
        # for align_trailing_edge_radially()'s theta_hub computation.
        # This is what fixes a real bug: if theta_hub were computed from
        # an already-shortened backbone (e.g. after corner rounding has
        # trimmed the hub streamline), it would silently correspond to a
        # point before the hub's true design endpoint rather than at it.
        X_mid_pristine, Y_mid_pristine, Z_mid_pristine = X_mid.copy(), Y_mid.copy(), Z_mid.copy()

        # --- Industrial post-processing adjustments (optional, off by default) ---
        # Applied to the mean-surface backbone, in THIS ORDER, which
        # matters: leading-edge sweep (Kaplan only, a pure translation
        # that changes nothing else) -> trailing-edge plane cut / wake
        # trim -> corner rounding LAST.
        #
        # Rounding must run AFTER the plane cut, never before - this was a
        # real bug in an earlier pipeline ordering, confirmed from the GUI
        # (rounded corners were invisible whenever the plane cut was also
        # enabled). Two things go wrong if rounding runs first: (1) the
        # plane cut re-derives its own target from the PRISTINE backbone
        # regardless of ordering (so that part is safe either way), but
        # (2) a streamline whose tip/root corner was already rounded may
        # still have far more chord left than the plane cut ultimately
        # keeps, in which case the plane cut trims straight through the
        # already-rounded region and removes it entirely, silently
        # erasing the rounding. Running the plane cut first and rounding
        # the resulting, already-cut edges afterward avoids this: rounding
        # then always acts on the actual final edge, never on a region
        # that a later step might delete.
        if self.check_industrial_postproc.isChecked():
            if is_kaplan and self.check_le_sweep.isChecked():
                X_mid, Y_mid, Z_mid = apply_leading_edge_sweep(
                    X_mid, Y_mid, Z_mid, radius_vec, self.spin_sweep_angle.value())

            # Available for both turbine types; see the "no free lunch"
            # framing in the checkbox tooltip and in
            # align_trailing_edge_radially()'s own docstring - neither
            # leaving this off nor turning it on is simply "the correct"
            # choice. The PRISTINE backbone (captured before any trimming
            # above) is passed as the reference for computing theta_hub -
            # see the comment where it was captured for why.
            if self.check_align_trailing_edge.isChecked():
                X_mid, Y_mid, Z_mid, _target_theta = align_trailing_edge_radially(
                    X_mid, Y_mid, Z_mid, X_mid_pristine, Y_mid_pristine, Z_mid_pristine)

            if (self.check_round_le_tip.isChecked() or self.check_round_te_tip.isChecked()
                    or self.check_round_le_root.isChecked() or self.check_round_te_root.isChecked()):
                X_mid, Y_mid, Z_mid = apply_corner_rounding(
                    X_mid, Y_mid, Z_mid, radius_vec, ref_chord,
                    self.check_round_le_tip.isChecked(), self.spin_le_tip_round_extent.value(),
                    self.check_round_te_tip.isChecked(), self.spin_te_tip_round_extent.value(),
                    self.check_round_le_root.isChecked(), self.spin_le_root_round_extent.value(),
                    self.check_round_te_root.isChecked(), self.spin_te_root_round_extent.value())

            res = dict(res)
            res['X_mid'], res['Y_mid'], res['Z_mid'] = X_mid, Y_mid, Z_mid

        X_ext, Y_ext, Z_ext, X_int, Y_int, Z_int = build_solid_blade(
            X_mid, Y_mid, Z_mid, profile_type, t_rel_input, N_cuerda, N_radios,
            hub_to_tip_ratio, m_hub, p_hub, m_tip, p_tip)

        # --- Save state ---
        self.X_mid, self.Y_mid, self.Z_mid = X_mid, Y_mid, Z_mid
        self.X_sol_ext, self.Y_sol_ext, self.Z_sol_ext = X_ext, Y_ext, Z_ext
        self.X_sol_int, self.Y_sol_int, self.Z_sol_int = X_int, Y_int, Z_int
        self.radius_vec = res['radius_vec']
        self.beta1_deg, self.beta2_deg = res['b1_vec'], res['b2_vec']
        self.Z_optimo_current = res['Z_optimo']
        self.R_color = res['R_color']
        self.last_is_kaplan = is_kaplan
        self.last_backbone = res
        self.is_computed = True
        self.btn_export_cad.setEnabled(True)
        self.btn_export_excel.setEnabled(True)

        # -----------------------------------------------------------------
        # DIRECT HYDRAULIC CALCULATIONS FOR THE RESULTS TABLE
        # -----------------------------------------------------------------
        rho = self.spin_rho.value()  # water density (kg/m3)
        eta_t = eta_h * eta_v * eta_o
        omega = (2.0 * np.pi * RPM) / 60.0

        # Shaft power and torque
        P_hyd = rho * g * Q0 * Hn
        P_mech = P_hyd * eta_t
        Torque = P_mech / omega if omega > 0 else 0.0

        # -----------------------------------------------------------------
        # BUILD RESULTS TABLE
        # -----------------------------------------------------------------
        rows = [
            ('Airfoil Selection', profile_type),
            ('Total Efficiency (eta_t)', f"{eta_t * 100:.2f} %"),
            ('Specific Speed (nq)', f"{res['nq']:.2f}"),
            ('Integrated 3D Chord (L_chord_ref)', f"{res['L_chord_reference']:.3f} m"),
            ('Number of Blades (Z)', f"{res['Z_optimo']}"),
            (res['rm_label'], res['rm_val']),
            ('Gross Hydraulic Power (P_hyd)', f"{P_hyd / 1e3:.2f} kW"),
            ('Net Mechanical Power (P_mech)', f"{P_mech / 1e3:.2f} kW"),
            ('Shaft Torque (T)', f"{Torque / 1e3:.2f} kN.m"),
            ('Angular Speed (omega)', f"{omega:.2f} rad/s"),
        ]

        self.table.setRowCount(len(rows))
        for row_idx, (param_name, param_value) in enumerate(rows):
            item_param = QTableWidgetItem(str(param_name))
            item_val = QTableWidgetItem(str(param_value))

            item_param.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            item_val.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)

            self.table.setItem(row_idx, 0, item_param)
            self.table.setItem(row_idx, 1, item_val)

        # Render in 3D
        try:
            self._render(res, is_kaplan)
        except Exception as e:
            self._set_status(f"3D render error: {e}", kind='error')
            QMessageBox.critical(
                self, "3D Render Error",
                f"The hydraulic/geometric computation succeeded (see results table below), "
                f"but the 3D viewer failed to render:\n\n{e}\n\n"
                f"This is usually a PyVista/VTK version incompatibility, not a problem "
                f"with your design parameters.")
            return

        self._set_status(
            f"Computation complete: {res['Z_optimo']} blades, {profile_type}.", kind='ok')
        self.last_computed_label.setText(
            f"Last computed: {datetime.datetime.now().strftime('%H:%M:%S')}")

    # -------------------------------------------------------------------
    # 3D RENDER (full resolution via PyVista; export uses the same mesh)
    # -------------------------------------------------------------------
    def _render(self, res, is_kaplan):
        Z_optimo = res['Z_optimo']
        delta_angle = (2 * np.pi) / Z_optimo
        metal_color = '#666666'

        # Preserve the user's current camera view (angle + zoom) across
        # recomputes; only fall back to the canonical default view the very
        # first time anything gets rendered, OR whenever the turbine type
        # itself changed since the last render (Kaplan <-> Deriaz have very
        # different geometry scale/position, so reusing the old camera would
        # leave the new geometry off-center or out of view).
        preserve_camera = (
            getattr(self, '_has_rendered_once', False)
            and getattr(self, '_last_rendered_is_kaplan', None) == is_kaplan
        )
        if preserve_camera:
            saved_cam_mid = self.plotter_mid.camera_position
            saved_cam_solid = self.plotter_solid.camera_position

        # -----------------------------------------------------------------
        # 0. COLOR LIMITS AND HUB GEOMETRY
        # -----------------------------------------------------------------
        if is_kaplan:
            R_cubo = res.get('R_hub', res.get('hub_radius', 0.3))
            L_z = res.get('L_z', 0.35)
            hub_mesh = pv.Cylinder(
                center=(0, 0, -L_z / 2.0),
                direction=(0, 0, 1),
                radius=R_cubo,
                height=L_z,
                resolution=60
            )
            # Kaplan color range
            r_min = np.min(res['R_color'])
            r_max = np.max(res['R_color'])
        else:
            Re_int = res.get('Re_int', 2.0)
            Re_ext = res.get('Re_ext', 2.6)
            gamma1 = res.get('gamma1', np.radians(30.0))
            gamma2 = res.get('gamma2', np.radians(60.0))

            # Seamless Deriaz hub mesh
            gg_mesh, tg_mesh = np.meshgrid(
                np.linspace(gamma1, gamma2, 30),
                np.linspace(0, 2 * np.pi, 61)
            )
            Xc = Re_int * np.cos(gg_mesh) * np.cos(tg_mesh)
            Yc = Re_int * np.cos(gg_mesh) * np.sin(tg_mesh)
            Zc = -Re_int * np.sin(gg_mesh)
            hub_mesh = pv.StructuredGrid(Xc, Yc, Zc)

            # Exact color range requested for Deriaz (projected cylindrical radius)
            r_min = Re_int * np.cos(gamma2)
            r_max = Re_ext * np.cos(gamma1)

        # -----------------------------------------------------------------
        # 1. RENDER A: MEAN SURFACE (concentric 2D grid)
        # -----------------------------------------------------------------
        self.plotter_mid.clear()
        self.plotter_mid.set_background(self._plot_bg_color())
        
        self.plotter_mid.enable_depth_peeling()
        
        self.plotter_mid.add_text(f"Mean Surface Design ({Z_optimo} Blades)", position='upper_left', font_size=10, color=self._plot_text_color())

        self.plotter_mid.add_mesh(hub_mesh, color=metal_color, opacity=0.85, smooth_shading=True, ambient=0.45, diffuse=0.65, specular=0.05)

        # Compute continuous projected cylindrical radius matrix
        R_matrix = np.sqrt(self.X_mid**2 + self.Y_mid**2)

        for k in range(Z_optimo):
            ang = k * delta_angle
            Xr = self.X_mid * np.cos(ang) - self.Y_mid * np.sin(ang)
            Yr = self.X_mid * np.sin(ang) + self.Y_mid * np.cos(ang)

            grid_mid = pv.StructuredGrid(Xr, Yr, self.Z_mid)

            # Correct index mapping to avoid streaks or seams
            if grid_mid.dimensions[0] == R_matrix.shape[0]:
                grid_mid['Radius'] = R_matrix.ravel(order='F')
            else:
                grid_mid['Radius'] = R_matrix.T.ravel(order='F')

            self.plotter_mid.add_mesh(
                grid_mid,
                scalars='Radius',
                cmap='viridis',
                clim=[r_min, r_max],
                smooth_shading=True,
                ambient=0.45,
                diffuse=0.65,
                specular=0.05,
                show_scalar_bar=False
            )

        self.plotter_mid.add_scalar_bar(title="Radius (m)", color=self._plot_text_color())
        if preserve_camera:
            self.plotter_mid.camera_position = saved_cam_mid
        else:
            self._set_default_camera_view(self.plotter_mid)
        self.plotter_mid.render()

        # -----------------------------------------------------------------
        # 2. RENDER B: SOLID BLADE (smooth normals to avoid faceted edges)
        # -----------------------------------------------------------------
        self.plotter_solid.clear()
        self.plotter_solid.set_background(self._plot_bg_color())

        naca_text = self.combo_naca.currentText() if hasattr(self, 'combo_naca') else "NACA"
        self.plotter_solid.add_text(f"Solid Runner ({naca_text})", position='upper_left', font_size=10, color=self._plot_text_color())

        self.plotter_solid.add_mesh(hub_mesh, color=metal_color, opacity=0.85, smooth_shading=True, ambient=0.45, diffuse=0.65, specular=0.05)

        # Rebuild the 3D solid with end caps
        F, V = build_export_mesh(
            self.X_sol_ext, self.Y_sol_ext, self.Z_sol_ext,
            self.X_sol_int, self.Y_sol_int, self.Z_sol_int,
            self.X_mid, self.Y_mid, self.Z_mid, is_solid=True
        )

        if len(F) > 0 and len(V) > 0:
            num_faces = F.shape[0]
            padding = np.full((num_faces, 1), 3, dtype=np.int32)
            faces_vtk = np.hstack((padding, F)).astype(np.int32).ravel()

            blade_base = pv.PolyData(V, faces_vtk)

            # SMOOTH NORMALS: removes dark facet lines and edge discontinuities
            blade_base = blade_base.compute_normals(
                cell_normals=False,
                point_normals=True,
                feature_angle=60.0,
                split_vertices=True
            )

            # Color gradient by projected radius
            R_verts = np.sqrt(blade_base.points[:, 0]**2 + blade_base.points[:, 1]**2)
            blade_base['Radius'] = R_verts

            for k in range(Z_optimo):
                angle_deg = np.degrees(k * delta_angle)
                blade_rotated = blade_base.rotate_z(angle_deg, inplace=False)

                self.plotter_solid.add_mesh(
                    blade_rotated,
                    scalars='Radius',
                    cmap='viridis',
                    clim=[r_min, r_max],
                    smooth_shading=True,
                    ambient=0.45,
                    diffuse=0.65,
                    specular=0.05,
                    show_edges=False,
                    show_scalar_bar=False
                )

        # Soft three-point lighting setup
        self.plotter_solid.remove_all_lights()
        self.plotter_solid.add_light(pv.Light(position=(0, 0, 10), focal_point=(0, 0, 0), color=[0.6, 0.6, 0.6]))
        self.plotter_solid.add_light(pv.Light(position=(0, 0, -10), focal_point=(0, 0, 0), color=[0.4, 0.4, 0.4]))
        self.plotter_solid.add_light(pv.Light(position=(5, -5, 2), focal_point=(0, 0, 0), color=[0.5, 0.5, 0.5]))

        self.plotter_solid.add_scalar_bar(title="Radius (m)", color=self._plot_text_color())
        if preserve_camera:
            self.plotter_solid.camera_position = saved_cam_solid
        else:
            self._set_default_camera_view(self.plotter_solid)
        self.plotter_solid.render()

        self._has_rendered_once = True
        self._last_rendered_is_kaplan = is_kaplan

    def _set_default_camera_view(self, plotter):
        """Sets the canonical default camera orientation (isometric, Z-up,
        no roll) and fits the view to the current geometry bounds."""
        plotter.view_isometric()
        plotter.reset_camera()

    def reset_3d_view(self):
        """Restores both 3D viewports to the default camera angle (not just
        re-centering/zooming, but also resetting the viewing orientation)."""
        if hasattr(self, 'plotter_mid'):
            self._set_default_camera_view(self.plotter_mid)
            self._set_default_camera_view(self.plotter_solid)
        self._set_status("3D view reset to default camera position.")

    # -------------------------------------------------------------------
    # CAD EXPORT (STL / IGES fallback)
    # -------------------------------------------------------------------
    def export_cad(self):
        """Exports the currently computed single-blade geometry to STL (via
        write_stl_ascii) or the point-cloud IGES fallback (via
        write_iges_fallback), respecting the 'Solid Blade (Full)' vs
        'Mean Surface Only' choice by passing is_solid through to
        build_export_mesh(). No-op with a warning if nothing has been
        computed yet."""
        if not self.is_computed:
            QMessageBox.warning(self, "Attention", "Compute geometry first.")
            return

        is_solid = (self.combo_export_type.currentText() == 'Solid Blade (Full)')
        is_stl = 'STL' in self.combo_export_format.currentText()

        default_name = 'TurbineBlade_Export.stl' if is_stl else 'TurbineBlade_Export.igs'
        filter_str = 'STL Files (*.stl)' if is_stl else 'IGES Files (*.igs *.iges)'
        file_path, _ = QFileDialog.getSaveFileName(self, "Save CAD Geometry", default_name, filter_str)
        if not file_path:
            return

        try:
            F, V = build_export_mesh(
                self.X_sol_ext, self.Y_sol_ext, self.Z_sol_ext,
                self.X_sol_int, self.Y_sol_int, self.Z_sol_int,
                self.X_mid, self.Y_mid, self.Z_mid, is_solid)

            if is_stl:
                write_stl_ascii(file_path, F, V)
            else:
                write_iges_fallback(file_path, F, V)

            QMessageBox.information(self, "Export Success",
                                     f"\u2713 Geometry exported successfully to:\n{file_path}")
            self._set_status(f"CAD geometry exported to: {file_path}", kind='ok')
        except Exception as e:
            QMessageBox.critical(self, "Export Failed", f"\u2717 Export Error: {e}")
            self._set_status(f"Export failed: {e}", kind='error')

    # -------------------------------------------------------------------
    # EXCEL EXPORT (blade angles)
    # -------------------------------------------------------------------
    def export_excel(self):
        """Exports the per-radius inlet/outlet blade angle table
        (beta1, beta2) computed for the current design to a two-column
        .xlsx workbook via openpyxl - useful for feeding the angle
        distribution into external CAD/CFD tools that expect tabular
        input rather than the STL mesh."""
        if not self.is_computed:
            QMessageBox.warning(self, "Attention", "Compute geometry first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Excel Report", "Blade_Angles_Data.xlsx", "Excel Files (*.xlsx)")
        if not file_path:
            return

        try:
            wb = openpyxl.Workbook()
            ws = wb.active
            ws.title = "Blade Angles"
            ws.append(['Radius_m', 'InletAngle_Beta1_deg', 'OutletAngle_Beta2_deg'])
            for r_val, b1, b2 in zip(self.radius_vec, self.beta1_deg, self.beta2_deg):
                ws.append([float(r_val), float(b1), float(b2)])
            wb.save(file_path)

            QMessageBox.information(self, "Export Success", f"\u2713 Excel report saved to:\n{file_path}")
            self._set_status(f"Excel report saved to: {file_path}", kind='ok')
        except Exception as e:
            QMessageBox.critical(self, "Export Failed", f"\u2717 Excel Export Error: {e}")
            self._set_status(f"Excel export failed: {e}", kind='error')

    # =====================================================================
    # PROFESSIONAL EDITION - ADDITIONAL UTILITIES
    # (Save/Load Configuration, Reset, Copy Table, About)
    # =====================================================================

    def gather_configuration(self):
        """Collects all design parameters into a JSON-serializable dict."""
        return {
            'TurbineType': self.combo_type.currentText(),
            'RPM': self.spin_rpm.value(),
            'RotationDirection': self.combo_rot.currentText(),
            'Q0': self.spin_q0.value(),
            'Hn': self.spin_hn.value(),
            'g': self.spin_g.value(),
            'EtaH': self.spin_eta_h.value(),
            'EtaV': self.spin_eta_v.value(),
            'EtaO': self.spin_eta_o.value(),
            'Sigma': self.spin_sigma.value(),
            'InterpolationScheme': self.combo_interp.currentText(),
            'R_hub': self.spin_rhub.value(),
            'R_tip': self.spin_rtip.value(),
            'L_z': self.spin_lz.value(),
            'Re_int': self.spin_re_int.value(),
            'Re_ext': self.spin_re_ext.value(),
            'gamma1_deg': self.spin_gamma1.value(),
            'gamma2_deg': self.spin_gamma2.value(),
            'HydrofoilProfile': self.combo_naca.currentText(),
            'MaxRelThickness': self.spin_thickness.value(),
            'HubToTipRatio': self.spin_hub_to_tip.value(),
            'm_hub': self.spin_m_hub.value(),
            'p_hub': self.spin_p_hub.value(),
            'm_tip': self.spin_m_tip.value(),
            'p_tip': self.spin_p_tip.value(),
            'N_radii': self.spin_nradios.value(),
            'N_chord': self.spin_ncuerda.value(),
            'ExportGeometry': self.combo_export_type.currentText(),
            'CADFormat': self.combo_export_format.currentText(),
        }

    def save_configuration(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Design Configuration", "TurbineDesign_Config.json", "JSON Files (*.json)")
        if not file_path:
            return
        try:
            cfg = self.gather_configuration()
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(cfg, f, indent=2)
            self._set_status(f"Configuration saved to: {file_path}", kind='ok')
        except Exception as e:
            QMessageBox.critical(self, "Save Failed", f"\u2717 Could not save configuration: {e}")
            self._set_status(f"Save failed: {e}", kind='error')

    def load_configuration(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Design Configuration", "", "JSON Files (*.json)")
        if not file_path:
            return
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                cfg = json.load(f)

            self.combo_type.setCurrentText(cfg['TurbineType'])
            self.spin_rpm.setValue(cfg['RPM'])
            self.combo_rot.setCurrentText(cfg['RotationDirection'])
            self.spin_q0.setValue(cfg['Q0'])
            self.spin_hn.setValue(cfg['Hn'])
            self.spin_g.setValue(cfg['g'])
            self.spin_rho.setValue(cfg['rho'])
            self.spin_eta_h.setValue(cfg['EtaH'])
            self.spin_eta_v.setValue(cfg['EtaV'])
            self.spin_eta_o.setValue(cfg['EtaO'])
            self.spin_sigma.setValue(cfg['Sigma'])
            self.combo_interp.setCurrentText(cfg['InterpolationScheme'])
            self.spin_rhub.setValue(cfg['R_hub'])
            self.spin_rtip.setValue(cfg['R_tip'])
            self.spin_lz.setValue(cfg['L_z'])
            self.spin_re_int.setValue(cfg['Re_int'])
            self.spin_re_ext.setValue(cfg['Re_ext'])
            self.spin_gamma1.setValue(cfg['gamma1_deg'])
            self.spin_gamma2.setValue(cfg['gamma2_deg'])
            self.combo_naca.setCurrentText(cfg['HydrofoilProfile'])
            self.spin_thickness.setValue(cfg['MaxRelThickness'])
            self.spin_hub_to_tip.setValue(cfg['HubToTipRatio'])
            self.spin_m_hub.setValue(cfg['m_hub'])
            self.spin_p_hub.setValue(cfg['p_hub'])
            self.spin_m_tip.setValue(cfg['m_tip'])
            self.spin_p_tip.setValue(cfg['p_tip'])
            self.spin_nradios.setValue(cfg['N_radii'])
            self.spin_ncuerda.setValue(cfg['N_chord'])
            self.combo_export_type.setCurrentText(cfg['ExportGeometry'])
            self.combo_export_format.setCurrentText(cfg['CADFormat'])

            # Refresh panel visibility (Kaplan/Deriaz, Customized) without
            # overwriting the values that were just loaded with factory defaults.
            is_kaplan = (self.combo_type.currentText() == 'Kaplan (Axial)')
            self.group_kaplan.setVisible(is_kaplan)
            self.group_deriaz.setVisible(not is_kaplan)
            self.on_profile_type_change()

            self._set_status(
                f"Configuration loaded from: {file_path} - press COMPUTE to regenerate.", kind='info')
        except Exception as e:
            QMessageBox.critical(self, "Load Failed", f"\u2717 Could not load configuration: {e}")
            self._set_status(f"Load failed: {e}", kind='error')

    def reset_to_defaults(self):
        """Restores every field, and the app state, to the pristine
        Kaplan-default configuration (matching a fresh app start)."""
        # Force the turbine type back to Kaplan first (the canonical
        # starting default) - previously this was never touched, so
        # clicking Reset while on Deriaz only refreshed Deriaz-specific
        # fields and looked like it "did nothing".
        self.combo_type.setCurrentText('Kaplan (Axial)')
        self.combo_rot.setCurrentText('Counter-Clockwise (Standard)')
        self.spin_g.setValue(9.81)
        self.spin_rho.setValue(1000)
        self.combo_interp.setCurrentText('Cubic (Standard)')

        # Explicit call, not just relying on the signal from setCurrentText
        # above (which is a no-op if it was already on Kaplan) - resets
        # RPM/Q0/Hn/eta/sigma/thickness/NACA for Kaplan.
        self.on_turbine_type_change()

        self.spin_rhub.setValue(0.30); self.spin_rtip.setValue(0.65); self.spin_lz.setValue(0.35)
        self.spin_re_int.setValue(2.00); self.spin_re_ext.setValue(2.60)
        self.spin_gamma1.setValue(30.0); self.spin_gamma2.setValue(60.0)

        self.spin_hub_to_tip.setValue(0.65)
        self.spin_m_hub.setValue(2.0); self.spin_p_hub.setValue(4.0)
        self.spin_m_tip.setValue(0.0); self.spin_p_tip.setValue(4.0)
        self.spin_nradios.setValue(100); self.spin_ncuerda.setValue(200)

        self.combo_export_type.setCurrentText('Solid Blade (Full)')
        self.combo_export_format.setCurrentText('STL (*.stl)')

        self.on_profile_type_change()

        # Clear any previous computation results so the reset is visually
        # unambiguous, instead of silently leaving a stale 3D render/table
        # from before that no longer matches the (now reset) input fields.
        self.is_computed = False
        self.table.setRowCount(0)
        self.btn_export_cad.setEnabled(False)
        self.btn_export_excel.setEnabled(False)
        self._has_rendered_once = False
        self._last_rendered_is_kaplan = None
        try:
            self.plotter_mid.clear()
            self.plotter_solid.clear()
            bg_col = self._plot_bg_color()
            self.plotter_mid.set_background(bg_col)
            self.plotter_solid.set_background(bg_col)
            self.plotter_mid.render()
            self.plotter_solid.render()
        except Exception:
            pass

        self.reset_3d_view()
        self._set_status("All parameters reset to defaults. Press COMPUTE to regenerate the geometry.")
        self.last_computed_label.setText("")

    def copy_results_to_clipboard(self):
        if not self.is_computed or self.table.rowCount() == 0:
            QMessageBox.warning(self, "Attention", "Compute geometry first.")
            return
        lines = []
        for row in range(self.table.rowCount()):
            param = self.table.item(row, 0).text()
            value = self.table.item(row, 1).text()
            lines.append(f"{param}\t{value}")
        QApplication.clipboard().setText("\n".join(lines))
        self._set_status("Results table copied to clipboard.", kind='ok')

    def show_about(self):
        QMessageBox.information(self, "About KaplanDeriaz3D", (
            "<b>KaplanDeriaz3D - Hydraulic Turbine Blade Designer</b><br>"
            "Professional Edition (Python / PySide6 / PyVista)<br><br>"
            "Solid 3D blade designer for Kaplan (axial) and Deriaz (diagonal) "
            "hydraulic turbine runners, with parametric hydrofoil profiles, "
            "CAD (STL/IGES) export and Excel blade-angle reports.<br><br>"
            "The underlying hydraulic and geometric formulation (velocity "
            "triangles and beta-angle interpolation laws) is unchanged from "
            "the original computational core.<br><br>"
            "All information can be found on https://github.com/KaplanDeriaz3D."))

# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == '__main__':
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = KaplanDeriaz3DApp()
    window.show()
    sys.exit(app.exec())