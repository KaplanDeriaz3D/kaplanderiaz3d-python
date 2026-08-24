# Contributor's Guide to Improving KaplanDeriaz3D

This document is a practical roadmap for anyone who wants to contribute to
**KaplanDeriaz3D**, an open-source parametric design tool for the 3D blade
geometry of Kaplan (axial-flow) and Deriaz (diagonal-flow) hydraulic
turbines ([github.com/KaplanDeriaz3D](https://github.com/KaplanDeriaz3D)).

The project currently ships three repositories:

- **`kaplanderiaz3d-matlab`** — the original MATLAB App Designer
  application, also distributed as a MATLAB Add-On (`.mlappinstall`).
- **`kaplanderiaz3d-python`** — a full Python port (PySide6 + PyVista),
  replicating the same functionality for users without a MATLAB licence.
- **`kaplanderiaz3d-theory`** — the calculation notebooks (MATLAB Live
  Scripts) and theoretical documentation bridging the underlying master's
  thesis and the published code.

All three are distributed under **GPLv3**: you are free to use, study,
modify, and redistribute the code, on the condition that any derivative
work remains open-source under the same licence.

Every improvement listed below comes directly from the "Future Work"
section of the thesis this project is built on (Section 6.4). Nothing here
is a vague wishlist — each item was identified as a genuine, scoped gap
during the original development, deliberately left out at that time only
for reasons of time and scope, not because it was considered unimportant
or already solved elsewhere. **None of the items below are implemented in
the currently published code** — the published version deliberately keeps
to a validated, minimal hydrodynamic core, and everything described here
is a real opportunity to extend it.

---

## How to use this guide

The six sections below mirror the six subsections of the thesis's future
work chapter, in the same order. Each item includes:

- **What it means** — explained for someone who has not read the thesis.
- **Why it matters** — the concrete gap or limitation it closes.
- **Where to start in the code** — the relevant function(s)/file(s).
- **Difficulty & prerequisites** — a rough, honest estimate.

**Be honest with yourself about prerequisites before picking an item.**
This is a hydraulic-engineering tool, and almost everything below requires
real domain knowledge — hydraulics/turbomachinery theory, CFD, structural
mechanics, or control systems — to contribute meaningfully, not just
general programming ability. There is no item here where "I can code but
don't know anything about turbines" is enough on its own to make useful
progress unsupervised. If that describes you, the most realistic way to
contribute is to pair with someone who has the domain background, or to
pick one of the few sub-items below that are genuinely closer to a pure
software/numerical-methods task once the physics has already been decided
by someone else:

- The **iterative/coupled convergence scheme** needed in Section 6.4.2 (the
  Euler-height, blockage-factor, and blade-count corrections depend on
  each other) is primarily a numerical-methods and software-engineering
  problem *once the underlying equations are settled* — turbomachinery
  background helps you sanity-check the result, but the core difficulty is
  designing and implementing a robust fixed-point/iterative solver, not
  deriving new physics.
- The **automatic hub/shroud radius estimation** in Section 6.4.6 applies
  an existing, closed-form relation between hydraulic inputs and geometry
  — it needs enough hydraulics to understand the formula you're given, but
  not open-ended domain research.
- The **3D rendering of the distributor** in Section 6.4.6 is largely a
  geometry/visualization task once someone else has specified what shape
  needs to be drawn.

Every other item — including, importantly, **Francis turbine support**,
which is one of the most demanding items in this entire guide, not one of
the easiest — genuinely requires the corresponding domain expertise
(hydraulics, CFD, structures, or controls) to contribute to in a
meaningful way. Section difficulty ratings below reflect this honestly.

---

## 6.4.1 — Numerical and Experimental Validation

The tool currently outputs geometry derived analytically from the Euler
turbomachinery equation. **None of this output has been validated against
CFD or experimental data.** This is, by a wide margin, the single most
valuable category of contribution available: everything else in this
guide builds *on top of* a design methodology whose real-world accuracy is
still unverified.

### What's needed

- **CFD simulation of the Deriaz geometry.** Run the runner geometry
  produced by the tool through a Navier–Stokes solver (OpenFOAM [21] is a
  free, open-source option; the [Leguizamón & Avellan Francis inverse
  design paper](https://doi.org/10.3390/en13082020) [23] is a good
  reference for how to structure this kind of validation study end to
  end) across a matrix of net head, flow rate, and rotational speed
  combinations, including computing the flow Reynolds number for each
  case, to characterize the efficiency the geometry actually achieves.
- **Optimal hub-to-shroud radius ratio study.** Use CFD to explore how
  efficiency varies with the ratio between `Re_int` and `Re_ext`
  (Deriaz), and with the range of inlet/outlet cone angles `gamma1`,
  `gamma2` — currently these are free user inputs with no optimality
  guidance behind them at all.
- **Alternative and validated hydrofoil profiles.** The tool currently
  offers five profile families (see `generate_hydro_profile()` in the
  source), all NACA-derived [28]. Contributions here fall into two
  categories: (a) implementing genuinely different profile families
  (e.g. profiles specifically developed for hydraulic — not aeronautical
  — applications), and (b) running CFD on the *existing* profiles to
  check whether their assumed behaviour (e.g. the "Anti-Cavitation" or
  "Low-Torque S-Camber" families) actually holds up.
- **Beta-angle interpolation laws that vary with radius.** The tool
  currently applies one interpolation law (cubic by default; see
  `evaluate_interpolation()`) uniformly across every streamline. A
  genuinely useful extension is a law that varies *by radius*, aimed at
  distributing blade loading more evenly and reducing the chord length
  at larger radii (a known practical concern — long, heavy blade tips are
  structurally and hydraulically undesirable). This would require
  modifying the interpolation call sites inside `compute_kaplan_backbone()`
  and `compute_deriaz_backbone()` so the interpolation weight can depend
  on the streamline index, not just the chordwise/gamma-wise position.
- **Blade thickness and maximum admissible width study.** Currently the
  relative thickness (`t/c`) is a free user input with a linear hub-to-tip
  taper (see `generate_hydro_profile()`). There is no study of what
  thickness is structurally necessary vs. hydraulically excessive.
- **Sensitivity to the model's simplifying assumptions** — most
  importantly, the complete absence of viscosity. The whole backbone
  derivation (`compute_kaplan_backbone()`, `compute_deriaz_backbone()`)
  assumes inviscid, loss-free flow. Quantifying how far this departs from
  real turbulent flow, and what correction factors would bring the model
  closer to reality, is an open, valuable research question.
- **Cavitation studies: NPSH available vs. required.** Net Positive
  Suction Head is not computed anywhere in the current tool. Implementing
  it requires both a plant-level input (installation elevation, suction
  head) and a runner-level output (local pressure distribution along the
  blade, which in turn requires at least a first-order pressure model —
  see Section 6.4.2 below for the Euler-height-based pressure work this
  would naturally build on).

### Where to start in the code

- `compute_kaplan_backbone()` / `compute_deriaz_backbone()` (the physics
  core) for anything touching the velocity triangle or interpolation law.
- `generate_hydro_profile()` for anything about hydrofoil families or
  thickness.
- None of this requires touching the GUI (`KaplanDeriaz3DApp` class) at
  all — it can be developed and validated as pure, standalone Python/CFD
  work against the exported geometry (STL export already exists via
  `write_stl_ascii()`), and only wired into the GUI afterward if desired.

### Difficulty

High for the CFD/validation work itself (requires a working CFD
toolchain, meshing knowledge, and turbomachinery CFD experience), but the
*deliverable* can be a written validation report or a set of correction
factors rather than new code — this is one of the few areas here where a
non-programming contribution is extremely valuable on its own.

---

## 6.4.2 — Refining the Hydrodynamic Model

The current model assumes an **infinite number of infinitely thin
blades** — the classic textbook idealization. This section is about
relaxing that idealization toward the finite-blade reality, and it is
where the project's most conceptually subtle open problem lives (see the
scope-vs-subtlety distinction under "Difficulty" below).

### What's needed

- **Finite-blade correction to the theoretical Euler head** (H<sub>z∞</sub>
  → H<sub>z</sub>). The infinite-blade assumption overestimates the head
  the runner actually extracts. A correction factor (several formulations
  exist in the classical turbomachinery literature, e.g. Stodola's slip
  factor) needs to be incorporated.
- **Contribution of the tangential velocity component (V<sub>θ</sub>) to
  the meridional velocity.** Currently the meridional/axial velocity is
  computed purely from continuity (`V_z = Q_real / Area_paso` in
  `compute_kaplan_backbone()`, and the analogous `Vm1`/`Vm2` terms in
  `compute_deriaz_backbone()`), ignoring any coupling with the swirl
  component. A more complete model would account for this coupling.
- **Blade blockage factor (τ)** in the meridional velocity calculation, as
  proposed by Morabito et al. [8] — the effective flow area is reduced by
  the blades' own thickness, which the current model does not account for
  at all (the backbone is computed for a zero-thickness blade, and
  thickness is only added afterward, purely geometrically, in
  `build_solid_blade()` — it never feeds back into the hydrodynamics).

### The circular-reference problem (the hard part)

**This is explicitly flagged in the thesis as the reason none of the
above three corrections are implemented yet**, and it is worth
understanding precisely, because it's the real technical obstacle a
contributor needs to solve, not just an implementation detail:

- The finite-blade (Euler) correction depends on the **number of blades**.
- The blade **blockage factor τ** also depends on the number of blades
  (more/thicker blades block more of the flow area).
- But the number of blades itself is derived from the **chord length**
  (via the target solidity — see `Z_optimo` in both backbone functions),
  and the chord length depends on the **backbone geometry**, which
  depends on the **beta angles**, which depend on the **corrected
  meridional velocity and corrected head** — which are exactly the
  quantities the two corrections above are trying to compute in the first
  place.

In short: *Z depends on the corrections, and the corrections depend on
Z.* This is a genuine fixed-point problem, not a simple sequential
calculation, and it's why the current code deliberately does *not*
implement any of these three corrections — doing so naively (e.g. just
picking an initial guess for Z and never revisiting it) would silently
produce an inconsistent design.

**What's needed:** an iterative or coupled solution scheme that resolves
this mutual dependency robustly — for example, starting from the current
infinite-blade Z as an initial guess, computing the corrected head,
meridional velocity, and blockage factor, deriving a new Z from that, and
repeating until convergence, with explicit checks for stability and
non-convergence (the loaded blade count could ping-pong indefinitely for
poorly chosen initial parameters if this isn't done carefully).

- **More rigorous, Deriaz-specific efficiency estimation.** The three
  efficiency inputs (`eta_h`, `eta_v`, `eta_o` — hydraulic, volumetric,
  mechanical/organic) are currently free numeric inputs the user sets
  directly (see `self.spin_eta_h` / `_eta_v` / `_eta_o` in
  `_build_tab_general()`), defaulting to generic reference values. A
  genuinely useful contribution would replace this with an estimation
  method specific to the Deriaz topology (published efficiency
  correlations, or a simplified loss-breakdown model), rather than
  requiring the user to already know a sensible value.

### Where to start in the code

`compute_kaplan_backbone()` and `compute_deriaz_backbone()` are the two
functions to modify — specifically the velocity-triangle and blade-count
sections. Because of the circular-dependency problem above, this is
**not** a localized, one-line change: it likely requires restructuring
these functions around an internal convergence loop. Budget real design
time before writing code.

### Difficulty

High. Requires solid turbomachinery theory background and comfort with
iterative/fixed-point numerical schemes. In terms of *scope*, this is more
contained than Francis support (Section 6.4.6) — it works within the
existing backbone functions rather than requiring new infrastructure — but
it is arguably the most conceptually **subtle** single problem in this
guide: the circular dependency is easy to get subtly wrong (e.g. a
convergence loop that looks like it works but silently settles on an
inconsistent blade count), so it demands real care, not just raw effort.
It is also arguably the most valuable single item here: it is the direct
bridge between the current idealized model and a design tool whose
*numbers* (not just its shape) can be trusted.

---

## 6.4.3 — Pump Mode Behaviour and Regulation

The Deriaz turbine's defining advantage is that it can run in **both
turbine and pump mode**, with adjustable blade pitch. None of this
reversible/regulation behaviour is modelled by the tool today — it only
ever designs a single, fixed operating point.

### What's needed

- **Detailed study of Deriaz pump-mode operation**, specifically
  **quantifying the real contribution of centrifugal force to the
  generated manometric head.** In pump mode, part of the head a
  centrifugal/mixed-flow pump develops comes from the centrifugal
  acceleration of the fluid as it moves outward through the impeller —
  the Euler equation's "U₁² − U₂²" term. The current model does not
  isolate or quantify this term at all in a pump-mode context; doing so
  properly requires deriving the pump-mode velocity triangles (which are
  not simply the turbine-mode ones reversed) from scratch.
- **How blade pitch angle regulates flow rate and power**, in both
  turbine and pump mode. The tool currently designs one blade angle
  distribution for one design point; it says nothing about how that
  distribution would need to change, dynamically, to regulate output
  away from the design point — which is precisely what makes the Deriaz
  topology valuable in the first place (see the thesis's Chapter 1
  motivation).
- **Optimal location of the blade pivot axis** (the mechanical hinge each
  blade rotates on to change its pitch), positioned so that the blade's
  own mass is balanced about that axis, minimizing the torque the
  actuation servomechanism has to supply. This is a rigid-body mechanics
  problem layered on top of the existing blade geometry — it needs the
  blade's mass distribution (which in turn needs a material/density
  assumption not currently part of the tool at all) and its centre of
  mass relative to a candidate pivot line.

### Where to start in the code

This is largely **new functionality**, not a modification of existing
functions — there is currently no pump-mode code path, no regulation
model, and no mechanical/mass model anywhere in the codebase. A
reasonable approach: start with a standalone analysis module (separate
from `compute_deriaz_backbone()`) that consumes an already-computed
backbone and blade solid, and adds the pump-mode velocity triangles, the
regulation curves, and the pivot-axis optimization as additional,
composable outputs.

### Difficulty

High for the pump-mode physics (genuinely new derivation, not a
refinement of existing equations); moderate for the pivot-axis
optimization (a self-contained rigid-body mechanics problem once a mass
model exists).

---

## 6.4.4 — Structural Analysis

The tool currently performs **zero structural analysis**. Blade thickness
comes entirely from the hydrofoil profile selection (Section 6.4.1 above)
with no stress, fatigue, or deflection check behind it.

### What's needed

- **Structural study of the blades**, in particular of the **wake
  generated at the trailing edge** — real Kaplan and Deriaz blades are
  never left with a sharp, tapering-to-zero trailing edge in practice;
  they are trimmed back to a manufacturable, structurally sound edge. The
  specific open question the thesis identifies: **is removing this wake
  via geometric trimming strictly necessary, or can it be improved upon**
  (a different trim strategy, a different edge treatment entirely,
  etc.)? This connects directly to the corner-rounding / edge-trimming
  discussion that has come up repeatedly during this project's
  development — any contribution here should engage with *why* a naive
  trim can look geometrically inconsistent across the span (streamlines
  at different radii don't sweep the same circumferential angle for the
  same physical trim), not just implement a cosmetic cut.

### Where to start in the code

There is no structural module in the codebase at all today. A sensible
starting point: a standalone stress-analysis script (even a simple beam-
theory or thin-shell approximation) that consumes the exported solid
geometry (`build_export_mesh()` / the STL export) and estimates stress
concentration at the trailing edge under a representative hydraulic
loading (which itself would need at least a first-order pressure
distribution — see the NPSH/pressure discussion in Section 6.4.1). For a
more rigorous treatment, exporting to an external FEA tool (FreeCAD's FEM
workbench, or a Python FEA library) driven by the existing STL/mesh
export is a reasonable integration point.

### Difficulty

Moderate to high, depending on the rigor targeted — a basic stress-
concentration estimate is approachable; a full fatigue/fracture analysis
under real hydraulic and centrifugal loading is a substantial mechanical
engineering undertaking on its own.

---

## 6.4.5 — Control and Automation

Entirely unaddressed today: the tool designs a runner, but says nothing
about how a real power plant would *operate* one.

### What's needed

- **PLC control system** that jointly regulates the distributor
  (wicket-gate) angle and the runner blade angle, based on reservoir
  level, available flow rate, and the power demanded by the electrical
  grid. This is a control-systems contribution (a regulation strategy —
  e.g. a governor curve or a coupled cam-based schedule, matching how
  real double-regulated units like the Deriaz are governed in practice —
  and potentially an actual PLC-targeted implementation, e.g. in
  structured text or ladder logic, or a simulation of one).

### Where to start in the code

This is fully outside the existing codebase's scope (the tool is a
*design* tool, not a *control-system* tool) — realistically a new,
separate component/repository that *consumes* a KaplanDeriaz3D design
(e.g. the exported blade-angle-vs-radius data already produced by
`export_excel()`) as one of its inputs, rather than a modification to the
existing application.

### Difficulty

Moderate, but requires control-systems / industrial-automation background
rather than hydraulics or CAD skills — a good entry point for a
contributor from a different engineering background than the rest of this
project.

---

## 6.4.6 — Extending the Software Itself

This is the **most self-contained, most approachable section** for a new
contributor, and the best place to start if you want a concrete first
pull request. Every item here is squarely a software-engineering task on
top of a physics model that's already implemented and working.

### What's needed

- **Automatic hub/shroud radius calculation.** Currently `R_hub`/`R_tip`
  (Kaplan) and `Re_int`/`Re_ext` (Deriaz) are free numeric inputs the user
  must already know reasonable values for (see the Kaplan/Deriaz Geometry
  groups in `_build_tab_general()`). The requested improvement: derive
  these automatically from net head, flow rate, rotational speed, and a
  user-chosen radius ratio ν (hub-to-shroud), removing the need to guess
  starting values. This is a genuinely well-scoped, self-contained
  addition: a new function (parallel to `compute_kaplan_backbone()`) that
  takes the hydraulic inputs plus ν and returns radius estimates, wired
  into the General tab with an "Estimate Radii" button or similar,
  leaving the fields editable afterward rather than silently overriding
  user choices.
- **Francis turbine support**, extending KaplanDeriaz3D beyond its current
  Kaplan/Deriaz scope. This is the single largest item in this entire
  guide. Francis runners follow a genuinely different design paradigm
  from Kaplan/Deriaz (radial-to-axial mixed flow, not a simple
  cylindrical or spherical streamline family), and the most complete
  published open-source methodology for it is the **inverse design
  method** described by Leguizamón & Avellan [23]
  (`doi.org/10.3390/en13082020`) — their MATLAB implementation is itself
  published open-source (linked from that paper's Supplementary
  Materials) and is the most promising concrete starting point for anyone
  taking this on, since it means the underlying mathematics does not need
  to be re-derived from a paper description alone. Be aware this is a
  substantial undertaking on the scale of a dedicated sub-project, not an
  incremental patch — a solver for two coupled PDEs (mean-flow stream
  function and periodic-flow potential) on a curvilinear mesh, plus a
  meridional-channel geometry generator, none of which exists anywhere in
  the current codebase.
- **Integrated distributor (wicket gate) design**, calculating the
  optimal wicket-gate blade angle as a function of flow rate and design
  head, directly inside the application — this would let the tool be used
  for the real-world control of an existing plant's distributor across
  varying operating flow rates, not just for new-runner design.
- **3D rendering of the distributor**, and calculation of the optimal
  number of guide vanes according to established design criteria (e.g. a
  solidity-based criterion, analogous to how the runner's own blade count
  `Z_optimo` is already derived from a target solidity in
  `compute_kaplan_backbone()` / `compute_deriaz_backbone()` — the same
  logic pattern could plausibly be adapted for the distributor).
- **Selection of the most appropriate hydrofoil profile for the
  distributor's own guide vanes** — a distinct question from the runner
  blade profile selection already implemented, since the guide vanes see
  a different flow regime (no rotation, primarily deflecting the
  approach flow angle).

### Where to start in the code

- **Radius estimation**: a new pure function alongside
  `compute_kaplan_backbone()` / `compute_deriaz_backbone()`, plus a small
  GUI addition in `_build_tab_general()`.
- **Francis support**: effectively a new, parallel backbone module (it
  will NOT fit the existing `compute_kaplan_backbone()` /
  `compute_deriaz_backbone()` pattern, since Francis geometry isn't a
  simple streamline family on a cylinder or sphere) plus new GUI tabs/
  turbine-type handling throughout `KaplanDeriaz3DApp`.
- **Distributor design/rendering**: a new module, most naturally
  integrated as an additional, optional stage after the existing runner
  computation in `compute_turbine()`, reusing the existing 3D-rendering
  infrastructure (`_render()`) as a template for how to add a new mesh to
  the PyVista viewports.

### Difficulty

Wide range, and **not uniformly easy** — this section has both the most
approachable item in the whole guide and the single hardest one. Radius
estimation is a good first pull request (low-to-moderate difficulty,
clearly scoped, applies an existing formula, immediately testable against
the existing GUI). **Francis support, by contrast, is realistically the
largest and most hydraulically/numerically demanding undertaking
described anywhere in this guide** — it requires solving coupled PDEs on
a curvilinear mesh (see the equations discussed above), genuine
turbomachinery design expertise, and should be treated as its own project
phase, not an entry-level contribution. The distributor items sit in
between: real hydraulic-design content (blade angle, vane count), but
self-contained and buildable incrementally, with the pure 3D-rendering
piece being comparatively approachable on its own.

---

## General Notes for Contributors

- **The mathematical core has no Qt/GUI dependency.** Every function
  described above that lives in the "mathematical core" section of the
  source file (`evaluate_interpolation()`, `compute_kaplan_backbone()`,
  `compute_deriaz_backbone()`, `generate_hydro_profile()`,
  `build_solid_blade()`, the mesh/export functions) is pure NumPy with no
  PySide6 or PyVista dependency. You can develop and test hydraulic/
  geometric contributions entirely from a plain Python script or notebook
  before touching the GUI at all — this is the recommended workflow for
  anything in Sections 6.4.1–6.4.4.
- **Physical infeasibility should fail loudly, not silently.** Both
  backbone functions already follow this pattern: if a velocity triangle
  becomes physically impossible (`V_theta1 >= U`), the function returns
  `{'error': <message>}` rather than producing invalid geometry. Any new
  physics you add should follow the same discipline — a wrong-looking
  number is far more dangerous than a clear error message.
- **When in doubt about scope, open an issue first.** Several items above
  (Section 6.4.2 in particular) are non-trivial design decisions, not
  just implementation work — get alignment on the approach before
  investing in an implementation.

---

## Reference List

Numbered to match the citations used above; full details in the thesis
bibliography.

- **[8]** A. Morabito, G. de Oliveira e Silva, P. Hendrick, "Deriaz
  pump-turbine for pumped hydro energy storage and micro applications,"
  *Journal of Energy Storage*, vol. 24, 2019.
  [doi:10.1016/j.est.2019.100788](https://doi.org/10.1016/j.est.2019.100788)
- **[21]** OpenCFD Ltd, *OpenFOAM: The Open Source CFD Toolbox*,
  [openfoam.com](https://www.openfoam.com/)
- **[23]** S. Leguizamón, F. Avellan, "Open-Source Implementation and
  Validation of a 3D Inverse Design Method for Francis Turbine Runners,"
  *Energies*, vol. 13, no. 8, p. 2020, 2020.
  [doi:10.3390/en13082020](https://doi.org/10.3390/en13082020)
- **[25]** S. L. Dixon, *Fluid Mechanics and Thermodynamics of
  Turbomachinery*, 4th ed., Butterworth-Heinemann, 1998.
- **[28]** I. H. Abbott, A. E. von Doenhoff, L. S. Stivers, "Summary of
  Airfoil Data," NACA Report 824, 1945.

---

*This guide is derived from Section 6.4 ("Trabajos futuros") of the
master's thesis underlying KaplanDeriaz3D. For the full mathematical
derivation of the current model, see the thesis document and the
`kaplanderiaz3d-theory` repository.*
