# Industrial Post-Processing Options — What This Variant Adds

**File:** `kaplan_deriaz_python_Post_Processing.py`

## A quick disclaimer before anything else

I am not a turbomachinery engineer, and none of what follows should be read as
validated industrial practice. This file is a variant of the main
KaplanDeriaz3D application with an extra, **optional** post-processing stage
bolted on top of the validated hydrodynamic core. It exists because the
raw, mathematically "ideal" blade geometry the core produces — every
streamline integrated exactly to its own Euler-equation design point — is not
what a real, manufacturable runner looks like. Real runners get trimmed,
rounded, and finished before they're cast or machined, and I wanted to see
what a first, honest attempt at automating that finishing step could look
like.

I don't know whether the specific choices made here (which surfaces to cut,
which corners to round, by how much) match what an actual manufacturer would
do. What I *can* say is that every option below was built with a
**"no free lunch" rule**: nothing is presented as a free improvement, and the
code and tooltips are explicit about what each option costs in exchange for
what it fixes. My hope is that this is useful less as a finished feature and
more as a **worked example of the kind of finishing logic a more rigorous
implementation should eventually include** — a sketch of the shape of the
problem, for someone with the right background to correct, extend, or
replace.

---

## The problem this tries to address

The core hydrodynamic model (shared with the base application, unchanged
here) integrates each streamline independently, from its own leading edge to
its own trailing edge, purely from the local Euler-equation velocity
triangle. Nothing in that derivation forces different streamlines (different
radii) to agree on *where in space* their trailing edges end up. In practice
they don't: the total circumferential angle a streamline sweeps through
(theta) can differ by tens of degrees between the hub and the tip. Visually,
this produces the untrimmed blade's characteristic "tail" — a stretched,
twisted-looking trailing edge that doesn't correspond to anything a real,
finished runner would have. Left alone, it's not just a cosmetic issue: an
irregular discharge edge is a plausible source of stress concentration,
vibration, and wake irregularity downstream.

Everything in this file's post-processing panel is an attempt at one piece
of the finishing work a real design process would apply on top of that raw
geometry.

---

## The options, one by one

### 1. Apply Leading-Edge Sweep (Axial Translation) — Kaplan only

A rigid translation of each streamline along the shaft axis, growing
linearly with distance from the hub, at a rate set by a sweep angle (default
10°, range 0–45°). This is the one option in the whole panel that is
genuinely free: because it moves a whole streamline by the same constant
offset, it changes neither its chord length, nor its local blade angle, nor
its total swept angle — only its position in space. It's disabled for the
Deriaz geometry: a Kaplan streamline lives on a cylinder, which a Z
translation preserves exactly, but a Deriaz streamline lives on a sphere,
which the same translation would pull it off of (verified numerically —
up to several tenths of a metre of drift at the tip for a 10° sweep).

### 2. Align Trailing Edge (Plane Cut) — both turbine types

This is the main attempt at actually closing the "tail" problem described
above. It cuts the blade with a flat plane and keeps only the material
behind it. The plane is defined by two lines that both pass through the
rotation axis: the axis itself, and the line from the origin to the hub's
own, naturally-computed trailing-edge point. Two lines sharing a point
define a plane; because both of these pass through the axis, the resulting
plane sits at a single, constant circumferential angle — exactly the hub's
own discharge angle — and cutting with it is equivalent to trimming every
other streamline back until its own accumulated angle first reaches that
same value.

**The cost, stated plainly:** every streamline other than the hub is cut off
before it reaches the point where its exit tangential velocity was designed
to hit zero. It delivers less specific energy than the design intended, and
discharges with some residual swirl. This is not a rounding error — depending
on the geometry, a large fraction of a streamline's chord near the tip can
be removed. The option is off by default for exactly this reason: it should
be a deliberate choice, not a silent default.

### 3–6. Corner rounding — four independent corners

The blade planform has four corners where an edge meets a span boundary:
leading edge/tip, trailing edge/tip, leading edge/root, trailing edge/root.
Each can be rounded independently, with its own on/off switch and its own
extent (default 0.15, i.e. 15% of the reference chord, all four enabled by
default). The construction is a genuine 2D fillet — the same geometry a CAD
corner fillet uses — scaled to the blade's own chord length rather than its
span, applied to the mean surface before thickness is added so the hydrofoil
profile naturally tapers to nothing at the new edge.

The tip corners are the closest thing here to an established practice: they
approximate the anti-cavitation lip documented on some commercial Kaplan
runners' suction side near the periphery. The root corners are offered for
the same geometric flexibility and are mechanically plausible (reducing
stress concentration at the root, like an ordinary fillet), but I could not
find a Kaplan/Deriaz-specific source confirming them as standard practice —
treat their default extent as illustrative, not as a validated figure.

**A pipeline-order detail that matters:** rounding is applied *after* the
plane cut above, not before. Rounding a corner and then cutting it away
again with a separate operation would silently undo the rounding; cutting
first and rounding the resulting, already-final edge afterward avoids that.

---

## A worked example of a subtlety, for whoever picks this up next

One thing this exercise made concrete for me: the same-looking operation can
behave very differently depending on the underlying geometry. The plane-cut
idea above only trims arc length along a streamline's own existing curve —
it never displaces a point off of it — so it stays exactly on a cylinder
(Kaplan) or a sphere (Deriaz) by construction, verified numerically. A cone
cut based on the streamline's meridional (radius/axial) position, by
contrast, controls something completely different from the plane cut (it
never touches the circumferential angle at all) and turned out to be
essentially degenerate on the Kaplan geometry specifically, because every
Kaplan streamline already shares the same axial trailing-edge position by
construction — there's no meridional "tail" there for a cone to cut back in
the first place. An earlier version of this file included such a cone-cut
option for the Deriaz case; it was removed after further testing suggested
it wasn't behaving reliably enough to keep. I mention it here in case it
saves the next person some time re-deriving why a seemingly natural idea
along these lines needs care.

---

## What I think is still missing

Framed honestly, as a non-expert: this post-processing stage only addresses
the *shape* of the finishing cut, not the *physics* of what that cut does to
the flow. A more complete treatment would presumably want to quantify the
actual efficiency loss and residual swirl each option introduces (rather
than describing it only qualitatively, as the tooltips do here), and to
check the resulting geometry against some structural criterion, not just a
geometric fillet radius. I'd also guess a real design process iterates this
finishing step *together with* the hydrodynamic design rather than bolting
it on afterward, which this file deliberately does not attempt.

If any of this is useful as a starting sketch for someone who actually knows
the subject, that's exactly the spirit it's offered in.
