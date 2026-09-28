# Stereographic Projection Tool (Wulff Net)

A Python/matplotlib tool for plotting crystallographic stereographic
projections (poles and zones/great circles), in the spirit of
[pycotem's StereoProj](https://mompiou.github.io/pycotem/stereoproj/).

## Status

Early prototype. Core projection math and the `StereoProjector` class work
and have been sanity-checked visually (see "Verified so far" below). A
Streamlit UI (`app.py`) now sits on top of it. No crystal-system
generalization, no symmetry expansion yet.

## Files

- `stereoproj.py` — the projection engine. Single class `StereoProjector`.
  Run directly (`python stereoproj.py`) to regenerate `stereoproj_demo.png`,
  a quick visual sanity check.
- `app.py` — Streamlit UI on top of `StereoProjector`. Run with
  `streamlit run app.py`. Lets a student add poles (hkl) and zones [uvw]
  one at a time via forms, see them appear on the projection immediately,
  and "Undo last" / "Clear all" a single combined history (poles and
  zones interleaved in the order they were added, so undo removes
  whichever came last regardless of kind). Each added item gets the next
  color from matplotlib's default cycle so successive entries are easy
  to tell apart; a "Currently plotted" list below the plot echoes each
  item's label and color swatch.
- `requirements.txt` — `numpy`, `matplotlib`, `adjustText`, `streamlit`.

## Dependencies

- **adjustText** (https://github.com/Phlya/adjustText): auto-repels
  overlapping pole/zone labels and draws a thin leader line back to the
  point/mid-arc when a label had to move. Used in `plot()`; pass
  `adjust_labels=False` to fall back to raw, unadjusted label placement.
- **streamlit**: powers `app.py`. Pole/zone labels use matplotlib mathtext
  (`$...$`) internally, which Streamlit also renders as inline LaTeX in
  markdown, so the same label strings are reused verbatim in the UI's
  "Currently plotted" list.

## How the tool is meant to work (agreed design)

1. Always start from an **empty projection circle** (the primitive/equatorial
   boundary circle) — this is the canvas everything else is drawn on.
2. From there the user adds, in any combination:
   - **Zones** — entered as a zone axis `[uvw]`. The tool draws the great
     circle 90° from that axis (the locus of poles of every plane containing
     that direction). This was chosen over "pick two poles" or "click two
     points" as the input method.
   - **Poles** — entered as Miller indices `(hkl)`, plotted as points.

## Design decisions made so far (don't re-litigate without reason)

- **Platform**: plain Python script/matplotlib, not a web artifact or GUI —
  meant to be run/extended like pycotem, not embedded in chat.
- **Pole input**: Miller indices only, for now. They are currently treated
  as **direct Cartesian components** — exact for a cubic lattice, an
  approximation otherwise. The user explicitly deferred general crystal
  parameters + symmetry to a *later* addition (see "Next steps").
- **Zone input**: a zone axis `[uvw]` entered directly, not derived from
  two clicked/selected poles (that may still be worth adding later as a
  convenience, but isn't the primary input path).
- **Projection convention** (matches Stereoproj / DoITPoMS / orix):
  - Right-handed `(x, y, z)` frame, unit vectors.
  - Equal-angle (stereographic) projection. **Upper-hemisphere** poles
    (`z >= 0`) project from the south pole `(0,0,-1)` onto the `z = 0`
    plane and are drawn as **filled** markers:
    `X0 = x / (1 + z)`, `Y0 = y / (1 + z)`.
  - **Lower-hemisphere** poles (`z < 0`) project from the north pole
    `(0,0,+1)` instead, using their own coordinates (no antipode flip),
    and are drawn as **open** markers, larger by default than the filled
    ones: `X0 = x / (1 - z)`, `Y0 = y / (1 - z)`. This way a pole and its
    through-the-plane counterpart `(h,k,l)` / `(h,k,-l)` land at the same
    point, and the open ring visibly encloses the filled dot rather than
    one marker hiding the other.
  - Both are then rotated 90° in-plane, `(X, Y) = (Y0, -X0)`, so that
    `[100]` plots at the south (bottom) of the diagram and `[-100]` at
    the north (top).

## Verified so far

Ran the `__main__` demo (poles `(001)`, `(100)`, `(111)`, `(110)`;
zones `[001]`, `[100]`, `[110]`) and visually confirmed:
- `(001)` projects to the exact center.
- `(100)`, `(110)` land on the boundary circle (both have `z = 0`).
- Zone `[001]` reproduces the boundary circle exactly (axis along `z`).
- Zones `[100]` and `[110]` (axes in the equatorial plane) each project as
  a straight diameter through the center, as expected — any zone axis
  lying in the equatorial plane produces a great circle that itself
  contains the polar axis, which projects as a straight line, not an arc.

Also confirmed after the north/south rotation + open/filled hemisphere
change: `(100)` plots at the bottom (south) of the diagram and `(-100)`
at the top (north); a coincident pair like `(111)`/`(11-1)` renders as
an open ring (lower hemisphere) enclosing a smaller filled dot (upper
hemisphere) at the same point.

## Next steps / open questions (raised but not yet decided)

- **Symmetry-equivalent poles**: pycotem has an "add all symmetric poles"
  option. Worth adding for cubic point-group symmetry even before general
  crystal-parameter support lands.
- **Style sheet**: marker/line kwargs are already exposed per-call
  (`add_pole(..., marker=, color=, size=, open_size=...)`,
  `add_zone(..., color=, linestyle=, linewidth=...)`). Fill is now
  automatic (filled = upper hemisphere, open = lower), not a manual
  `fill=` flag. Undecided whether to add default style presets (e.g.
  planes vs. directions get different default markers).
- **General crystal systems**: planned future addition. Convert Miller
  indices via the direct/reciprocal metric matrices, as in Stereoproj:
  - Direct space: `[x,y,z]^T = D [u,v,w]^T`
  - Reciprocal space (pole normals): `[x,y,z]^T = D* [h,k,l]^T`,
    with `D* = (D^-1)^T`
  - `D` depends on lattice parameters `a,b,c,α,β,γ` (see reference below
    for the explicit matrix form). This also unlocks the metric tensor
    `G* = D*ᵀD*` for d-spacing calculations, if ever wanted.
  - Would also need per-crystal-system point-group symmetry tables to
    support "add all symmetric" for non-cubic systems.

## Non-goals (unless explicitly requested later)

Diffraction patterns, Kikuchi lines, Euler-angle/orientation tracking,
Schmid factor, tilt-holder coordinate systems — all present in pycotem,
none requested for this tool. Don't add speculatively.

## Reference material (equations/conventions this was built from)

- pycotem StereoProj docs — construction, projection formula, metric
  matrices `D`/`D*`, hexagonal 4-index notation:
  https://mompiou.github.io/pycotem/stereoproj/
- DoITPoMS — identifying poles via zone intersections / Weiss zone law:
  https://www.doitpoms.ac.uk/tlplib/stereographic/wulff_identify_poles.php
- orix — Wulff net grid conventions (meridian/parallel spacing), reference
  implementation in a mature crystallography library:
  https://orix.readthedocs.io/en/latest/examples/stereographic_projection/wulff_net.html
