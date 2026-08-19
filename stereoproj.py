"""
stereoproj.py
=============

A minimal equal-angle (Wulff) stereographic projection tool.

Workflow:
    1. Start from an empty projection circle (the equator / primitive circle).
    2. Add zones (great circles), each defined by a zone axis [uvw].
    3. Add poles (hkl), plotted as points.
    4. plot() draws everything together.

Convention (matches Stereoproj / DoITPoMS / orix):
    - Directions live in a right-handed (x, y, z) frame.
    - Projection is from the south pole (0, 0, -1) onto the z = 0 plane.
    - Only the upper hemisphere (z >= 0) is kept; a direction with z < 0
      is replaced by its antipode -[uvw]/-(hkl) so it can still be shown
      (this mirrors what Stereoproj does, and is flagged in the label).
    - A unit direction (x, y, z) with z >= 0 projects to:
          X = x / (1 + z)
          Y = y / (1 + z)

Current limitation: Miller indices [uvw] and (hkl) are treated as direct
Cartesian components. This is EXACT for cubic crystals, and a reasonable
first approximation otherwise. A general triclinic version (using the
direct/reciprocal metric matrices D and D* = (D^-1)^T, as in Stereoproj)
can be added later once crystal parameters + symmetry are wanted.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


class StereoProjector:

    def __init__(self, figsize=6.0):
        self.poles = []
        self.zones = []
        self.figsize = figsize

    # ------------------------------------------------------------------ #
    # geometry helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _normalize(v):
        v = np.array(v, dtype=float)
        n = np.linalg.norm(v)
        if n == 0:
            raise ValueError("Zero vector cannot be normalized.")
        return v / n

    @staticmethod
    def _project(vec):
        """Project a single upper-hemisphere unit vector to (X, Y)."""
        x, y, z = vec
        return x / (1.0 + z), y / (1.0 + z)

    @staticmethod
    def _fmt(n):
        """Format an index with an overbar for negative values, e.g. -1 -> '1̄'."""
        n = int(round(n))
        return f"{abs(n)}\u0305" if n < 0 else f"{n}"

    def _index_label(self, brackets, h, k, l):
        left, right = brackets
        return left + "".join(self._fmt(n) for n in (h, k, l)) + right

    # ------------------------------------------------------------------ #
    # public API
    # ------------------------------------------------------------------ #
    def add_pole(self, h, k, l, label=None, marker='o', color='black',
                 fill=True, size=60):
        """Add a pole (hkl) / direction [uvw] from Miller indices."""
        v = self._normalize([h, k, l])
        flipped = False
        if v[2] < -1e-9:
            v = -v
            h, k, l = -h, -k, -l
            flipped = True
        if label is None:
            label = self._index_label("()", h, k, l)
            if flipped:
                label += "*"   # marks that the antipode was plotted instead
        self.poles.append(dict(vec=v, label=label, marker=marker,
                                color=color, fill=fill, size=size))
        return self

    def add_zone(self, u, v, w, label=None, color='tab:red',
                 linestyle='-', linewidth=1.2, n_samples=361):
        """
        Add the great circle for zone axis [uvw]: the locus of poles of
        every plane that contains this direction (i.e. all points 90
        degrees away from the axis).
        """
        axis = self._normalize([u, v, w])
        if label is None:
            label = self._index_label("[]", u, v, w)
        self.zones.append(dict(axis=axis, label=label, color=color,
                                linestyle=linestyle, linewidth=linewidth,
                                n_samples=n_samples))
        return self

    def clear(self):
        self.poles.clear()
        self.zones.clear()
        return self

    # ------------------------------------------------------------------ #
    # great-circle sampling
    # ------------------------------------------------------------------ #
    def _great_circle_arcs(self, axis, n_samples):
        """
        Sample the great circle perpendicular to `axis`, keep only the
        upper-hemisphere (z >= 0) portion, and return it as a list of
        (N, 2) arrays of projected (X, Y) points -- one array per
        continuous arc (there is normally exactly one arc; it is the
        full boundary circle when axis = +/-z).
        """
        if abs(axis[2]) > 1 - 1e-9:
            # axis along z: the perpendicular great circle IS the equator
            e1 = np.array([1.0, 0.0, 0.0])
            e2 = np.array([0.0, 1.0, 0.0])
        else:
            ref = np.array([0.0, 0.0, 1.0])
            e1 = np.cross(ref, axis)
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(axis, e1)

        t = np.linspace(0.0, 2.0 * np.pi, n_samples)
        pts3d = np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2)

        upper = pts3d[:, 2] >= -1e-9
        arcs, current = [], []
        for keep, p in zip(upper, pts3d):
            if keep:
                current.append(self._project(p))
            elif current:
                arcs.append(np.array(current))
                current = []
        if current:
            arcs.append(np.array(current))
        # stitch the arc that wraps across the t=0 / t=2*pi seam
        if len(arcs) >= 2 and upper[0] and upper[-1]:
            arcs[0] = np.vstack([arcs[-1], arcs[0]])
            arcs.pop()
        return arcs

    # ------------------------------------------------------------------ #
    # plotting
    # ------------------------------------------------------------------ #
    def plot(self, ax=None, show_labels=True, show_boundary=True):
        created_fig = ax is None
        if created_fig:
            fig, ax = plt.subplots(figsize=(self.figsize, self.figsize))

        ax.set_aspect('equal')
        ax.axis('off')

        if show_boundary:
            ax.add_patch(Circle((0, 0), 1.0, fill=False, color='black',
                                 linewidth=1.3))
            ax.plot(0, 0, '+', color='black', markersize=8, markeredgewidth=1.2)

        for z in self.zones:
            arcs = self._great_circle_arcs(z['axis'], z['n_samples'])
            for arc in arcs:
                ax.plot(arc[:, 0], arc[:, 1], color=z['color'],
                        linestyle=z['linestyle'], linewidth=z['linewidth'])
            if show_labels and arcs:
                mid = arcs[0][len(arcs[0]) // 2]
                ax.annotate(z['label'], mid, color=z['color'], fontsize=9,
                            xytext=(3, 3), textcoords='offset points')

        for p in self.poles:
            X, Y = self._project(p['vec'])
            fc = p['color'] if p['fill'] else 'none'
            ax.scatter([X], [Y], s=p['size'], marker=p['marker'],
                       facecolor=fc, edgecolor=p['color'], zorder=5)
            if show_labels:
                ax.annotate(p['label'], (X, Y), fontsize=9,
                            xytext=(5, 5), textcoords='offset points')

        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)

        if created_fig:
            plt.tight_layout()
        return ax


if __name__ == "__main__":
    # Sanity-check example: cubic poles and a couple of zones.
    sp = StereoProjector()
    sp.add_pole(1, 0, 0)
    sp.add_pole(0, 1, 0)
    #sp.add_pole(0, 0, 1)
    sp.add_pole(-1, 0, 0)
    sp.add_pole(0, -1, 0)
    sp.add_pole(0, 0, -1)
    #sp.add_pole(1, 1, 1)
    #sp.add_pole(1, 1, 0)
    #sp.add_zone(0, 0, 1)     # equator (boundary circle)
    #sp.add_zone(1, 0, 0)     # vertical great circle through the center
    #sp.add_zone(1, 1, 0)     # tilted great circle
    sp.plot()
    plt.show()
    #plt.savefig("stereoproj_demo.png", dpi=150)
    print("saved stereoproj_demo.png")