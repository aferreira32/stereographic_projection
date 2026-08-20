"""
app.py
======

A Streamlit front-end for stereoproj.StereoProjector: add poles (hkl) and
zones/great circles [uvw] one at a time, see them appear on a live
stereographic projection, and undo/clear as needed.

Run with:  streamlit run app.py
"""

import itertools

import matplotlib.pyplot as plt
import streamlit as st

from stereoproj import StereoProjector

st.set_page_config(page_title="Stereographic Projection", page_icon="\U0001F53A",
                    layout="centered")

if "history" not in st.session_state:
    # Ordered list of everything added so far, poles and zones interleaved,
    # so "undo" can remove whichever came last regardless of its kind.
    st.session_state.history = []

st.title("Stereographic Projection (Wulff Net)")
st.caption(
    "Build up a stereographic projection by adding poles and zones. "
    "Every plot starts from the empty projection circle."
)

with st.expander("What am I looking at?"):
    st.markdown(
        r"""
- A **pole** $(hkl)$ is a single point: the direction normal to the
  $(hkl)$ plane. Poles in the **upper hemisphere** ($z \ge 0$) are drawn
  as **filled** dots; poles in the **lower hemisphere** are drawn as
  **open** rings, so a coincident pair like $(hkl)$/$(hk\bar l)$ shows up
  as a ring enclosing a dot instead of one hiding the other.
- A **zone** $[uvw]$ is a direction (the zone axis). The tool draws the
  great circle $90°$ away from that axis -- the locus of poles of every
  plane whose normal is perpendicular to $[uvw]$, i.e. every plane that
  direction lies in.
- Miller indices are currently treated as direct Cartesian components,
  which is exact for a cubic lattice and an approximation otherwise.
        """
    )

col_pole, col_zone = st.columns(2)

with col_pole:
    st.subheader("Add a pole (hkl)")
    with st.form("pole_form", clear_on_submit=True):
        ph, pk, pl = st.columns(3)
        h = ph.number_input("h", value=None, step=1, format="%d",
                             placeholder="h", label_visibility="collapsed")
        k = pk.number_input("k", value=None, step=1, format="%d",
                             placeholder="k", label_visibility="collapsed")
        l = pl.number_input("l", value=None, step=1, format="%d",
                             placeholder="l", label_visibility="collapsed")
        add_pole = st.form_submit_button("Add pole", use_container_width=True)
    if add_pole:
        if h is None or k is None or l is None:
            st.error("Enter values for h, k, and l.")
        elif h == 0 and k == 0 and l == 0:
            st.error("(0 0 0) is not a valid pole.")
        else:
            st.session_state.history.append({"kind": "pole", "h": h, "k": k, "l": l})

with col_zone:
    st.subheader("Add a zone [uvw]")
    with st.form("zone_form", clear_on_submit=True):
        zu, zv, zw = st.columns(3)
        u = zu.number_input("u", value=None, step=1, format="%d",
                             placeholder="u", label_visibility="collapsed")
        v = zv.number_input("v", value=None, step=1, format="%d",
                             placeholder="v", label_visibility="collapsed")
        w = zw.number_input("w", value=None, step=1, format="%d",
                             placeholder="w", label_visibility="collapsed")
        add_zone = st.form_submit_button("Add zone", use_container_width=True)
    if add_zone:
        if u is None or v is None or w is None:
            st.error("Enter values for u, v, and w.")
        elif u == 0 and v == 0 and w == 0:
            st.error("[0 0 0] is not a valid zone axis.")
        else:
            st.session_state.history.append({"kind": "zone", "u": u, "v": v, "w": w})

undo_col, clear_col = st.columns(2)
with undo_col:
    if st.button("Undo last", use_container_width=True,
                  disabled=not st.session_state.history):
        st.session_state.history.pop()
with clear_col:
    if st.button("Clear all", use_container_width=True,
                  disabled=not st.session_state.history):
        st.session_state.history.clear()

# ---------------------------------------------------------------------- #
# Rebuild the projection from the current history and draw it.
# ---------------------------------------------------------------------- #
sp = StereoProjector()
palette = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
entries = []  # (label_html, color) for the sidebar-style summary list

for item in st.session_state.history:
    color = next(palette)
    if item["kind"] == "pole":
        sp.add_pole(item["h"], item["k"], item["l"], color=color)
        entries.append((sp.poles[-1]["label"], color, "pole"))
    else:
        sp.add_zone(item["u"], item["v"], item["w"], color=color)
        entries.append((sp.zones[-1]["label"], color, "zone"))

fig, ax = plt.subplots(figsize=(6.5, 6.5))
sp.plot(ax=ax)
st.pyplot(fig)
plt.close(fig)

if entries:
    st.subheader("Currently plotted")
    for label, color, kind in entries:
        dot_col, label_col = st.columns([1, 20])
        with dot_col:
            st.markdown(
                f'<div style="width:12px;height:12px;border-radius:50%;'
                f'background:{color};margin-top:8px;"></div>',
                unsafe_allow_html=True,
            )
        with label_col:
            # Plain markdown (no raw HTML mixed in) so the $...$ mathtext
            # in the label renders as LaTeX instead of showing literally.
            st.markdown(f"{label} — {kind}")
else:
    st.info("Add a pole or a zone above to get started.")
