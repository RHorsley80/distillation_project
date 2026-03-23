# app.py
import streamlit as st
import matplotlib.pyplot as plt
from dist_plots import plot_mccabe_thiele, build_equilibrium_callable
import pandas as pd
from dist_formulas import (
    describe_feed_condition,
    calc_feedline_equilibrium_intersection,
    calc_minimum_reflux,
    fenske_min_theoretical_stages,
    kirkbride_rectifying_stripping_tray_ratio,
    kirkbride_tray_counts
)
import numpy as np

# --- User inputs (sidebar) ---
st.sidebar.header("System")

st.sidebar.header("Equilibrium Data")
data_source = st.sidebar.radio(
    "Equilibrium curve source",
    ["Constant alpha", "Upload x-y data"]
)

alpha_values = None  # default before branching

if data_source == "Constant alpha":
    alpha_top = st.sidebar.number_input("Alpha (top)", min_value=1.0, value=2.5)
    alpha_bottom = st.sidebar.number_input("Alpha (bottom)", min_value=1.0, value=2.0)
    alpha_values = [alpha_top, alpha_bottom]
    equil_func = build_equilibrium_callable(alpha=alpha_values)

else:
    uploaded_file = st.sidebar.file_uploader("Upload CSV (columns: x, y)", type="csv")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        if "x" not in df.columns or "y" not in df.columns:
            st.error("CSV must have columns named 'x' and 'y'.")
            st.stop()
        equil_func = build_equilibrium_callable(
            x_data=df["x"].tolist(),
            y_data=df["y"].tolist()
        )
    else:
        st.sidebar.markdown("**Or enter data manually:**")
        default_df = pd.DataFrame({"x": [0.0, 0.5, 1.0], "y": [0.0, 0.7, 1.0]})
        edited_df = st.sidebar.data_editor(
            default_df,
            num_rows="dynamic",
            use_container_width=True
        )
        equil_func = build_equilibrium_callable(
            x_data=edited_df["x"].tolist(),
            y_data=edited_df["y"].tolist()
        )

st.sidebar.header("Separation Targets")
x_d = st.sidebar.number_input("x_D (distillate)", min_value=0.0, max_value=1.0, value=0.95)
x_b = st.sidebar.number_input("x_B (bottoms)", min_value=0.0, max_value=1.0, value=0.05)
z_f = st.sidebar.number_input("z_F (feed composition)", min_value=0.0, max_value=1.0, value=0.50)

st.sidebar.header("Feed Condition")

feed_input_type = st.sidebar.radio(
    "Feed condition input",
    ["Enter f (vapor fraction)", "Enter q (liquid fraction)"]
)

if feed_input_type == "Enter f (vapor fraction)":
    f = st.sidebar.number_input(
        "f (vapor fraction in feed)",
        min_value=-0.5, max_value=1.5, value=0.5
    )
else:
    q = st.sidebar.number_input(
        "q (liquid fraction in feed)",
        min_value=-0.5, max_value=1.5, value=0.5
    )
    f = 1 - q

# Display feed condition description
feed_info = describe_feed_condition(f)
st.sidebar.info(
    f"**{feed_info['condition']}**\n\n{feed_info['description']}"
)

st.sidebar.header("Operating Parameters")
reflux_ratio = st.sidebar.number_input("Reflux Ratio (R)", min_value=0.0, value=2.0)

st.sidebar.header("Condenser / Reboiler")
condenser_type = st.sidebar.selectbox("Condenser type", ["total", "partial"])
reboiler_type = st.sidebar.selectbox("Reboiler type", ["partial", "total"])

st.sidebar.header("Molar Flow Rates")
dist_molar_drawrate = st.sidebar.number_input(
    "Distillate molar draw rate (D)", min_value=0.01, value=1.0)
bottoms_molar_drawrate = st.sidebar.number_input(
    "Bottoms molar draw rate (B)", min_value=0.01, value=1.0)

# --- Pre-calculate R_min for validation (constant alpha path only) ---
r_min_check = None
if data_source == "Constant alpha" and alpha_values is not None:
    n = len(alpha_values)
    alpha_mean = float(np.prod(alpha_values) ** (1 / n))
    try:
        x_prime, y_prime = calc_feedline_equilibrium_intersection(alpha_mean, f, z_f)
        r_min_check = calc_minimum_reflux(x_d, x_prime, y_prime)
    except Exception:
        r_min_check = None

if r_min_check is not None and reflux_ratio <= r_min_check:
    st.error(
        f"Reflux ratio R = {reflux_ratio:.3f} is at or below "
        f"R_min = {r_min_check:.3f}. Separation is not achievable. "
        "Increase the reflux ratio."
    )
    st.stop()

# --- Build figure and axes ---
fig, ax = plt.subplots(figsize=(7, 7))

# --- Call coordinator, receive intermediate values ---
results = plot_mccabe_thiele(
    ax, x_d, x_b, f, z_f, reflux_ratio,
    alpha=alpha_values if data_source == "Constant alpha" else None,
    equil_func=equil_func,
    condenser_type=condenser_type,
    reboiler_type=reboiler_type
)

# --- Display diagram ---
# --- Display diagram ---
st.pyplot(fig)

# --- Display intermediate values ---
st.subheader("Operating Line Intersection")
col1, col2 = st.columns(2)
col1.metric("x_I", f"{results['x_i']:.4f}")
col2.metric("y_I", f"{results['y_i']:.4f}")

st.subheader("Minimum Reflux")
if data_source == "Constant alpha" and results['r_min'] is not None:
    st.metric("R_min", f"{results['r_min']:.4f}")
else:
    st.info("R_min not available — requires constant alpha input.")

# --- Stage results ---
st.subheader("Theoretical Stages")
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total Stages", results["total_stages"])
col2.metric("Feed Stage", results["feed_stage"])
col3.metric("Rectifying Stages", results["rectifying_stages"])
col4.metric("Stripping Stages", results["stripping_stages"])

# --- Fenske Minimum Stages ---
st.subheader("Fenske Equation (Minimum Stages at Total Reflux)")
if data_source == "Constant alpha" and alpha_values is not None:
    n_min = fenske_min_theoretical_stages(x_d, x_b, alpha_values)
    #st.metric("N_min (Fenske)", f"{n_min:.2f}")
    col1, col2 = st.columns(2)
    col1.metric("N_min (Fenske, total reflux)", f"{n_min:.2f}")
    col2.metric(f"N (McCabe-Thiele, R={reflux_ratio:.2f})",
                results["total_stages"])
else:
    st.info("Fenske equation requires constant alpha input.")

# --- Kirkbride Feed Tray Location ---
st.subheader("Kirkbride Equation (Feed Tray Location)")
from dist_formulas import kirkbride_rectifying_stripping_tray_ratio, kirkbride_tray_counts

# Binary system simplification
z_hk = 1 - z_f
z_lk = z_f
x_blk = x_b
x_dhk = 1 - x_d
st.caption(
    "Note: McCabe-Thiele stage counts are determined graphically by the "
    "stepping algorithm. Kirkbride provides an independent estimate based "
    "on feed and product compositions and molar draw rates. Discrepancies "
    "between the two methods are normal and reflect different underlying "
    "assumptions. Engineering judgment should be applied when selecting "
    "the feed tray location for detailed design."
)

try:
    nr_ns_ratio = kirkbride_rectifying_stripping_tray_ratio(
        z_hk, z_lk, x_blk, x_dhk,
        bottoms_molar_drawrate, dist_molar_drawrate
    )
    kirkbride_counts = kirkbride_tray_counts(
        results["total_stages"], nr_ns_ratio
    )
    col1, col2, col3 = st.columns(3)
    col1.metric("N_R / N_S ratio", f"{nr_ns_ratio:.3f}")
    col2.metric("Rectifying Stages (Kirkbride)",
                f"{kirkbride_counts['N_R']:.1f}")
    col3.metric("Stripping Stages (Kirkbride)",
                f"{kirkbride_counts['N_S']:.1f}")
except ValueError as e:
    st.warning(f"Kirkbride calculation not available: {e}")