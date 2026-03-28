# dist_plots.py
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt
import warnings
from dist_formulas import (
    calc_rectifying_ol_params,
    calc_feedline_intercept,
    calc_operating_feedline_intersection,
    calc_feedline_equilibrium_intersection,
    calc_minimum_reflux,
    run_mccabe_thiele_stepper,
    calc_murphree_pseudo_equilibrium,
    build_pseudo_equil_callable
)

def build_equilibrium_callable(x_data=None, y_data=None, alpha=None):
    """
    Returns a callable f(x) -> y representing the equilibrium curve.
    Accepts either paired x-y data or alpha value(s).

    Parameters
    ----------
    x_data : list or array, optional
        Liquid phase mole fractions of light component.
    y_data : list or array, optional
        Vapor phase mole fractions of light component.
    alpha : float or list of float, optional
        Relative volatility. Accepts:
            - Single value  [alpha]                          — constant alpha
            - Two values    [alpha_top, alpha_bottom]        — geometric mean
            - Three values  [alpha_top, alpha_mid, alpha_bottom] — cubic mean

    Returns
    -------
    callable
        Function that takes x (float or array) and returns y.
    """
    if x_data is not None and y_data is not None:
        # Drop NaN rows that data_editor appends for new row entry
        x_arr = np.array(x_data, dtype=float)
        y_arr = np.array(y_data, dtype=float)
        valid = ~(np.isnan(x_arr) | np.isnan(y_arr))
        x_arr = x_arr[valid]
        y_arr = y_arr[valid]

        if len(x_arr) < 2:
            raise ValueError("At least 2 valid data points are required.")
        if not all(0 <= x <= 1 for x in x_arr):
            raise ValueError("All x_data values must be between 0 and 1.")
        if not all(0 <= y <= 1 for y in y_arr):
            raise ValueError("All y_data values must be between 0 and 1.")
        x_arr = np.array(x_data)
        y_arr = np.array(y_data)
        if x_arr[0] != 0:
            x_arr = np.insert(x_arr, 0, 0.0)
            y_arr = np.insert(y_arr, 0, 0.0)
        if x_arr[-1] != 1:
            x_arr = np.append(x_arr, 1.0)
            y_arr = np.append(y_arr, 1.0)
        return CubicSpline(x_arr, y_arr)

    elif alpha is not None:
        # Handle single float or list uniformly
        if isinstance(alpha, (int, float)):
            alpha_values = [alpha]
        else:
            alpha_values = list(alpha)

        if any(a < 1 for a in alpha_values):
            raise ValueError("All relative volatility values must be >= 1.")
        if len(alpha_values) > 3:
            warnings.warn(
                f"{len(alpha_values)} alpha values provided. "
                "Standard practice uses 1-3 values. Verify this is intentional.",
                UserWarning
            )

        # Geometric mean regardless of how many values provided
        n = len(alpha_values)
        alpha_mean = np.prod(alpha_values) ** (1 / n)

        return lambda x: (alpha_mean * x) / (1 + (alpha_mean - 1) * x)

    else:
        raise ValueError("Must provide either x_data/y_data or alpha.")


def plot_equilibrium_curve(ax, equil_func, n_points=200):
        """
        Plots the equilibrium curve on the provided axes.

        Parameters
        ----------
        ax : matplotlib Axes
        equil_func : callable
            Equilibrium curve function from build_equilibrium_callable.
        n_points : int
            Number of points for smooth curve rendering.
        """
        x_vals = np.linspace(0, 1, n_points)
        y_vals = equil_func(x_vals)
        ax.plot(x_vals, y_vals, color="blue", linewidth=1.5, label="Equilibrium curve")

def plot_pseudo_equilibrium_curve(ax, x_vals, y_pseudo):
    """
    Plots the Murphree pseudo-equilibrium curve on the provided axes.
    Sits between the operating lines and the true equilibrium curve,
    representing actual separation achieved per tray.

    Parameters
    ----------
    ax : matplotlib Axes
    x_vals : array
        x values from calc_murphree_pseudo_equilibrium.
    y_pseudo : array
        y values from calc_murphree_pseudo_equilibrium.
    """
    ax.plot(x_vals, y_pseudo, color="purple", linewidth=1.0,
            linestyle="--", label="Pseudo-equilibrium curve")

def plot_diagonal(ax):
    """
    Plots the y = x diagonal (45-degree line) on the provided axes.
    Represents total reflux and serves as the reference line for
    McCabe-Thiele stage stepping.

    Parameters
    ----------
    ax : matplotlib Axes
    """
    ax.plot([0, 1], [0, 1], color="black", linewidth=1.0,
            linestyle="--", label="y = x")

def plot_feedline(ax, f: float, z_f: float, x_range=(0, 1)):
    info = calc_feedline_intercept(f, z_f)
    line_type = info["line_type"]
    label = "Feed line"

    if line_type == "vertical":
        ax.axvline(x=info["x_intercept"], linestyle="--", color="green", label=label)

    elif line_type == "horizontal":
        ax.axhline(y=info["y_intercept"], linestyle="--", color="green", label=label)

    else:  # normal
        x_min, x_max = x_range
        x_vals = np.linspace(x_min, x_max, 200)
        y_vals = (z_f - (1 - f) * x_vals) / f
        mask = (y_vals >= x_min) & (y_vals <= x_max)
        ax.plot(x_vals[mask], y_vals[mask], linestyle="--", color="green", label=label)


def plot_rectifying_operating_line(ax, x_d, x_i, y_i, reflux_ratio):
    """
    Plots the rectifying section operating line from (x_I, y_I) to (x_D, x_D).

    Parameters
    ----------
    ax : matplotlib Axes
    x_d : float
        Mole fraction of light component in the distillate (0 to 1).
    x_i : float
        x-coordinate of q-line / rectifying OL intersection.
    y_i : float
        y-coordinate of q-line / rectifying OL intersection.
    reflux_ratio : float
        Reflux ratio R = L/D. Must be >= 0.
    """
    x_vals = np.array([x_i, x_d])
    y_vals = np.array([y_i, x_d])

    ax.plot(x_vals, y_vals, color="red", linewidth=1.5,
            label=f"Rectifying OL (R={reflux_ratio:.2f})")

def plot_stripping_operating_line(ax, x_b, x_i, y_i):
    """
    Plots the stripping section operating line from (x_B, x_B) to (x_I, y_I).

    Parameters
    ----------
    ax : matplotlib Axes
    x_b : float
        Mole fraction of light component in the bottoms (0 to 1).
    x_i : float
        x-coordinate of q-line / rectifying OL intersection.
    y_i : float
        y-coordinate of q-line / rectifying OL intersection.

    Returns
    -------
    tuple
        (slope, x_intercept) for use by coordinator and stepper.
    """
    if x_b < 0 or x_b > 1:
        raise ValueError("Bottoms mole fraction must be between 0 and 1.")
    if not 0 <= x_i <= 1 or not 0 <= y_i <= 1:
        raise ValueError("Intersection point must be within [0,1] bounds.")
    if x_i <= x_b:
        raise ValueError(
            "Intersection point x_I must be greater than x_B. "
            "Check feed condition and operating line inputs."
        )

    # Stripping OL runs from (x_B, x_B) on the diagonal up to (x_I, y_I)
    x_vals = np.array([x_b, x_i])
    y_vals = np.array([x_b, y_i])

    # Calculate slope for stepper use
    slope = (y_i - x_b) / (x_i - x_b)

    ax.plot(x_vals, y_vals, color="orange", linewidth=1.5,
            label="Stripping OL")

    return slope, x_b

def plot_minimum_reflux_line(ax, x_d, x_prime, y_prime):
    """
    Plots the minimum reflux operating line as a visual reference.
    Runs from (x_D, x_D) to the feed line / equilibrium curve intersection.

    Parameters
    ----------
    ax : matplotlib Axes
    x_d : float
        Distillate mole fraction.
    x_prime : float
        x-coordinate of feed line / equilibrium curve intersection.
    y_prime : float
        y-coordinate of feed line / equilibrium curve intersection.
    """
    r_min = calc_minimum_reflux(x_d, x_prime, y_prime)

    x_vals = np.array([x_prime, x_d])
    y_vals = np.array([y_prime, x_d])

    ax.plot(x_vals, y_vals, color="darkred", linewidth=1.0,
            linestyle=":", label=f"Min reflux OL (R_min={r_min:.2f})")

    return r_min

def plot_mccabe_thiele(ax, x_d, x_b, f, z_f, reflux_ratio,
                           alpha, equil_func,
                           condenser_type="total",
                           reboiler_type="partial",
                           murphree_efficiency=1.0):
    """
    Coordinator function for the McCabe-Thiele diagram.
    Calls individual plot functions in dependency order and passes
    calculated intermediate values between them.

    Parameters
    ----------
    ax : matplotlib Axes
    x_d : float
        Mole fraction of light component in the distillate (0 to 1).
    x_b : float
        Mole fraction of light component in the bottoms (0 to 1).
    f : float
        Mole fraction of vapor in the feed.
    z_f : float
        Feed composition, mole fraction of light component (0 to 1).
    reflux_ratio : float
        Reflux ratio R = L/D.
    equil_func : callable
        Equilibrium curve function from build_equilibrium_callable.

    Returns
    -------
    dict
        Intermediate calculated values for use by stepper and Streamlit display.
    """
    # --- Layer 1: Reference lines, no dependencies ---
    plot_equilibrium_curve(ax, equil_func)
    plot_diagonal(ax)

    # --- Layer 2: Feed line ---
    plot_feedline(ax, f, z_f)

    # --- Layer 3: Rectifying OL ---
    # Coordinator calculates x_i, y_i internally...
    rect_slope, rect_y_int = calc_rectifying_ol_params(x_d, reflux_ratio)
    x_i, y_i = calc_operating_feedline_intersection(rect_slope, rect_y_int, f, z_f)

    # ...then passes them down to the plot function
    plot_rectifying_operating_line(ax, x_d, x_i, y_i, reflux_ratio)

    # --- Layer 3b: Minimum reflux reference line ---
    # --- Layer 3b: Minimum reflux reference line ---
    if alpha is not None:
        n = len(alpha)
        alpha_mean = np.prod(alpha) ** (1 / n)
        x_prime, y_prime = calc_feedline_equilibrium_intersection(alpha_mean, f, z_f)
        r_min = plot_minimum_reflux_line(ax, x_d, x_prime, y_prime)
    else:
        r_min = None  # R_min not available without alpha

    # --- Layer 4: Feed line / rectifying OL intersection ---
    x_i, y_i = calc_operating_feedline_intersection(
                    rect_slope, rect_y_int, f, z_f)

    # --- Layer 5: Stripping OL, anchored to intersection and x_B ---
    strip_slope, strip_x_int = plot_stripping_operating_line(ax, x_b, x_i, y_i)

    # --- Layer 5b: Murphree pseudo-equilibrium curve ---
    def rect_ol_func(x):
        return rect_slope * x + rect_y_int

    def strip_ol_func(x):
        return strip_slope * (x - x_b) + x_b

    if murphree_efficiency < 1.0:
        x_pseudo, y_pseudo = calc_murphree_pseudo_equilibrium(
            equil_func, rect_ol_func, strip_ol_func,
            x_i, murphree_efficiency
        )
        pseudo_callable = build_pseudo_equil_callable(x_pseudo, y_pseudo)
        plot_pseudo_equilibrium_curve(ax, x_pseudo, y_pseudo)
        active_equil_func = pseudo_callable
    else:
        active_equil_func = equil_func

    # --- Layer 6: Stepper --- uses active_equil_func, not equil_func directly
    stepper_results = run_mccabe_thiele_stepper(
        x_d, x_b, f, z_f,
        rect_slope, rect_y_int,
        strip_slope, x_b,
        active_equil_func,  # ← pseudo-curve when efficiency < 1, else true curve
        x_i,
        condenser_type=condenser_type,
        reboiler_type=reboiler_type
    )
    plot_staircase(ax, stepper_results["stage_coords"])

    # --- Axes formatting --- (must be BEFORE return)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("x (liquid mole fraction)")
    ax.set_ylabel("y (vapor mole fraction)")
    ax.set_title("McCabe-Thiele Diagram")
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.12),
              ncol=3, frameon=True)
    ax.grid(True, linestyle="--", alpha=0.4)

    # --- Return intermediate values ---
    return {
        "rect_slope": rect_slope,
        "rect_y_int": rect_y_int,
        "strip_slope": strip_slope,
        "strip_x_int": x_b,
        "x_i": x_i,
        "y_i": y_i,
        "r_min": r_min,
        "total_stages": stepper_results["total_stages"],
        "feed_stage": stepper_results["feed_stage"],
        "rectifying_stages": stepper_results["rectifying_stages"],
        "stripping_stages": stepper_results["stripping_stages"],
    }

def plot_staircase(ax, stage_coords):
    """
    Plots the McCabe-Thiele staircase from stepper output.

    Parameters
    ----------
    ax : matplotlib Axes
    stage_coords : list of (x, y) tuples
        Coordinate pairs defining the staircase path.
    """
    x_vals = [coord[0] for coord in stage_coords]
    y_vals = [coord[1] for coord in stage_coords]
    ax.plot(x_vals, y_vals, color="black", linewidth=1.0, label="Stages")

