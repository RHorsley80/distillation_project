import numpy as np
import warnings

def calc_rel_volatility_alpha(mole_frac_in_vapor, mole_frac_in_liquid):
    #Flag errors.
    if mole_frac_in_vapor > 1 or mole_frac_in_vapor < 0:
        raise ValueError("The mole fraction of vapor must between zero and one.")
    if mole_frac_in_liquid > 1 or mole_frac_in_liquid < 0:
        raise ValueError("The mole fraction of liquid must between zero and one.")
    # calculate relative volatility alpha for a substance
    alpha = (mole_frac_in_vapor/mole_frac_in_liquid)/((1-mole_frac_in_vapor)/(1-mole_frac_in_liquid))
    return  alpha

#The following equation is also the equilibrium line equation.
def calc_mole_frac_in_vapor_phase(rel_volatility_alpha, mole_frac_in_liquid):
    #Flag errors.
    if rel_volatility_alpha < 1:
        raise ValueError("Relative volatility must be >= 1 (alpha < 1 means no separation is possible).")
    if mole_frac_in_liquid > 1 or mole_frac_in_liquid < 0:
        raise ValueError("The mole fraction of liquid must between zero and one.")
    #Solve for y, the mole fraction of the light component in the vapor phase
    y = (rel_volatility_alpha*mole_frac_in_liquid)/(1+(rel_volatility_alpha-1)*mole_frac_in_liquid)
    return y

#Equations for the Rectifying Section
def calc_rectifying_slope(reflux_ratio):
    if reflux_ratio < 0:
        raise ValueError("The reflux ratio cannot be negative.")
    # R = 0 is a valid edge case (no reflux, slope = 0, horizontal OL)
    return reflux_ratio / (reflux_ratio + 1)

def calc_rectifying_y_intercept(distillate_mole_frac_liquid, reflux_ratio):
    #Flag errors.
    if reflux_ratio < 0:
        raise ValueError("The reflux ratio cannot be negative.")
    if distillate_mole_frac_liquid <0 or distillate_mole_frac_liquid >1:
        raise ValueError("The liquid mole fraction of x in the distillate must be between zero and one.")
    #Solve for y intercept
    y_int = distillate_mole_frac_liquid/(reflux_ratio+1)
    return y_int

def calc_reflux_ratio(reflux_molrate_to_column, distillate_draw_molrate):
    if reflux_molrate_to_column < 0:
        raise ValueError("The reflux rate to the column cannot be negative.")
    if distillate_draw_molrate <= 0:
        raise ValueError("The distillate draw rate must be greater than zero.")
    return reflux_molrate_to_column / distillate_draw_molrate

#Equations for the Stripping Section
def calc_stripping_slope(boilup_ratio):
    if boilup_ratio < 0:
        raise ValueError("The boilup ratio cannot be negative.")
    if boilup_ratio == 0:
        # No boilup: stripping line is vertical, slope undefined
        return float('inf')
    return (boilup_ratio + 1) / boilup_ratio

def calc_stripping_x_intercept(bottoms_mole_frac_liquid, boilup_ratio):
    if bottoms_mole_frac_liquid < 0 or bottoms_mole_frac_liquid > 1:
        raise ValueError("The mole fraction of liquid in the bottoms must be between zero and one.")
    if boilup_ratio < 0:
        raise ValueError("The boilup ratio cannot be negative.")
    if boilup_ratio == 0:
        # No boilup: stripping line is vertical at x = x_B, intercept undefined in normal sense
        return None
    return (bottoms_mole_frac_liquid + boilup_ratio) / (boilup_ratio + 1)

def calc_boilup_ratio(boilup_molrate_to_column, bottoms_draw_molrate):
    if boilup_molrate_to_column < 0:
        raise ValueError("The boilup rate cannot be negative.")
    if bottoms_draw_molrate <= 0:
        raise ValueError("The bottoms draw rate must be greater than zero.")
    return boilup_molrate_to_column / bottoms_draw_molrate

#Equations for the feed line.
def calc_feedline_slope(mole_frac_feed_vapor):
    if mole_frac_feed_vapor < -0.5 or mole_frac_feed_vapor > 1.5:
        warnings.warn(            f"f = {mole_frac_feed_vapor:.3f} is outside the expected physical range "
            "[-0.5, 1.5]. Verify feed condition.",
            UserWarning)
    if mole_frac_feed_vapor == 0:
        return float('inf')   # Vertical line, slope is undefined
    if mole_frac_feed_vapor == 1:
        return 0.0            # Horizontal line
    return (mole_frac_feed_vapor - 1) / mole_frac_feed_vapor

def calc_feedline_intercept(f: float, z_f: float) -> dict:
    """
    Calculate the feed line (q-line) intercept(s) for a McCabe-Thiele diagram.

    Parameters
    ----------
    f : float
        Mole fraction of vapor in the feed (-0.5 to 1.5).
    z_f : float
        Mole fraction of the more volatile component in the feed (0 to 1).

    Returns
    -------
    dict with keys:
        'line_type'   : str   — 'vertical', 'horizontal', or 'normal'
        'x_intercept' : float or None
        'y_intercept' : float or None
    """
    if not 0 <= f <= 1:
        warnings.warn(f"f = {f:.3f} is outside the expected physical range [-0.5, 1.5]. "
            "Subcooled (f<0) and superheated (f>1) feeds are supported but verify input.",
            UserWarning)
    if not 0 < z_f <= 1:
        raise ValueError("Feed component mole fraction (z_f) must be between 0 and 1.")

    # --- Degenerate cases ---
    if f == 0:
        # Saturated liquid: q-line is vertical at x = z_f, intercepts are undefined
        return {"line_type": "vertical", "x_intercept": z_f, "y_intercept": None}

    if f == 1:
        # Saturated vapor: q-line is horizontal at y = z_f, intercepts are undefined
        return {"line_type": "horizontal", "x_intercept": None, "y_intercept": z_f}

    # --- Normal case ---
    result = {"line_type": "normal", "x_intercept": None, "y_intercept": None}

    if z_f >= (1 - f):
        result["x_intercept"] = z_f / (1 - f)

    if z_f <= (1 - f):
        result["y_intercept"] = (z_f + f - 1) / f

    return result

#Equations to find minimum reflux

#TODO: Find reference problem(s) to validate the calc_feedline_equilibrium_intersection formula.

def calc_feedline_equilibrium_intersection(alpha, f, z_f):
    """
    Finds (x', y') where the q-line intersects the equilibrium curve.
    Needed for minimum reflux calculation.

    Parameters
    ----------
    alpha : float
        Relative volatility (geometric mean if multiple values).
    f : float
        Mole fraction of vapor in feed.
    z_f : float
        Feed composition, mole fraction of light component.

    Returns
    -------
    tuple
        (x_prime, y_prime) intersection coordinates.
    """
    if f == 0:
        # Vertical q-line: x_prime = z_f, y_prime from equilibrium curve
        x_prime = z_f
        y_prime = (alpha * x_prime) / (1 + (alpha - 1) * x_prime)
        return x_prime, y_prime

    if f == 1:
        # Horizontal q-line: y_prime = z_f, solve equilibrium curve for x_prime
        y_prime = z_f
        x_prime = y_prime / (alpha - (alpha - 1) * y_prime)
        return x_prime, y_prime

    # Normal case — quadratic solution
    A_term1 = 1 / (alpha - 1)
    A_term2 = z_f / (f - 1)
    A_term3 = (alpha * f) / ((alpha - 1) * (f - 1))
    A = A_term1 + A_term2 - A_term3

    B = z_f / ((alpha - 1) * (f - 1))

    x_prime = -0.5 * A + (0.25 * A ** 2 - B) ** 0.5
    y_prime = x_prime * ((f - 1) / f) + (z_f / f)

    return x_prime, y_prime

def calc_minimum_reflux(x_D, x_prime, y_prime):
    return (x_D - y_prime) / (y_prime - x_prime)

def describe_feed_condition(f: float) -> dict:
    """
    Identifies the feed thermal condition from the vapor fraction f.

    Parameters
    ----------
    f : float
        Mole fraction of vapor in the feed. May be negative (subcooled)
        or greater than 1 (superheated).

    Returns
    -------
    dict with keys:
        'condition'   : str  — short label
        'description' : str  — operational meaning
        'q'           : float — q value (q = 1 - f)
        'line_type'   : str  — matches calc_feedline_intercept output
    """
    q = 1 - f

    if f < 0:
        condition = "Subcooled Liquid"
        description = ("Feed is below its bubble point. Column must supply heat to reach "
                       "saturation before vaporization begins. Feed tray sees increased "
                       "liquid traffic — raises internal L/V ratio in stripping section.")
        line_type = "normal"

    elif f == 0:
        condition = "Saturated Liquid (Bubble Point)"
        description = ("Feed enters exactly at its bubble point. Most common design basis. "
                       "Feed tray liquid and vapor loads are well-defined.")
        line_type = "vertical"

    elif 0 < f < 1:
        condition = "Partially Vaporized"
        description = (f"Feed is a two-phase mixture ({f*100:.1f}% vapor). "
                       "Both liquid and vapor traffic increase at the feed tray. "
                       "Feed tray design must accommodate mixed phase entry.")
        line_type = "normal"

    elif f == 1:
        condition = "Saturated Vapor (Dew Point)"
        description = ("Feed enters exactly at its dew point. Vapor traffic increases "
                       "in rectifying section. Feed tray sees no additional liquid load.")
        line_type = "horizontal"

    else:  # f > 1
        condition = "Superheated Vapor"
        description = ("Feed is above its dew point. Excess heat condenses liquid on "
                       "the feed tray, reducing liquid flow in stripping section. "
                       "Can significantly impact tray hydraulics near feed point.")
        line_type = "normal"

    return {
        "condition": condition,
        "description": description,
        "q": q,
        "line_type": line_type
    }

#Theoretical stage formulas
def get_condenser_stage_credit(condenser_type: str) -> int:
    """Returns number of theoretical stages credited to the condenser."""
    options = {"total": 0, "partial": 1}
    if condenser_type not in options:
        raise ValueError(f"condenser_type must be 'total' or 'partial', got '{condenser_type}'")
    return options[condenser_type]

def get_reboiler_stage_credit(reboiler_type: str) -> int:
    """Returns number of theoretical stages credited to the reboiler."""
    options = {"partial": 1, "total": 0}
    if reboiler_type not in options:
        raise ValueError(f"reboiler_type must be 'partial' or 'total', got '{reboiler_type}'")
    return options[reboiler_type]

def fenske_min_theoretical_stages(x_d: float, x_b: float,
                                   alpha_values: list) -> float:
    """
    Calculate the minimum number of theoretical stages at total reflux (Fenske equation).

    Parameters
    ----------
    x_d : float
        Mole fraction of light component in the distillate (0 to 1).
    x_b : float
        Mole fraction of light component in the bottoms (0 to 1).
    alpha_values : list of float
        Relative volatility values. Accepts:
            - Single value  [alpha]            — constant alpha assumption
            - Two values    [alpha_top, alpha_bottom]        — geometric mean (Kister default)
            - Three values  [alpha_top, alpha_mid, alpha_bottom] — cubic mean (PE Handbook)
        All values must be >= 1.

    Returns
    -------
    float
        Minimum number of theoretical stages N_min at total reflux.
    """
    if x_b < 0 or x_b > 1:
        raise ValueError("The mole fraction of x in the bottoms must be between 0 and 1.")
    if x_d < 0 or x_d > 1:
        raise ValueError("The mole fraction of x in the distillate must be between 0 and 1.")
    if not alpha_values:
        raise ValueError("At least one alpha value must be provided.")
    if any(a < 1 for a in alpha_values):
        raise ValueError("All relative volatility values must be >= 1.")
    if len(alpha_values) > 3:
        warnings.warn(
            f"{len(alpha_values)} alpha values provided. Standard practice uses 1-3 values. "
            "Verify this is intentional.",
            UserWarning
        )

    # Geometric mean of however many alpha values are provided
    n = len(alpha_values)
    alpha_mean = np.prod(alpha_values) ** (1 / n)

    n_min = np.log((x_d / (1 - x_d)) * ((1 - x_b) / x_b)) / np.log(alpha_mean)

    return n_min

def kirkbride_rectifying_stripping_tray_ratio(z_hk, z_lk,x_blk, x_dhk, bottom_molar_drawrate, dist_molar_drawrate):
    """
    The Kirkbride equation is N{R}/N{S} = [(z{HK}/z{LK})*((x{B,LK}/x{D,HK})**2)*(B/D)]**0.206, where
    the solution is the ratio of the number or rectifying to number of stripping trays.  The total tray count can
    be used to determine the number of rectifying and stripping trays in the column.

    z_hk: float, mole fraction (or rate) of the heavy key in the feed stream
    z_lk: float, mole fraction (or rate) of the light key in the feed stream
    x_blk: float, the mole fraction of the light key in the bottoms draw
    x_dhk: float, the mole fraction of the heavy key in the distillate draw
    bottom_molar_drawrate: float, the total bottoms molar rate
    dist_molar_drawrate: float, the total distillate draw molar rate

    return: {R}/N{S}
    """
    if z_hk <=0 or z_hk >1:
        raise ValueError("The mole fraction of the heavy key in the feed stream must be between 0 and 1.")
    if z_lk <=0 or z_lk >1:
        raise ValueError("The mole fraction of the light key in the feed stream must be between 0 and 1.")
    if x_blk <=0 or x_blk >1:
        raise ValueError("The mole fraction of the light key in the bottom stream must be between 0 and 1.")
    if x_dhk <=0 or x_dhk >1:
        raise ValueError("The mole fraction of the heavy key in the distillate stream must be between 0 and 1.")
    if bottom_molar_drawrate <=0:
        raise ValueError("The bottoms molar draw rate cannot be less than or equal to zero.")
    if dist_molar_drawrate <=0:
        raise ValueError("The distillate molar draw rate cannot be less than or equal to zero.")

    return ((z_hk/z_lk)*((x_blk/x_dhk)**2)*(bottom_molar_drawrate/dist_molar_drawrate))**0.206

def kirkbride_tray_counts(n_min: float, nr_ns_ratio: float) -> dict:
    """
    Given total minimum stages from Fenske and the N_R/N_S ratio from Kirkbride,
    return the number of rectifying and stripping stages.

    Parameters
    ----------
    n_min : float
        Total minimum theoretical stages from Fenske equation.
    nr_ns_ratio : float
        Ratio N_R/N_S from Kirkbride equation.

    Returns
    -------
    dict with keys:
        'N_R' : float — rectifying stages
        'N_S' : float — stripping stages
    """
    if n_min <= 0:
        raise ValueError("Total stage count must be greater than zero.")
    if nr_ns_ratio <= 0:
        raise ValueError("N_R/N_S ratio must be greater than zero.")

    # N_R + N_S = N_min and N_R/N_S = ratio
    # Solving: N_R = N_min * ratio / (1 + ratio)
    n_r = n_min * nr_ns_ratio / (1 + nr_ns_ratio)
    n_s = n_min / (1 + nr_ns_ratio)

    return {"N_R": n_r, "N_S": n_s}

def calc_operating_feedline_intersection(rect_slope, rect_y_int, f, z_f):
    """
    Finds intersection of rectifying operating line and feed line (q-line).
    This point anchors the stripping operating line.

    Parameters
    ----------
    rect_slope : float
        Slope of rectifying operating line.
    rect_y_int : float
        y-intercept of rectifying operating line.
    f : float
        Mole fraction of vapor in feed.
    z_f : float
        Feed composition, mole fraction of light component.

    Returns
    -------
    tuple
        (x_I, y_I) intersection coordinates.
    """
    if f == 0:
        # Vertical q-line: x is fixed at z_f, solve rectifying OL for y
        x_i = z_f
        y_i = rect_slope * x_i + rect_y_int
        return x_i, y_i

    if f == 1:
        # Horizontal q-line: y is fixed at z_f, solve rectifying OL for x
        y_i = z_f
        x_i = (y_i - rect_y_int) / rect_slope
        return x_i, y_i

    # General case: solve simultaneously
    # q-line:        y = feed_slope * x + feed_y_int
    # rectifying OL: y = rect_slope * x + rect_y_int
    feed_slope = calc_feedline_slope(f)
    feed_y_int = z_f / f   # y-intercept of q-line: set x=0 → y = z_f/f

    # Setting equal and solving for x:
    x_i = (rect_y_int - feed_y_int) / (feed_slope - rect_slope)
    y_i = rect_slope * x_i + rect_y_int

    return x_i, y_i

def calc_rectifying_ol_params(x_d: float, reflux_ratio: float) -> tuple:
    """
    Calculates rectifying operating line slope and y-intercept.
    Separated from plotting to allow the coordinator to use these values
    before the line is drawn (specifically for feed line intersection).

    Parameters
    ----------
    x_d : float
        Mole fraction of light component in the distillate (0 to 1).
    reflux_ratio : float
        Reflux ratio R = L/D. Must be >= 0.

    Returns
    -------
    tuple
        (slope, y_intercept)
    """
    slope = calc_rectifying_slope(reflux_ratio)
    y_int = calc_rectifying_y_intercept(x_d, reflux_ratio)
    return slope, y_int

def run_mccabe_thiele_stepper(x_d, x_b, f, z_f, rect_slope, rect_y_int,
                                  strip_slope, x_b_anchor, equil_func,
                                  x_i, condenser_type="total",
                                  reboiler_type="partial",
                                  max_stages=100):
    """
    Executes the McCabe-Thiele stepping algorithm.

    Parameters
    ----------
    x_d : float
        Distillate mole fraction.
    x_b : float
        Bottoms mole fraction.
    f : float
        Mole fraction of vapor in feed.
    z_f : float
        Feed composition.
    rect_slope : float
        Slope of rectifying operating line.
    rect_y_int : float
        y-intercept of rectifying operating line.
    strip_slope : float
        Slope of stripping operating line.
    x_b_anchor : float
        x-intercept anchor of stripping operating line (= x_b).
    equil_func : callable
        Equilibrium curve callable from build_equilibrium_callable.
    x_i : float
        x-coordinate of feed line / operating line intersection.
        Used to determine when to switch operating lines.
    condenser_type : str
        'total' or 'partial'. Partial condenser counts as one stage.
    max_stages : int
        Safety limit to prevent infinite loops.

    Returns
    -------
    dict with keys:
        'stage_coords'       : list of (x, y) tuples defining the staircase
        'total_stages'       : int
        'feed_stage'         : int
        'rectifying_stages'  : int
        'stripping_stages'   : int
    """
    from scipy.optimize import brentq

    def equil_x_from_y(y_target):
        """Invert equilibrium curve: given y, find x."""
        if y_target <= 0:
            return 0.0
        if y_target >= 1:
            return 1.0
        return brentq(lambda x: equil_func(x) - y_target, 0.0, 1.0)

    def rect_ol(x):
        return rect_slope * x + rect_y_int

    def strip_ol(x):
        return strip_slope * (x - x_b_anchor) + x_b_anchor

    # Starting point
    if condenser_type == "total":
        y_current = x_d      # Start on diagonal at (x_d, x_d)
    else:
        y_current = x_d      # Partial condenser handled via stage credit

    stage_coords = [(x_d, x_d)]   # Starting point for staircase plot
    stage_count = 0
    feed_stage = None
    on_stripping = False

    for _ in range(max_stages):
        # --- Horizontal step: find x on equilibrium curve at y_current ---
        x_eq = equil_x_from_y(y_current)
        stage_coords.append((x_eq, y_current))   # horizontal line endpoint
        stage_count += 1

        # --- Check termination: have we passed x_B? ---
        if x_eq <= x_b:
            break

        # --- Determine active operating line ---
        # Switch to stripping OL once x drops below feed intersection
        if not on_stripping and x_eq <= x_i:
            on_stripping = True
            if feed_stage is None:
                feed_stage = stage_count

        # --- Vertical step: drop to active operating line ---
        if on_stripping:
            y_new = strip_ol(x_eq)
        else:
            y_new = rect_ol(x_eq)

        stage_coords.append((x_eq, y_new))   # vertical line endpoint
        y_current = y_new

    else:
        raise ValueError(
            f"Stepper reached {max_stages} stages without achieving x_B separation. "
            "Possible causes:\n"
            "  - Specified reflux ratio is below R_min\n"
            "  - Operating lines cross the equilibrium curve\n"
            "  - Separation targets are infeasible for this system"
        )

    # Feed stage defaults to last stage if switch never triggered
    if feed_stage is None:
        feed_stage = stage_count

    rect_stages = feed_stage - 1
    strip_stages = stage_count - feed_stage + 1

    #Condenser and reboiler type may change the stage count.

    condenser_credit = get_condenser_stage_credit(condenser_type)
    reboiler_credit = get_reboiler_stage_credit(reboiler_type)

    adjusted_total = stage_count - condenser_credit - reboiler_credit
    rect_stages = feed_stage - 1 - condenser_credit
    strip_stages = stage_count - feed_stage + 1 - reboiler_credit

    #print(f"DEBUG stepper: condenser_type={condenser_type}, reboiler_type={reboiler_type}")
   # print(f"DEBUG stepper: stage_count={stage_count}, feed_stage={feed_stage}")
   # print(f"DEBUG stepper: condenser_credit={condenser_credit}, reboiler_credit={reboiler_credit}")
   # print(f"DEBUG stepper: adjusted_total={adjusted_total}, rect_stages={rect_stages}, strip_stages={strip_stages}")

    return {
            "stage_coords"      : stage_coords,
            "total_stages"      : adjusted_total,
            "feed_stage"        : feed_stage,
            "rectifying_stages" : rect_stages,
            "stripping_stages"  : strip_stages,
        }