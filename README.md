# McCabe-Thiele Distillation Diagram Generator

An interactive web application for binary distillation column analysis using 
the McCabe-Thiele graphical method, built with Python and Streamlit.

## What This Tool Does

The **McCabe-Thiele method** is a graphical technique for determining the 
number of theoretical stages required to achieve a specified separation in a 
binary distillation column. It combines three key elements on a single diagram:

- **Equilibrium curve** — the vapor-liquid equilibrium relationship between 
  the more volatile (light) component in the liquid and vapor phases
- **Operating lines** — material balance relationships for the rectifying 
  (above feed) and stripping (below feed) sections of the column
- **Feed line (q-line)** — describes the thermal condition of the feed stream

Theoretical stages are counted by "stepping" between the equilibrium curve 
and the operating lines until the bottoms composition is reached. The feed 
tray is identified at the point where the stepping crosses the feed line.

The **Fenske equation** provides an independent estimate of the minimum 
number of theoretical stages required at total reflux — a useful lower bound 
for comparison against the McCabe-Thiele result at finite reflux.

The **Kirkbride equation** provides an independent estimate of the optimal 
feed tray location based on feed and product compositions and molar draw 
rates. Discrepancies between Kirkbride and McCabe-Thiele feed tray locations 
are normal and reflect different underlying assumptions — engineering judgment 
should be applied when selecting the feed tray for detailed design.

## Features

- Equilibrium curve from constant relative volatility (geometric mean of 
  top/bottom alpha values) or user-supplied x-y VLE data (CSV upload or 
  manual entry)
- Automatic feed condition identification (subcooled liquid, bubble point, 
  partial vaporization, dew point, superheated vapor)
- f or q feed condition input
- Minimum reflux ratio (R_min) calculation and visual reference line
- Rectifying and stripping operating lines
- McCabe-Thiele stepping algorithm with automatic feed tray identification
- Total and partial condenser/reboiler stage credit handling
- Fenske minimum stage calculation (constant alpha path only)
- Kirkbride feed tray location estimate

## Installation
```bash
pip install streamlit matplotlib scipy numpy pandas
streamlit run app.py
```

## File Structure
```
dist_formulas.py   — All thermodynamic and distillation calculations
dist_plots.py      — Matplotlib plot functions and McCabe-Thiele coordinator
app.py             — Streamlit user interface
test_dist_formulas.py — Unit tests (see TODO below)
dist_visuals.py    — Placeholder for future column schematic visualization
```

## Assumptions and Limitations

- Binary systems only
- Constant relative volatility assumed for R_min and Fenske calculations
- Real VLE data supported for equilibrium curve and stepping, but Fenske 
  and R_min are unavailable on this path
- Murphree tray efficiency not yet implemented — all stages are theoretical
- Heat duty calculations not in scope for this version

## TODO

- **Verification and validation cases** — locate and implement well-vetted 
  textbook reference problems (Geankoplis, Seader & Henley, or Kister) with 
  known stage counts for use in `test_dist_formulas.py`
- Unit tests for all functions in `dist_formulas.py`

## Potential Future Functionality

- **Column schematic visualization** — dynamic SVG diagram showing feed, 
  condenser, reboiler, and stream compositions (`dist_visuals.py`)
- **Thermodynamic feed condition calculation** — compute f from feed 
  temperature, pressure, Cp, and enthalpy of vaporization
- **Heat duty estimation** — condenser and reboiler duties from molar flow 
  rates and latent heat
- **Multicomponent extension** — Hengstebeck diagram for pseudo-binary 
  representation of multicomponent systems
- **Export** — download diagram as PNG and results as CSV

## References

- Kister, H.Z. *Distillation Design*. McGraw-Hill, 1992.
- Smith, J.M., Van Ness, H.C., Abbott, M.M. *Introduction to Chemical 
  Engineering Thermodynamics*. McGraw-Hill.
- Geankoplis, C.J. *Transport Processes and Separation Process Principles*.
