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

## References

- Kister, H.Z. *Distillation Design*. McGraw-Hill, 1992.
- Smith, J.M., Van Ness, H.C., Abbott, M.M. *Introduction to Chemical 
  Engineering Thermodynamics*. McGraw-Hill.
- Geankoplis, C.J. *Transport Processes and Separation Process Principles*.
