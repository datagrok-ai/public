[//]: # (doc max 10 lines: public/help/visualize/viewers/grid#summary-columns: confidence interval)

# Confidence Interval Bar

## Model

The renderer supports multiple input modes, specified in the 'type' column setting.
The type dictates which settings are shown and taken into account.

- "Value + Min Max": estimate, lower bound, upper bound
- "Value + Margin": estimate + margin of error (symmetric CI)

Center mark: dot, diamond, or vertical line for the point estimate
Whiskers: horizontal line with optional end caps (serifs) for CI bounds
Optional fill: semi-transparent band between lower and upper bounds

## Scale & Axis

Global scale (default): all cells in the column share the same min/max → values are visually comparable across rows
Per-row scale: each cell scales independently (useful when magnitudes vary wildly)
Custom fixed scale: user-defined min/max (e.g., 0–100 for percentages)
Symmetric around zero: for difference/change metrics
Reference line: optional vertical line at a meaningful value (zero, target, control mean)
Log scale option for data spanning orders of magnitude (common in pharma — IC50s, concentrations)

Color options: 
- single user-chosen color for all cells (default option)
- a column that is used for color-coding

## Interactivity

Tooltip on hover: show exact values — point estimate, lower, upper, CI width, confidence level, and N if available

## Edge Cases & Data Quality

Missing bounds: render point estimate only (as a dot with no whiskers)
Missing point estimate: render range only (bar without center mark)
All missing: render empty cell or a "no data" indicator
Inverted bounds (lower > upper): render with warning indicator (e.g., dashed line + orange)
Clipped whiskers: if bounds extend beyond the scale, show an arrow/chevron at the edge to indicate truncation
Infinite bounds: one-sided CI (e.g., lower bound only) — render open-ended whisker with arrow
Negative values: handle gracefully when scale includes negative range