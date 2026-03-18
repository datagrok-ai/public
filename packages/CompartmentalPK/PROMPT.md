# Prompt: Multi-Compartment Pharmacokinetic Model with Brute-Force Parameter Optimization

Create an interactive application of a two- and three-compartment PK model with built-in brute-force parameter fitting.

## Mathematical Model

### Two-compartment model (toggle):
```
dA_c/dt  = -k₁₀·A_c - k₁₂·A_c + k₂₁·A_p + R_in(t)
dA_p/dt  =  k₁₂·A_c - k₂₁·A_p
```

### Three-compartment model (toggle):
```
dA_c/dt  = -k₁₀·A_c - k₁₂·A_c - k₁₃·A_c + k₂₁·A_p + k₃₁·A_d + R_in(t)
dA_p/dt  =  k₁₂·A_c - k₂₁·A_p
dA_d/dt  =  k₁₃·A_c - k₃₁·A_d
```

Where:
- A_c — drug amount in the central compartment
- A_p — peripheral compartment (rapid exchange)
- A_d — deep compartment (slow exchange, three-compartment only)
- R_in(t) — input function (bolus, infusion, or oral absorption)
- Plasma concentration: C(t) = A_c(t) / V_c

Micro-constants derived from physiological parameters:
- k₁₀ = CL / V_c
- k₁₂ = Q / V_c,  k₂₁ = Q / V_p
- k₁₃ = Q₂ / V_c,  k₃₁ = Q₂ / V_d  (three-compartment only)

Default parameters (two-compartment):
- CL = 5 L/h, V_c = 20 L, V_p = 40 L, Q = 10 L/h

## Input Modes (R_in)

1. **IV Bolus**: instantaneous dose D added to A_c at t = 0 (and every τ hours for multiple dosing).
2. **IV Infusion**: R_in = D / T_inf during infusion, 0 afterwards.
3. **Oral (first-order absorption)**: add equation dA_gut/dt = −k_a·A_gut, and R_in = k_a·A_gut for the central compartment.

## ODE Solver

- Use the **MRT (Modified Rosenbrock Triple)** method from the **diff-grok** library for numerical integration.
- Step size: adaptive, initial h = 0.01.
- Integration interval: determined by dosing regimen (up to 120 h for multiple-dose schedules).
- **Critical for optimization**: the solver is called hundreds/thousands of times during brute-force search — MRT's efficiency with adaptive stepping is key here.

## Visualization

### Line Chart — main (concentration–time profile)
- **Line chart**: C(t) = A_c(t) / V_c — plasma concentration over time.
- Toggle between linear and semi-log Y-axis.
- On linear scale, show:
  - C_max and t_max marked with a labeled point.
  - Horizontal dashed lines for MEC (minimum effective concentration) and MTC (minimum toxic concentration).
  - Area between MEC and MTC shaded as the "therapeutic window".
- On semi-log scale, the distribution (α) and elimination (β) phases appear as straight lines with different slopes.

### Line Chart — multiple dosing
- When "repeated doses" is selected: show accumulation toward steady state.
- Label C_ss_max and C_ss_min.
- Annotate the point where ~90% of steady state is reached.

### Line Chart — compartment amounts
- Drug amount in each compartment over time: A_c(t), A_p(t), A_d(t) as separate lines.

### Animated compartment diagram
- SVG schematic: rectangles labeled "Central", "Peripheral", "Deep" with directional arrows and rate constant labels.
- Arrow thickness proportional to current mass flow at time t (animated).

## Brute-Force Parameter Optimization

### Concept
The system generates synthetic observed concentration–time data by running the model with current parameters and adding realistic noise. Then it performs a **grid search** over selected PK parameters, solving the ODE at every grid point via MRT, and finds the parameter combination that minimizes the objective function. This demonstrates how brute-force optimization recovers the "true" parameters from noisy observations.

### Synthetic Data Generation
- **"Generate Noisy Observations"** button: runs the model with current parameters, samples at predefined clinical time points (0.5, 1, 2, 4, 8, 12, 24 h), and adds log-normal noise (CV adjustable via slider, default 15%).
- **Noise CV slider** (5–40%) — controls the spread of synthetic noise. Higher CV = harder optimization problem.
- **Sampling schedule**: toggle between "Standard PK" (0.5, 1, 2, 4, 8, 12, 24 h), "Rich" (adds 0.25, 0.75, 1.5, 3, 6, 10, 16, 20 h), and "Sparse" (1, 4, 12, 24 h) — to show how data density affects parameter identifiability.
- Observed data shown as **scatter plot** (circles) overlaid on the model line chart.

### Optimization Setup
- **Parameter selection**: checkboxes to choose which parameters to optimize. At least 1, at most 3 simultaneously (for visualization and speed).
  - Available: CL, V_c, V_p, Q, k_a (if oral), F (if oral).
- **Search range**: for each selected parameter, min/max sliders defining the search grid.
  - Defaults: ±50% around the current parameter value.
- **Grid resolution**: slider for number of grid points per parameter axis (5–50, default 20).
  - Show estimated number of ODE solves: N_solves = grid_points^(number_of_parameters).
  - Show warning if N_solves > 10,000: "This may take a few seconds."

### Objective Function
- **Weighted Sum of Squared Residuals (WSSR)**:
  `WSSR = Σ wᵢ · (C_obs,i − C_pred,i)²`
- Weight options (toggle):
  - Uniform: wᵢ = 1
  - 1/C_obs²: proportional weighting (standard for PK)
  - 1/C_pred²: iteratively reweighted

### Optimization Execution & Results
- **"Run Optimization"** button with progress bar showing % of grid evaluated.
- On completion, display:
  - **Best-fit parameters** with values and units.
  - **Best-fit curve** overlaid on observed data (distinct style from the manual model curve).
  - **Goodness of fit**: WSSR value, R², AIC.

### Optimization Landscape Visualization (key differentiator!)
- **If 1 parameter optimized**: **Line chart** — objective function (WSSR) vs parameter value. Minimum marked with a vertical line. Shows whether the minimum is sharp (well-identified) or flat (poorly identified).
- **If 2 parameters optimized**: **Heatmap / contour plot** — WSSR as color over the 2D parameter grid. Minimum marked with a crosshair. Contour lines show confidence regions. Correlation between parameters visible as elongated valleys.
- **If 3 parameters optimized**: three **2D heatmaps** (pairwise projections, marginalizing over the third parameter). Each shows the minimum WSSR across the third dimension.

### Compare: Slider Parameters vs Optimized Fit
- Toggle to overlay both curves: the current slider-set parameters and the optimized best-fit.
- Residual plot: **scatter plot** of (C_obs − C_pred) vs time for both parameter sets.

## Interactivity

- **Toggle**: two-compartment / three-compartment model.
- **Toggle input mode**: IV Bolus / IV Infusion / Oral.
- **PK parameter sliders**: CL (0.1–50 L/h), V_c (1–100 L), V_p (1–200 L), Q (0.1–50 L/h).
- **For three-compartment**: V_d (1–200 L), Q₂ (0.01–5 L/h).
- **For oral**: k_a (0.1–5 h⁻¹), F (bioavailability 0–1).
- **Dosing**: Dose (1–1000 mg), interval τ (1–48 h), number of doses (1–20), infusion duration (0.5–24 h).
- **MEC and MTC sliders**: set therapeutic window boundaries.
- **Disease presets**: buttons "Vancomycin", "Theophylline", "Metformin" — load real-world PK parameters.

## Tooltips (mandatory, meaningful)

- **CL (clearance)**: "Volume of plasma completely cleared of drug per unit time. Determines the slope of the terminal phase on a semi-log plot. Half-life t₁/₂ = 0.693 · V_ss / CL."
- **V_c (central volume)**: "Apparent volume of distribution in blood/plasma. Determines the initial concentration after an IV bolus: C₀ = Dose / V_c."
- **Q (intercompartmental clearance)**: "Rate of drug exchange between central and peripheral compartments. Large Q → fast distribution → short α-phase."
- **Therapeutic window**: "The range between MEC (below = no effect) and MTC (above = toxicity). The goal of dosing is to keep the concentration within this corridor."
- **Semi-log plot**: "On a log scale, exponential processes appear as straight lines. Two distinct slopes = two phases: fast (α, distribution) and slow (β, elimination)."
- **Steady state**: "With repeated dosing, concentration reaches steady state after ~5 half-lives. C_ss depends only on CL and dose/interval: C_ss_avg = F · Dose / (CL · τ)."
- **k_a (absorption rate)**: "Rate of absorption from the GI tract. When k_a >> k₁₀ (no flip-flop kinetics), the peak comes quickly. When k_a ≈ k₁₀, the peak is blunted and delayed."
- **Bioavailability F**: "Fraction of the dose reaching systemic circulation after oral administration. F = 1.0 for IV; typically 0.1–0.9 for oral due to first-pass metabolism."
- **Compartment diagram**: "Arrow thickness is proportional to current mass flow. Right after a bolus — strong flow from central to peripheral. Later — reverse flow (redistribution)."
- **Generate synthetic data**: "Runs the model at current parameters, samples at clinical time points, and adds log-normal noise. Try recovering the true parameters with the optimizer — higher CV makes it harder!"
- **Grid resolution**: "Number of points per parameter axis. Total ODE solves = grid_points ^ number_of_parameters. For 2 parameters at 30 grid points = 900 solves — usually takes < 2 seconds with MRT."
- **WSSR objective**: "Weighted Sum of Squared Residuals — the quantity being minimized. 1/C² weighting is standard in PK because measurement error is typically proportional to concentration."
- **Optimization landscape (1D)**: "If the curve has a sharp V-shape, the parameter is well-identified by the data. A flat or multi-modal curve means the data cannot uniquely determine this parameter."
- **Optimization landscape (2D heatmap)**: "Elongated valleys indicate correlation between parameters — many parameter combinations fit equally well. Round contours mean parameters are independently identifiable."
- **Residual plot**: "Residuals should be randomly scattered around zero. Systematic patterns (e.g., all positive at early times) suggest the model structure is wrong, not just the parameter values."

## Meaningful UI Labels (mandatory)

All headings, axis labels, panel titles, slider labels, and legend entries must be **informative and self-explanatory**.

Examples of correct naming:
- Main chart title: **"Plasma Concentration–Time Profile"** — not "C(t)", not "PK Plot"
- Main chart Y-axis: **"Concentration (mg/L)"** — not "C", not "conc"
- Main chart X-axis: **"Time after first dose (hours)"** — not "t", not "hours"
- Semi-log toggle: **"Switch to semi-log scale (reveals α and β phases)"** — not "Log Y"
- Therapeutic window legend: **"Therapeutic window (MEC to MTC)"** — not "window"
- MEC line: **"MEC — minimum effective concentration"** — not "MEC"
- Slider CL: **"CL — clearance (L/h)"** — not "CL", not "param1"
- Slider V_c: **"V_c — central volume of distribution (L)"** — not "Vc"
- Dosing interval: **"τ — dosing interval (hours)"** — not "tau"
- Compartment diagram title: **"Drug Distribution Between Compartments"** — not "Diagram"
- Arrows in diagram: **"CL = 5 L/h (elimination)"**, **"Q = 10 L/h (distribution)"** — not just "k₁₀", "k₁₂"
- Preset buttons: **"Vancomycin (CL=4.6, V_c=28, t₁/₂≈6h)"**, **"Theophylline (CL=3.5, V_c=30, t₁/₂≈8h)"** — not "Preset 1"
- Optimization panel title: **"Brute-Force Parameter Estimation"** — not "Optimize", not "Fit"
- Parameter checkboxes: **"Optimize CL (current: 5 L/h, search: 2.5–7.5 L/h)"** — not "CL ☑"
- Grid resolution slider: **"Grid density per axis (total solves: 400)"** — not "N", not "resolution"
- Run button: **"Run Grid Search (estimated: 400 ODE solves)"** — not "Run", not "Go"
- Results panel: **"Best-Fit Parameters"** — not "Results"
- Heatmap title: **"Optimization Landscape — WSSR over CL × V_c grid"** — not "Heatmap"
- Heatmap colorbar: **"WSSR (lower = better fit)"** — not "value"
- Residual plot title: **"Residuals — Observed minus Predicted"** — not "Residuals"
- Residual Y-axis: **"C_obs − C_pred (mg/L)"** — not "residual"
- Compare toggle: **"Overlay: current sliders vs optimized fit"** — not "Compare"
- Synthetic data button: **"Generate Noisy Observations (CV = 15%)"** — not "Generate Data"
- Noise slider: **"Observation noise CV (coefficient of variation, %)"** — not "noise", not "CV"
- Sampling toggle: **"Sampling schedule: Standard PK / Rich / Sparse"** — not "Schedule"

**Principle**: a first-time user — even one without PK training — should understand every UI element without consulting any documentation. A PK-trained user should find the labels precise and unambiguous.
