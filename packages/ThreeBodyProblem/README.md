# ThreeBodyProblem

`ThreeBodyProblem` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform

This package is created using the following CC skills:

1. `create-interactive-scientific-application-spec`

User input:
```bash
/create-interactive-scientific-application-spec use the prompt from PROMPT.md. Save result to SPEC.md
```

Output: [SPEC.md](SPEC.md)

2. `implement-interactive-scientific-application-from-spec`

Manual updates:

* bumped version of `diff-grok`


## Source Prompt: Circular Restricted Three-Body Problem

Create spec of an interactive application the Circular Restricted Three-Body Problem (CR3BP).

### Mathematical Model

Circular Restricted Three-Body Problem in the co-rotating (synodic) frame:

```
dx/dt = vx
dy/dt = vy
dvx/dt = 2·vy + x - (1-μ)·(x+μ)/r₁³ - μ·(x-1+μ)/r₂³
dvy/dt = -2·vx + y - (1-μ)·y/r₁³ - μ·y/r₂³
```

Where:
- μ = m₂/(m₁+m₂) — mass parameter
- r₁ = √((x+μ)² + y²) — distance to body 1
- r₂ = √((x-1+μ)² + y²) — distance to body 2

Default parameters: μ = 0.01215 (Earth–Moon system).
Initial conditions: set interactively.

### ODE Solver

- Use the **MRT (Modified Rosenbrock Triple)** method from the **diff-grok** library for numerical integration.
- Step size: adaptive, initial h = 0.001.
- Integration interval: t ∈ [0, T], where T is controlled by a slider (up to 100 dimensionless time units).

### Visualization

#### Scatter Plot — main view (orbit)
- **Scatter plot**: x vs y — orbit of the test body in the rotating frame.
- Two massive bodies shown as large circles: m₁ (blue, "Earth"), m₂ (gray, "Moon").
- Trajectory colored by velocity magnitude: slow segments — cool tones, fast segments — warm tones.
- Display all 5 Lagrange points (L1–L5) as small labeled markers.

#### Line Chart — Jacobi constant (conservation check)
- **Line chart**: C_J(t) — Jacobi integral over time.
- Should remain constant — any drift indicates integration error.

#### Scatter Plot — zero-velocity surfaces
- **Scatter plot** or contour plot: display Zero Velocity Curves for the current C_J value.
- Shade the forbidden regions (Hill regions) with a semi-transparent fill.

### Interactivity

- **Slider μ** (0.001–0.5): presets for "Earth–Moon" (0.01215), "Sun–Jupiter" (0.000953), "Pluto–Charon" (0.1).
- **Click on scatter plot**: sets the initial position (x₀, y₀). Drag to set the direction and magnitude of the initial velocity (shown as an arrow).
- **Sliders vx₀, vy₀** (−2 to 2) — alternative way to set initial velocity.
- **Slider T** (1–100) — integration time.
- Button **"Preset: Halo Orbit around L1"** — automatically loads initial conditions for a characteristic periodic orbit.
- Button **"Preset: Free-Return Trajectory"** — loads initial conditions for a figure-eight-type free-return path.

### Tooltips (mandatory, meaningful)

- **μ (mass parameter)**: "Ratio of the smaller body's mass to the total. For Earth–Moon, μ ≈ 0.012. Increasing μ enlarges the gravitational sphere of influence of the second body — the topology of forbidden regions changes."
- **Lagrange points L1–L3**: "Collinear equilibrium points along the x-axis. Unstable — a small perturbation drives the test body away. L1 is used for solar observatories (SOHO, JWST near Sun–Earth L2)."
- **Lagrange points L4, L5**: "Triangular equilibrium points. Stable when μ < 0.0385… This is where Jupiter's Trojan asteroids reside."
- **Jacobi integral C_J**: "The only integral of motion in CR3BP. Determines which spatial regions are accessible to the test body — higher C_J means more forbidden zones."
- **Zero-velocity surfaces**: "Boundaries of Hill regions: the body cannot cross a forbidden region at a given C_J. As C_J decreases, the 'necks' around L1 and L2 open up, allowing transit between the two massive bodies."
- **Click to launch**: "Click to place the test body. Drag to set the velocity vector. The orbit will be computed instantly."
- **Halo orbit**: "A periodic orbit around L1 or L2 — used by real space missions (e.g., JWST). Requires precise initial conditions."

### Meaningful UI Labels (mandatory)

All headings, axis labels, panel titles, slider labels, and legend entries must be **informative and self-explanatory** — not raw variable names or technical abbreviations, but phrases that let the user immediately understand what they are looking at and why.

Examples of correct naming:
- Scatter plot title: **"Orbit in Rotating Frame (x, y)"** — not "Plot 1", not "Scatter"
- Line chart title: **"Jacobi Constant Drift (integration accuracy)"** — not "C_J(t)", not "Line Chart"
- Scatter plot X-axis: **"x — distance along primary axis"** — not just "x"
- Scatter plot Y-axis: **"y — distance perpendicular to axis"** — not just "y"
- Slider label for μ: **"μ — mass ratio m₂/(m₁+m₂)"** — not "mu", not "param1"
- Slider label for T: **"Integration time (orbital periods)"** — not "T", not "time"
- Slider label for vx₀: **"Initial velocity vx₀ (rotating frame)"** — not "vx0"
- Preset buttons: **"Earth–Moon L1 Halo Orbit"**, **"Free-Return Trajectory"** — not "Preset 1", "Preset 2"
- Forbidden regions panel title: **"Forbidden Regions (Hill surfaces at current C_J)"** — not "ZVC plot"
- Lagrange point legend entries: **"L1 — inner equilibrium"**, **"L4 — leading triangle point"** — not just "L1", "L4"
- Control panel heading: **"Initial Conditions & System Parameters"** — not "Controls", not "Settings"

**Principle**: a first-time user opening the simulation should understand every UI element without consulting any documentation.


## Issues

### 1. App

```bash
logger.ts:105 ReferenceError: Cannot access 'ctrlMu' before initialization
    at getInputs (package.js:13477:23)
    at updateLagrangePanel (package.js:13396:20)
    at threeBodyProblemApp (package.js:13413:5)
    at Module.threeBodyProblemAppEntry (package.js:14184:68)
    at run (eval at <anonymous> (login.dart.js_1.part.js:228734:6), <anonymous>:3:41)
    at Module.Un (js-api.js:1:330776)
    at login.dart.js_1.part.js:228514:11
    at fAp.a (login.dart.js:3886:71)
    at dart.fAp.$2 (login.dart.js:4128:24)
    at dart.eX6.$1 (login.dart.js:4124:31)
```

### 2. Failed tests

A single test failed - fixed automatically
