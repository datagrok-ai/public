# Application Specification: Golf Ball Flight Simulator

<!-- TOC for partial reading -->
<!-- Section 1: General Architecture — lines 15–330 -->
<!-- Section 2: Main View — lines 331–339 -->
<!-- Section 3: Controls — lines 340–430 -->
<!-- Section 4: Display Elements — lines 431–510 -->
<!-- Section 5: Layout — lines 511–575 -->
<!-- Section 6: User Feedback — lines 576–625 -->
<!-- Section 7: Validation — lines 626–670 -->
<!-- Section 8: Pipeline — lines 671–770 -->
<!-- Section 9: Reactivity — lines 771–790 -->
<!-- Sections 10–15: Data, Errors, Resources, Closure, UX, Testing — lines 791–1000 -->

## 1. General Architecture

### 1.0. General Information

| Field | Value |
|---|---|
| Application name | Golf Ball Flight Simulator |
| Package | GolfBall |
| Entry function | `golfBallApp()` |
| Brief description | Interactive 2D ballistic simulation of a golf ball with air drag, trajectory visualization, height optimization, and sensitivity analysis. |
| Main view type | `DG.TableView` |

### 1.1. Core

#### Task List

| Task ID | Name | Pipeline type | Trigger | Synchronicity | Execution environment | Parallelization |
|---|---|---|---|---|---|---|
| `task_primary` | Golf ball flight simulation | Primary (reactive) | Any input change | Sync | Main thread | No |
| `task_maximize_height` | Maximize height (grid search) | Secondary (on demand) | Button `btn_maximize` | Async | WebWorkers (parallel) | Yes — 1000 grid points distributed across worker pool |
| `task_sensitivity` | Sensitivity analysis | Secondary (on demand) | Button `btn_sensitivity` | Async | WebWorkers (parallel) | Yes — 6 parameter sweeps run concurrently |

#### Dependencies Between Tasks

```
task_primary          — independent
task_maximize_height  — does not depend on task_primary results;
                        after completion, triggers a single run of task_primary
task_sensitivity      — independent; opens results in a new TableView
```

#### Task: `task_primary`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Golf ball flight simulation |
| Description | Numerical solution of the 2D ballistic ODE system with air drag using MRT, returning trajectory x(t), y(t), vx(t), vy(t), v(t), and summary statistics |
| Dependency on other tasks | No |

**Computation Formulas and Model**

**Level 1 — required minimum:**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| `x` | Horizontal position | m | `x ≥ 0` |
| `y` | Vertical position (height) | m | `y ≥ 0` during flight |
| `vx` | Horizontal velocity | m/s | any |
| `vy` | Vertical velocity | m/s | any |
| `v` | Total speed | m/s | `v ≥ 0` |
| `v0` | Initial speed | m/s | `v0 > 0` |
| `theta` | Launch angle | degrees | `0 < theta < 90` |
| `m` | Ball mass | kg | `m > 0` |
| `d` | Ball diameter | m | `d > 0` |
| `Cd` | Drag coefficient | dimensionless | `Cd > 0` |
| `rho` | Air density | kg/m³ | `rho > 0` |
| `g` | Gravitational acceleration | m/s² | `g > 0` |

Relationships:

```
A = π·(d/2)²                          (cross-section area)
v = √(vx² + vy²)                      (total speed)
drag = Cd·ρ·A·v / (2·m)               (drag factor)

ODE system:
  dx/dt  = vx
  dy/dt  = vy
  dvx/dt = −drag · vx
  dvy/dt = −g − drag · vy

Initial conditions:
  x(0)  = 0
  y(0)  = 0
  vx(0) = v0 · cos(θ·π/180)
  vy(0) = v0 · sin(θ·π/180)

Flight ends when y(t) ≤ 0 for t > 0 (ball hits the ground).
```

Output properties (invariants):

| Property | Description |
|---|---|
| `x(t) ≥ 0` for all `t` | Horizontal position is non-negative |
| `y(0) = 0` | Ball starts on the ground |
| `vx(0) = v0·cos(θ)` | Initial horizontal velocity preserved |
| `vy(0) = v0·sin(θ)` | Initial vertical velocity preserved |
| `v(t) = √(vx(t)² + vy(t)²)` | Total speed is consistent |
| `maxHeight ≥ 0` | Maximum height is non-negative |
| `flightDistance ≥ 0` | Total horizontal distance is non-negative |

Reference examples:

| # | Inputs | Expected output | Source |
|---|---|---|---|
| 1 | `v0=70, θ=45, m=0.0459, d=0.0427, Cd=0.25, ρ=1.225, g=9.81` | `vx(0)≈49.497, vy(0)≈49.497` | Manual: `70·cos(45°)`, `70·sin(45°)` |
| 2 | `v0=70, θ=45, Cd=0` (no drag) | `maxHeight ≈ v0²·sin²(θ)/(2g) ≈ 124.8 m`, `distance ≈ v0²·sin(2θ)/g ≈ 499.5 m` | Analytical projectile formula |
| 3 | At `t=0`, any params | `dvx/dt = −drag·v0·cos(θ)`, `dvy/dt = −g − drag·v0·sin(θ)` | Manual calculation of RHS |

**Level 2 — full formalization:**

- Complete mathematical formulation: 2D ballistic ODE with quadratic drag, IVP from ground level.
- Analytical properties: without drag, parabolic trajectory; with drag, range reduced, optimal angle < 45°.
- Numerical method: MRT (adaptive multi-rate), suitable for both stiff and non-stiff ODEs.

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `v0` | `number` | m/s | `> 0` | Initial speed |
| `theta` | `number` | degrees | `(0, 90)` | Launch angle |
| `m` | `number` | kg | `> 0` | Ball mass |
| `d` | `number` | m | `> 0` | Ball diameter |
| `Cd` | `number` | dimensionless | `> 0` | Drag coefficient |
| `rho` | `number` | kg/m³ | `> 0` | Air density |
| `g` | `number` | m/s² | `> 0` | Gravitational acceleration |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `t` | `Float64Array` | Time values |
| `x` | `Float64Array` | Horizontal position |
| `y` | `Float64Array` | Vertical position (height) |
| `v` | `Float64Array` | Total speed |
| `maxHeight` | `number` | Maximum height achieved |
| `maxHeightTime` | `number` | Time of maximum height |
| `flightDistance` | `number` | Horizontal distance at landing |
| `flightTime` | `number` | Total flight time |
| `landingSpeed` | `number` | Speed at landing |
| `landingAngle` | `number` | Angle at landing (degrees, positive down) |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | External library | `diff-grok` v1.2.0, function `mrt(task: ODEs)`. MRT is an adaptive multi-rate method suitable for both stiff and non-stiff ODEs. | [diff-grok README](https://github.com/datagrok-ai/diff-grok) |
| 2 | Custom method | Truncate solution arrays at ground impact (`y ≤ 0`) via linear interpolation for exact landing point | — |
| 3 | Custom method | Compute summary statistics: `maxHeight`, `maxHeightTime`, `flightDistance`, `flightTime`, `landingSpeed`, `landingAngle` | — |

**Library call:**

```typescript
import { ODEs, mrt } from 'diff-grok';

const thetaRad = theta * Math.PI / 180;
const A = Math.PI * (d / 2) ** 2;

const task: ODEs = {
  name: 'GolfBallFlight',
  arg: { name: 't', start: 0, finish: 30, step: 0.01 },
  initial: [0, 0, v0 * Math.cos(thetaRad), v0 * Math.sin(thetaRad)],
  func: (_t, y, out) => {
    const vx = y[2], vy = y[3];
    const speed = Math.sqrt(vx * vx + vy * vy);
    const drag = Cd * rho * A * speed / (2 * m);
    out[0] = vx;           // dx/dt
    out[1] = vy;           // dy/dt
    out[2] = -drag * vx;   // dvx/dt
    out[3] = -g - drag * vy; // dvy/dt
  },
  tolerance: 1e-6,
  solutionColNames: ['x(t)', 'y(t)', 'vx(t)', 'vy(t)'],
};
const solution = mrt(task);
```

**ODE right-hand side reference examples:**

| # | v0 | θ | m | d | Cd | ρ | g | vx | vy | Expected dvx/dt | Expected dvy/dt | Derivation |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 70 | 45 | 0.0459 | 0.0427 | 0.25 | 1.225 | 9.81 | 49.497 | 49.497 | −15.17 | −25.0 | A=0.001432, v=69.99, drag=Cd·ρ·A·v/(2m)=0.3064, dvx=−0.3064·49.497=−15.17, dvy=−9.81−0.3064·49.497=−25.0 |
| 2 | 50 | 30 | 0.0459 | 0.0427 | 0.25 | 1.225 | 9.81 | 43.301 | 25.0 | −10.59 | −15.92 | v=50, drag=0.2442, dvx=−0.2442·43.301=−10.57, dvy=−9.81−0.2442·25=−15.92 |

**Execution environment constraint:** main thread (synchronous, < 100 ms).

---

#### Task: `task_maximize_height`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Maximize Height (grid search) |
| Description | Brute-force grid search over a user-selected parameter (θ, v0, or Cd) with 1000 steps. For each value, solves the ODE and finds the peak height. Returns the value that maximizes peak altitude. |
| Dependency on other tasks | Uses current parameter values from primary pipeline controls |

**Task input parameters:**

| Parameter | Type | Description |
|---|---|---|
| `paramName` | `OptimizeParam` | Selected parameter to optimize (`theta`, `v0`, or `Cd`) |
| `paramMin, paramMax` | `number` | Range of the selected parameter |
| `steps` | `number` | Number of grid steps (1000) |
| `baseParams` | `GolfBallParams` | Current values of all parameters (snapshot) |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `optimalValue` | `number` | Value of the parameter that maximizes height |
| `maxHeight` | `number` | Maximum height achieved at the optimal value |

**Computation implementation:**

| Step | Implementation method | Details |
|---|---|---|
| 1 | Custom method | Generate 1000 values of the selected parameter via linspace |
| 2 | External library in workers | `diff-grok` v1.2.0, `mrt(task)` — for each grid point in a WebWorker |
| 3 | Custom method | Find grid point with maximum peak height |
| 4 | Datagrok API | Write optimal value to corresponding control via batch update |

**Parallelization strategy:**

```
Number of workers = Math.max(1, navigator.hardwareConcurrency - 2)

1000 grid points
  → distributed round-robin to worker pool
  → each worker receives: { paramName, values[], baseParams }
  → each worker returns: { value, maxHeight }[] for each point
  → after all complete: find max(maxHeight) → optimal value
```

---

#### Task: `task_sensitivity`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Sensitivity Analysis |
| Description | For each of 6 model parameters independently (v0, θ, Cd, ρ, m, d), sweeps across its full range (100 steps), computes total flight distance. Opens results in a new TableView with line charts. |
| Dependency on other tasks | Uses current parameter values from primary pipeline controls |

**Task input parameters:**

| Parameter | Type | Description |
|---|---|---|
| `baseParams` | `GolfBallParams` | Current values of all parameters (snapshot) |
| `sweepParams` | `SweepParam[]` | Array of 6 parameter specs: `{ name, min, max, steps }` |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `results` | `SweepResult[]` | Array of 6 results, each containing `{ paramName, paramValues[], distances[] }` |

**Computation implementation:**

| Step | Implementation method | Details |
|---|---|---|
| 1 | Custom method in workers | For each of 6 parameters: generate 100 values, solve ODE for each, record flight distance |
| 2 | Datagrok API | Create DataFrame with columns `[param1_values, param1_distance, param2_values, ...]`, open new TableView, add line charts |

**Parallelization strategy:**

```
6 parameter sweeps × 100 steps each = 600 total ODE solves
  → 6 workers (one per parameter sweep)
  → each worker receives: { paramName, paramMin, paramMax, steps, baseParams }
  → each worker returns: { paramName, paramValues[], distances[] }
  → results assembled into DataFrame
```

---

### 1.2. Ports

#### Task ports: `task_primary`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `GolfBallParams` | `{ v0, theta, m, d, Cd, rho, g }` | 7 validated numeric parameters |
| Output | `GolfBallSolution` | `{ t, x, y, v, maxHeight, maxHeightTime, flightDistance, flightTime, landingSpeed, landingAngle }` | Solution arrays + derived stats |

#### Task ports: `task_maximize_height`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `OptimizeTask` | `{ paramName, paramMin, paramMax, steps, baseParams }` | Parameter to optimize + base values |
| Output | `OptimizeResult` | `{ value, maxHeight, error? }` | Optimal parameter value |

#### Task ports: `task_sensitivity`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `SensitivityTask` | `{ paramName, paramMin, paramMax, steps, baseParams }` | Parameter to sweep + base values |
| Output | `SensitivityResult` | `{ paramName, paramValues[], distances[], error? }` | Distance vs parameter |

#### Application-level ports

| Port | Used | Description |
|---|---|---|
| Progress | Yes | `DG.TaskBarProgressIndicator` for `task_maximize_height` and `task_sensitivity` (indeterminate, cancelable) |
| Cancellation | Yes | `pi.onCanceled` → terminates workers, resolves promises |
| Data | No | No external data loading |

### 1.3. Adapters

| Adapter | Used | Implementation |
|---|---|---|
| UI adapter | Yes | Datagrok inputs (`ui.input.float` × 7), `ui.bigButton`, `ui.iconFA`, `ui.icons.help`, `ui.input.choice` |
| Display adapter | Yes | Datagrok viewers (line chart × 2, scatter plot), custom `comp_stats` |
| Worker adapter | Yes | Web Worker wrappers for `task_maximize_height` (`workers/optimize-worker.ts`) and `task_sensitivity` (`workers/sensitivity-worker.ts`) |
| Progress adapter | Yes | `DG.TaskBarProgressIndicator` (indeterminate, cancelable) |
| Data adapter | No | No external data loading |

### 1.4. Coordinator

The coordinator is the `golfBallApp()` function in `app.ts`:

- **Input listening:** each `ui.input.float` calls `debouncedRun()` from `onValueChanged`.
- **Reactivity management:** all 7 controls trigger the same primary pipeline; no cascading dependencies between controls.
- **Validation trigger:** `runPrimary()` calls `validate(inputs)` before computation; validators are also attached to each input via `addValidator()` for inline hints.
- **Control state during computations:** `computationsBlocked` flag prevents `runPrimary()` during batch updates.
- **Computation blocking:** `computationsBlocked = true` before batch writes, `= false` after, then single `runPrimary()`.
- **Results to display:** `updateDataFrame()` creates new `DG.DataFrame` and assigns to `view.dataFrame` and all viewers; `updateStatsPanel()` updates labels.
- **Resource lifecycle:** `subs[]` collects subscriptions; `activeWorkers[]` tracks workers; both cleaned up in `onViewRemoved` handler.

### 1.5. Independence Principle

Confirmed. Input controls and their reactivity (debounce, validation hints) function independently of the computational core. The core receives a ready `GolfBallParams` object; it has no knowledge of UI elements. The `validate()` function is a pure function in `core.ts`.

---

## 2. Main View

| Field | Value |
|---|---|
| View type | `DG.TableView` |
| Description | Table view with docked trajectory chart, time series chart, grid, and left-panel form with controls and stats |

---

## 3. Controls (Inputs)

### 3.1. Primary Pipeline Controls

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_v0` | Initial speed v₀ | `ui.input.float` | `number` | `70` | `10` | `120` | `0.0` | No | Initial ball speed off the clubface. A typical driver produces 50–75 m/s. | Launch |
| `ctrl_theta` | Launch angle θ | `ui.input.float` | `number` | `45` | `1` | `89` | `0.0` | No | Launch angle relative to the ground. Optimal for distance is usually 10–15° with drag, ~45° without. | Launch |
| `ctrl_m` | Ball mass m | `ui.input.float` | `number` | `0.0459` | `0.01` | `0.1` | `0.0000` | No | Mass of the golf ball. Regulation balls weigh no more than 45.93 g (0.0459 kg). | Ball |
| `ctrl_d` | Ball diameter d | `ui.input.float` | `number` | `0.0427` | `0.02` | `0.08` | `0.0000` | No | Diameter of the ball. Regulation minimum is 42.67 mm (0.0427 m). | Ball |
| `ctrl_Cd` | Drag coefficient Cd | `ui.input.float` | `number` | `0.25` | `0.05` | `0.6` | `0.00` | No | Aerodynamic drag coefficient. A smooth ball is ~0.5, a dimpled golf ball ~0.25. | Environment |
| `ctrl_rho` | Air density ρ | `ui.input.float` | `number` | `1.225` | `0.5` | `1.5` | `0.000` | No | Air density. Sea level ≈ 1.225 kg/m³. Higher altitude or warmer air → lower density → longer flight. | Environment |
| `ctrl_g` | Gravity g | `ui.input.float` | `number` | `9.81` | `1.0` | `25.0` | `0.00` | No | Gravitational acceleration. Earth ≈ 9.81 m/s². Moon ≈ 1.62, Mars ≈ 3.72, Jupiter ≈ 24.79. | Environment |

### 3.2. Secondary Task Triggers

| ID | Label / icon | Launches task | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_maximize` | `ui.bigButton('Maximize Height')` | `task_maximize_height` (via dialog) | Find the value of the selected parameter that produces the highest ball apex. | Always (disabled during optimization and sensitivity) |
| `btn_sensitivity` | `ui.iconFA('chart-line')` on ribbon | `task_sensitivity` | Show how each parameter independently affects the total flight distance. | Always |

### 3.3. Secondary Task Controls

#### Task controls: `task_maximize_height`

UI type: Datagrok dialog (`ui.dialog('Maximize Height')`).

| ID | Label | Control type | Data type | Default | Nullable | Tooltip text |
|---|---|---|---|---|---|---|
| `dlg_optimize_param` | Parameter to optimize | `ui.input.choice` | `OptimizeParam` | `'theta'` | No | Choose which parameter to vary in the search for maximum ball height. |

Dialog buttons:
- **Cancel** — closes dialog
- **Run** — closes dialog, launches `task_maximize_height`. Tooltip: "Run brute-force grid search to find the parameter value that maximizes peak altitude."

### 3.4. Other Buttons and Actions

| ID | Label / icon | Action | Tooltip text | Ribbon group | Availability condition |
|---|---|---|---|---|---|
| `btn_reset` | `ui.iconFA('undo')` | Reset all controls to default values | Reset all parameters to default values | Actions | Always |
| `btn_help` | `ui.icons.help` | Open Wikipedia article on projectile motion in a new tab | Learn more about projectile motion with drag. Opens in a separate page. | Help | Always |

### 3.5. Custom UI Components

| Component ID | Brief description | Role |
|---|---|---|
| `comp_stats` | Display panel showing flight statistics. Built with `ui.divV`/`ui.label` components; values updated via `textContent`. Styled via `.golf-ball-app-stats-panel`. | Display |

---

## 4. Result Display Elements

### 4.1. Primary Pipeline Display Elements

| ID | Type | Associated output data | Docking location |
|---|---|---|---|
| `view_trajectory` | Datagrok viewer `scatter plot` | `x`, `y` | Right area, top |
| `view_timeseries` | Datagrok viewer `line chart` | `Time`, `x`, `y`, `v` | Right area, bottom |
| `view_table` | `DG.TableView` grid | `Time`, `x`, `y`, `v` | Center (default) |
| `comp_stats` | Custom HTMLElement | Flight statistics | Left panel, below controls |

**Details for `view_trajectory`:**

| Property | Value |
|---|---|
| X axis | `x` |
| Y axis | `y` |
| Title | Trajectory |

**Details for `view_timeseries`:**

| Property | Value |
|---|---|
| X axis | `Time` |
| Y axis | `x`, `y`, `v` (three series) |
| Title | Time Series |

### 4.2. Secondary Task Display Elements

`task_maximize_height` writes optimal value back to primary pipeline control (batch update).

#### Task display: `task_sensitivity`

Results are displayed in a **new TableView** (separate from main view).

| ID | Type | Associated output data | Placement |
|---|---|---|---|
| `sens_view` | `DG.TableView` | DataFrame with 12 columns (6 param value columns + 6 distance columns) | New view |
| `sens_charts` | Datagrok viewers `line chart` × 6 | Distance vs each parameter | Docked in sensitivity view |

---

## 5. Layout and UI Element Placement

### 5.1. Control Placement

| Area | Content |
|---|---|
| Left panel | `ui.form` with grouped controls + stats display |
| Ribbon (Actions group) | `btn_reset` |
| Ribbon (Analyze group) | `btn_sensitivity` |
| Ribbon (Help group) | `btn_help` |

**Structure of left panel form:**

```
ui.h2('Launch')
  ctrl_v0
  ctrl_theta

ui.h2('Ball')
  ctrl_m
  ctrl_d

ui.h2('Environment')
  ctrl_Cd
  ctrl_rho
  ctrl_g

btn_maximize (button)

ui.h2('Flight Statistics')
  comp_stats
```

### 5.2. Display Element Placement

| Element ID | Docking area | Position / ratio |
|---|---|---|
| `form` (left panel) | `DG.DOCK_TYPE.LEFT` relative to root | ratio `0.6` |
| `view_trajectory` | `DG.DOCK_TYPE.RIGHT` relative to root | ratio `0.7` |
| `view_timeseries` | `DG.DOCK_TYPE.DOWN` relative to trajectory node | ratio `0.5` |
| `view_table` (grid) | Default position (center) | — |

### 5.3. Styles

CSS file: `css/golf-ball.css`. Import: `import '../../css/golf-ball.css'`. All custom classes use the `golf-ball-app-` prefix for isolation.

**Static styles:**

| Element | CSS class(es) | Description |
|---|---|---|
| Stats panel | `.golf-ball-app-stats-panel` | Background, padding, border-radius, gap |
| Stats section header | `.golf-ball-app-stats-panel > label` | Bold, grey color |
| Stats value | `.golf-ball-app-stats-value` | Monospace font, right-aligned |
| Help icon | `.golf-ball-app-help-icon` | font-size, font-weight, margin |

**Dynamic styles:**

| Element | CSS class(es) | Condition | Description |
|---|---|---|---|
| Maximize button | `.golf-ball-app-btn--disabled` | During optimization/sensitivity | Disables pointer events, reduces opacity |

---

## 6. User Feedback

### 6.1. Control Tooltips

| Control ID | Tooltip text | Mechanism |
|---|---|---|
| `ctrl_v0` | Initial ball speed off the clubface. A typical driver produces 50–75 m/s. | `tooltipText` property |
| `ctrl_theta` | Launch angle relative to the ground. Optimal for distance is usually 10–15° with drag, ~45° without. | `tooltipText` property |
| `ctrl_m` | Mass of the golf ball. Regulation balls weigh no more than 45.93 g (0.0459 kg). | `tooltipText` property |
| `ctrl_d` | Diameter of the ball. Regulation minimum is 42.67 mm (0.0427 m). | `tooltipText` property |
| `ctrl_Cd` | Aerodynamic drag coefficient. A smooth ball is ~0.5, a dimpled golf ball ~0.25. | `tooltipText` property |
| `ctrl_rho` | Air density. Sea level ≈ 1.225 kg/m³. Higher altitude or warmer air → lower density → longer flight. | `tooltipText` property |
| `ctrl_g` | Gravitational acceleration. Earth ≈ 9.81 m/s². Moon ≈ 1.62, Mars ≈ 3.72, Jupiter ≈ 24.79. | `tooltipText` property |
| `btn_maximize` | Find the value of the selected parameter that produces the highest ball apex. | `ui.bigButton` third argument |
| `btn_reset` | Reset all parameters to default values | `ui.iconFA` third argument |
| `btn_sensitivity` | Show how each parameter independently affects the total flight distance. | `ui.iconFA` third argument |
| `btn_help` | Learn more about projectile motion with drag. Opens in a separate page. | `ui.icons.help` second argument |

### 6.2. Validators as Feedback

| Input ID | Validation source | Description |
|---|---|---|
| `ctrl_v0` … `ctrl_g` | `validate()` from core | Each input has `addValidator()` calling `validate(getInputs())` and returning the error for its own ID, or `null`. |

### 6.3. Progress Bar

| Task | Progress bar | Type | Cancellation support |
|---|---|---|---|
| `task_primary` | No | — | — |
| `task_maximize_height` | Yes | Indeterminate (`DG.TaskBarProgressIndicator`) | Yes (`pi.onCanceled` → terminates workers) |
| `task_sensitivity` | Yes | Indeterminate (`DG.TaskBarProgressIndicator`) | Yes (`pi.onCanceled` → terminates workers) |

---

## 7. Validation

### 7.1. Primary Pipeline Validation

| Rule ID | Condition (invalid) | Affected inputs (ID) | Error message |
|---|---|---|---|
| `val_01` | `v0 ≤ 0` | `ctrl_v0` | "Initial speed must be positive" |
| `val_02` | `theta ≤ 0 or theta ≥ 90` | `ctrl_theta` | "Launch angle must be between 0° and 90° (exclusive)" |
| `val_03` | `m ≤ 0` | `ctrl_m` | "Ball mass must be positive" |
| `val_04` | `d ≤ 0` | `ctrl_d` | "Ball diameter must be positive" |
| `val_05` | `Cd ≤ 0` | `ctrl_Cd` | "Drag coefficient must be positive" |
| `val_06` | `rho ≤ 0` | `ctrl_rho` | "Air density must be positive" |
| `val_07` | `g ≤ 0` | `ctrl_g` | "Gravitational acceleration must be positive" |

Validation order: val_01 through val_07 (all independent, checked in sequence).

Returned map format: `Map<InputId, string>` where `InputId = 'ctrl_v0' | ... | 'ctrl_g'`.

### 7.2. Secondary Task Validation

`task_maximize_height` reuses primary pipeline validation. If `validate(inputs).size > 0`, optimization is rejected with `grok.shell.error('Invalid parameters. Fix inputs before optimizing.')`.

`task_sensitivity` reuses primary pipeline validation similarly.

---

## 8. Main Pipeline

### 8.1. Primary Pipeline

| Step | Description |
|---|---|
| Parameter input | User sets values through 7 `ui.input.float` controls → `getInputs()` converts to `GolfBallParams` |
| Validation | `validate(inputs)` returns `Map<InputId, string>`. On errors: mark invalid inputs with `d4-invalid` class, call `clearResults()` |
| Computation | `solve(inputs)` — synchronous ODE solution via `mrt()`, truncated at ground impact |
| Result display | `updateDataFrame(result)` creates new DataFrame, assigns to view and viewers; `updateStatsPanel(result)` updates labels |

Reactive trigger: `onValueChanged` with debounce 50 ms.

Error behavior on validation failure: clear results (empty DataFrame).

Error behavior on computation failure: clear results + `grok.shell.error(msg)`.

### 8.2. Secondary Pipelines

#### Task pipeline: `task_maximize_height`

| Step | Description |
|---|---|
| Trigger | `btn_maximize` click → opens dialog |
| Custom UI | Dialog with `ui.input.choice` for parameter selection |
| Validation | `validate(getInputs())` — reuses primary validation |
| Computation | Generate 1000 grid points → distribute across WebWorkers → collect results → find max height |
| Result display | N/A — no dedicated display |
| Feedback to primary | Optimal value written to corresponding control via batch update |

#### Task pipeline: `task_sensitivity`

| Step | Description |
|---|---|
| Trigger | `btn_sensitivity` click |
| Custom UI | N/A |
| Validation | `validate(getInputs())` — reuses primary validation |
| Computation | 6 WebWorkers: each sweeps one parameter (100 steps), solves ODE for each, records flight distance |
| Result display | New `DG.TableView` with DataFrame and line charts |
| Feedback to primary | No |

### 8.3. Common Pipeline Aspects

#### Control Behavior During Computations

| Pipeline / task | Controls blocked | Which controls |
|---|---|---|
| `task_primary` | No | — (computation < 100ms) |
| `task_maximize_height` | `btn_maximize` only | `btn_maximize` disabled via CSS class |
| `task_sensitivity` | `btn_maximize` + all inputs | All `ctrl_*` disabled via `input.enabled = false`; `btn_maximize` disabled via CSS class |

#### Computation Error Handling

| Pipeline / task | Strategy | Notification method |
|---|---|---|
| `task_primary` | Clear results | `grok.shell.error` |
| `task_maximize_height` (partial failures) | Skip failed points, use valid results | `grok.shell.warning` |
| `task_maximize_height` (all fail) | Abort | `grok.shell.error` |
| `task_sensitivity` (partial failures) | Skip failed sweeps | `grok.shell.warning` |

### 8.4. Computation Blocking and Batch Input Updates

| Scenario | Source (task) | Target controls (ID) | Blocked pipelines |
|---|---|---|---|
| Optimization result write | `task_maximize_height` | One of `ctrl_v0`, `ctrl_theta`, `ctrl_Cd` | Primary (blocked) |
| Reset to defaults | `btn_reset` | All `ctrl_*` | Primary (blocked) |
| Format setting on init | Coordinator | All `ctrl_*` | Primary (blocked) |

Reactivity mode during batch update:

| Scenario | Reactivity mode |
|---|---|
| All above | `computationsBlocked = true` → writes → `computationsBlocked = false` → single `runPrimary()` |

---

## 9. Reactivity and Dependencies Between Inputs

### 9.1. Dependency Graph

| Source (input ID) | Target | Reaction type | Logic |
|---|---|---|---|
| Any `ctrl_*` | `comp_stats` | Display update | Recalculate flight stats |
| Any `ctrl_*` | All viewers | Reactive update | Rerun `task_primary`, update charts and table |

### 9.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| All `ctrl_*` | debounce | 50 |

---

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via `ui.form` controls.

Initial state: all controls initialized with defaults. `task_primary` runs automatically on initialization.

### 10.2. Loading from Resources

N/A — no external data loading.

### 10.3. Results Table Lifecycle

Update strategy: **DataFrame replacement** — on each recomputation a new `DG.DataFrame` is created and assigned to `view.dataFrame` and all viewers.

```
1. Application initialization
   → solve with DEFAULTS → new DG.DataFrame with columns [Time, x, y, v]
   → DataFrame added to TableView

2. task_primary completion
   → new DataFrame created from solution arrays
   → assigned to view.dataFrame, trajectory.dataFrame, timeseries.dataFrame

3. task_maximize_height completion
   → Optimal value written to control (batch update)
   → task_primary runs once → new DataFrame

4. btn_reset press
   → All controls reset to defaults
   → task_primary runs → new DataFrame
```

---

## 11. Error Handling Beyond Computations

| Error type | Strategy | Notification method |
|---|---|---|
| ODE solver error | Clear results, show message | `grok.shell.error` |
| Worker creation error | Abort task | `grok.shell.error` |
| Partial worker errors | Skip failed points, use valid results | `grok.shell.warning` |
| All workers fail | Abort task | `grok.shell.error` |

---

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

| Subscription | Event | Cleanup mechanism |
|---|---|---|
| View removal | `grok.events.onViewRemoved` | `sub.unsubscribe()` in cleanup handler |
| Optimization cancel | `pi.onCanceled` | `cancelSub.unsubscribe()` after completion |

### 12.2. Worker Termination

| Worker pool | Created in | Termination mechanism |
|---|---|---|
| Optimization workers | `runMaximizeHeight()` | `w.terminate()` for each worker in `terminateWorkers()` |
| Sensitivity workers | `runSensitivity()` | `w.terminate()` for each worker in `terminateWorkers()` |

---

## 13. Application Closure

On view close (`grok.events.onViewRemoved`), the coordinator performs:

- [x] All event subscriptions unsubscribed (`subs[]` iteration)
- [x] All web workers terminated (`terminateWorkers()`)
- [x] Pending debounce timer cleared (`clearTimeout(debounceTimer)`)

---

## 14. Accessibility and UX

### 14.1. Keyboard Shortcuts

N/A — no custom keyboard shortcuts.

### 14.2. Context Menus

N/A — no custom context menus.

### 14.3. Undo / Redo

Not supported. The `btn_reset` button resets all controls to default values as the only rollback mechanism.

---

## 15. Testing

### 15.1. Computational Part (Core)

| File | Categories | Test count | Description |
|---|---|---|---|
| `src/tests/golf-ball-math-tests.ts` | Math: GB func, Math: GB Solve properties, Math: MRT solver | 14 | ODE RHS verification, solve properties, MRT reference problems |
| `src/tests/golf-ball-api-tests.ts` | API: Validation | 18 | Validation rules for all 7 inputs, defaults, boundaries, multiple errors |

Tests are run via `grok test` (entry point: `src/package-test.ts`).

### 15.2. Inputs

| Category | Coverage | Description |
|---|---|---|
| Boundary values | `v_bnd_v0` (`v0=10`), `v_bnd_theta` (`theta=1`) | Lower allowed bounds |
| Invalid values | Zero and negative/out-of-range for each of 7 parameters (14 tests) | |
| Multiple simultaneous errors | `v0=0, theta=0, m=0, Cd=0` → ≥ 4 errors | |
| Valid defaults | All defaults pass validation | |

### 15.3. Mathematical Verification

#### Level 1 Verification

**Formula/equation verification:**

| Test category | Test count | What is verified | Reference source |
|---|---|---|---|
| Math: GB func | 4 | ODE RHS: concrete (v0, θ, m, d, Cd, ρ, g, vx, vy) → expected (dvx/dt, dvy/dt) | Manual calculation |

**Output property verification:**

| Test category | Test count | Properties verified |
|---|---|---|
| Math: GB Solve properties | 8 | Non-empty arrays, non-negativity of y, initial conditions, v consistency, summary stats, zero-drag analytical match |

#### Level 2 Verification

**Numerical method verification:**

| Test category | Test count | Reference problems | Tolerance | Source |
|---|---|---|---|---|
| Math: MRT solver | 2 | Non-stiff 1D, Stiff 1D | Max error < 0.1 | Chapra & Canale textbook |
