# Application Specification: Lotka-Volterra Predator-Prey Simulation

<!-- TOC for partial reading -->
<!-- Section 1: General Architecture — lines 15–324 -->
<!-- Section 2: Main View — lines 325–333 -->
<!-- Section 3: Controls — lines 334–385 -->
<!-- Section 4: Display Elements — lines 386–437 -->
<!-- Section 5: Layout — lines 438–502 -->
<!-- Section 6: User Feedback — lines 503–536 -->
<!-- Section 7: Validation — lines 537–560 -->
<!-- Section 8: Pipeline — lines 561–636 -->
<!-- Section 9: Reactivity — lines 637–653 -->
<!-- Sections 10–15: Data, Errors, Resources, Closure, UX, Testing — lines 654–827 -->

## 1. General Architecture

### 1.0. General Information

| Field | Value |
|---|---|
| Application name | Lotka-Volterra Predator-Prey Simulation |
| Package | LotkaVolterraGuided |
| Entry function | `lotkaVolterraApp()` |
| Brief description | Interactive ODE simulation of the Lotka-Volterra predator-prey model with real-time visualization, phase portrait, and brute-force parameter optimization via web workers. |
| Main view type | `DG.TableView` |

### 1.1. Core

#### Task List

| Task ID | Name | Pipeline type | Trigger | Synchronicity | Execution environment | Parallelization |
|---|---|---|---|---|---|---|
| `task_primary` | Lotka-Volterra MRT solution | Primary (reactive) | Any input change | Sync | Main thread | No |
| `task_optimize` | Optimize max prey (grid search) | Secondary (on demand) | Button `btn_optimize` | Async | WebWorkers (parallel) | Yes — grid points distributed across worker pool of size `Math.max(1, navigator.hardwareConcurrency - 2)` |
| `task_equilibrium` | Explore equilibrium point | Secondary (on demand) | Button `btn_analyze` → dialog → Run | Async | WebWorker (single) | No |

#### Dependencies Between Tasks

```
task_primary     — independent
task_optimize    — does not depend on task_primary results;
                   after completion, triggers a single run of task_primary
task_equilibrium — independent; opens results in a new TableView
```

#### Task: `task_primary`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Lotka-Volterra MRT solution |
| Description | Numerical solution of the 2D ODE system `dx/dt = αx − βxy`, `dy/dt = δxy − γy` using MRT, returning trajectories x(t) and y(t), equilibrium points, and summary statistics |
| Dependency on other tasks | No |

**Computation Formulas and Model**

**Level 1 — required minimum:**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| `x` | Prey population | individuals | `x ≥ 0` |
| `y` | Predator population | individuals | `y ≥ 0` |
| `α` (alpha) | Prey birth rate | 1/time | `α > 0` |
| `β` (beta) | Predation rate | 1/(individuals·time) | `β > 0` |
| `δ` (delta) | Predator growth efficiency | 1/(individuals·time) | `δ > 0` |
| `γ` (gamma) | Predator death rate | 1/time | `γ > 0` |
| `x₀` | Initial prey population | individuals | `x₀ > 0` |
| `y₀` | Initial predator population | individuals | `y₀ > 0` |
| `T` | Total simulation time | time | `T > 0` |

Relationships:

```
dx/dt = α·x − β·x·y
dy/dt = δ·x·y − γ·y
```

Equilibrium points:
- Trivial: (0, 0)
- Non-trivial: (x*, y*) = (γ/δ, α/β)

Output properties (invariants):

| Property | Description |
|---|---|
| `x(t) ≥ 0` for all `t` | Prey population is non-negative |
| `y(t) ≥ 0` for all `t` | Predator population is non-negative |
| `x(0) = x₀` | Initial prey condition preserved |
| `y(0) = y₀` | Initial predator condition preserved |
| Phase portrait is a closed orbit (for default params) | Conservation of Hamiltonian `H = δx + βy − γ ln(x) − α ln(y)` |

Reference examples:

| # | Inputs | Expected output | Source |
|---|---|---|---|
| 1 | `α=1.0, β=0.1, δ=0.075, γ=1.5, x₀=10, y₀=5` | `dx/dt = 5.0`, `dy/dt = −3.75` | Manual calculation |
| 2 | At equilibrium `x*=γ/δ=20, y*=α/β=10` | `dx/dt = 0`, `dy/dt = 0` | Manual calculation |
| 3 | `α=1.0, β=0.1, δ=0.075, γ=1.5, x=30, y=4` | `dx/dt = 18`, `dy/dt = 3.0` | Manual calculation |

**Level 2 — full formalization:**

- Complete mathematical formulation: classical Lotka-Volterra predator-prey ODE system with constant coefficients. No boundary conditions (IVP).
- Analytical properties: non-trivial equilibrium at (γ/δ, α/β); solutions are periodic orbits around the equilibrium; conserved quantity H = δx + βy − γ ln(x) − α ln(y).
- Numerical method: MRT (adaptive multi-rate), suitable for both stiff and non-stiff problems.

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `alpha` | `number` | 1/time | `> 0` | Prey birth rate |
| `beta` | `number` | 1/(ind·time) | `> 0` | Predation rate |
| `delta` | `number` | 1/(ind·time) | `> 0` | Predator growth efficiency |
| `gamma` | `number` | 1/time | `> 0` | Predator death rate |
| `x0` | `number` | individuals | `> 0` | Initial prey population |
| `y0` | `number` | individuals | `> 0` | Initial predator population |
| `T` | `number` | time | `> 0` | Total simulation time |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `t` | `Float64Array` | Time values |
| `x` | `Float64Array` | Prey population values |
| `y` | `Float64Array` | Predator population values |
| `xStar` | `number` | Equilibrium prey: γ/δ |
| `yStar` | `number` | Equilibrium predator: α/β |
| `maxPrey` | `number` | Maximum prey population: max(x) |
| `maxPredators` | `number` | Maximum predator population: max(y) |
| `stepCount` | `number` | Number of time steps: t.length |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | External library | `diff-grok` v1.2.0, function `mrt(task: ODEs)`. MRT is an adaptive multi-rate method suitable for both stiff and non-stiff ODEs. | [diff-grok README](https://github.com/datagrok-ai/diff-grok) |
| 2 | Custom method | Equilibrium: `xStar = gamma / delta`, `yStar = alpha / beta` | — |
| 3 | Custom method | Summary stats: `maxPrey = max(x)`, `maxPredators = max(y)`, `stepCount = t.length` | — |

**Library call:**

```typescript
import { ODEs, mrt } from 'diff-grok';

const task: ODEs = {
  name: 'LotkaVolterra',
  arg: { name: 't', start: 0, finish: T, step: 0.05 },
  initial: [x0, y0],
  func: (t, y, out) => {
    out[0] = alpha * y[0] - beta * y[0] * y[1];
    out[1] = delta * y[0] * y[1] - gamma * y[1];
  },
  tolerance: 1e-6,
  solutionColNames: ['x(t)', 'y(t)'],
};
const solution = mrt(task);
```

**ODE right-hand side reference examples:**

| # | α | β | δ | γ | x | y | Expected dx/dt | Expected dy/dt | Derivation |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.0 | 0.1 | 0.075 | 1.5 | 10 | 5 | 5.0 | −3.75 | `1.0·10−0.1·10·5=5`, `0.075·10·5−1.5·5=−3.75` |
| 2 | 1.0 | 0.1 | 0.075 | 1.5 | 20 | 10 | 0.0 | 0.0 | At equilibrium (γ/δ=20, α/β=10) |
| 3 | 1.0 | 0.1 | 0.075 | 1.5 | 30 | 4 | 18.0 | 3.0 | `1.0·30−0.1·30·4=18`, `0.075·30·4−1.5·4=3` |
| 4 | 2.0 | 0.5 | 0.3 | 1.0 | 5 | 2 | 5.0 | 1.0 | `2.0·5−0.5·5·2=5`, `0.3·5·2−1.0·2=1` |
| 5 | 0.5 | 0.02 | 0.01 | 0.3 | 50 | 10 | 15.0 | 2.0 | `0.5·50−0.02·50·10=15`, `0.01·50·10−0.3·10=2` |

**Execution environment constraint:** main thread (synchronous, < 100 ms). The core does not use Datagrok API — it depends only on `diff-grok` and custom methods.

---

#### Task: `task_optimize`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Optimize Max Prey (grid search) |
| Description | Brute-force grid search over 4 coefficients (α, β, δ, γ) with 10% step size. For each combination, solves the ODE and finds the peak prey population. Returns the combination that maximizes peak prey. |
| Dependency on other tasks | Uses current x₀, y₀, T from primary pipeline controls |

**Task input parameters:**

| Parameter | Type | Description |
|---|---|---|
| `alpha_min, alpha_max` | `number` | Current slider range for α |
| `beta_min, beta_max` | `number` | Current slider range for β |
| `delta_min, delta_max` | `number` | Current slider range for δ |
| `gamma_min, gamma_max` | `number` | Current slider range for γ |
| `x0, y0, T` | `number` | Taken from current state of main UI controls (snapshot) |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `alpha_opt` | `number` | Optimal α |
| `beta_opt` | `number` | Optimal β |
| `delta_opt` | `number` | Optimal δ |
| `gamma_opt` | `number` | Optimal γ |
| `maxPrey` | `number` | Maximum peak prey achieved |

**Computation implementation:**

| Step | Implementation method | Details |
|---|---|---|
| 1 | Custom method | Generate grid: 11 values per parameter (0%, 10%, ..., 100% of range), total 11⁴ = 14641 points |
| 2 | External library in workers | `diff-grok` v1.2.0, `mrt(task)` — for each grid point in a WebWorker |
| 3 | Custom method | Find grid point with maximum peak prey |
| 4 | Datagrok API | Write optimal values to sliders via batch update |

**Parallelization strategy:**

```
Number of workers = Math.max(1, navigator.hardwareConcurrency - 2)

14641 grid points (11^4)
  → worker pool
  → tasks distributed to workers as they become available (queue)
  → each worker receives: { alpha, beta, delta, gamma, x0, y0, T }
  → each worker returns: { alpha, beta, delta, gamma, maxPrey }
  → after all complete: find max(maxPrey) → optimal params
```

---

#### Task: `task_equilibrium`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Explore Equilibrium Point |
| Description | Varies a user-selected model parameter across its full range (100000 steps) while keeping other parameters fixed. For each value, computes the non-trivial equilibrium point (x* = γ/δ, y* = α/β). Opens results in a new TableView with a line chart. |
| Dependency on other tasks | Uses current parameter values from primary pipeline controls |

**Task input parameters:**

| Parameter | Type | Description |
|---|---|---|
| `paramName` | `ParamName` | Selected parameter to vary (`alpha`, `beta`, `delta`, `gamma`, `x0`, `y0`, `T`) |
| `paramMin, paramMax` | `number` | Range of the selected parameter (from `RANGES`) |
| `steps` | `number` | Number of steps (100000) |
| `baseParams` | `LotkaVolterraParams` | Current values of all parameters (snapshot) |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `paramValues` | `number[]` | Values of the varied parameter |
| `xStar` | `number[]` | Equilibrium prey for each parameter value |
| `yStar` | `number[]` | Equilibrium predator for each parameter value |

**Computation implementation:**

| Step | Implementation method | Details |
|---|---|---|
| 1 | Custom method in worker | Generate 100000 values of the selected parameter via linspace |
| 2 | Custom method in worker | For each value: `x* = gamma / delta`, `y* = alpha / beta` (with the varied param substituted) |
| 3 | Datagrok API | Create DataFrame with columns `[paramName, Prey*, Predator*]`, open new TableView, add line chart |

---

### 1.2. Ports

#### Task ports: `task_primary`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `LotkaVolterraParams` | `{ alpha, beta, delta, gamma, x0, y0, T }` | 7 validated numeric parameters |
| Output | `LotkaVolterraSolution` | `{ t, x, y, xStar, yStar, maxPrey, maxPredators, stepCount }` | Solution arrays + derived stats |

#### Task ports: `task_optimize`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `WorkerTask[]` | `{ alpha, beta, delta, gamma, x0, y0, T }[]` | Grid of parameter combinations |
| Output | `WorkerResult[]` | `{ alpha, beta, delta, gamma, maxPrey, error? }[]` | Per-point results |

#### Task ports: `task_equilibrium`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `EquilibriumTask` | `{ paramName, paramMin, paramMax, steps, baseParams }` | Parameter to vary + base values |
| Output | `EquilibriumResult` | `{ paramValues[], xStar[], yStar[], error? }` | Equilibrium arrays |

#### Application-level ports

| Port | Used | Description |
|---|---|---|
| Progress | Yes | `DG.TaskBarProgressIndicator` for `task_optimize` and `task_equilibrium` (indeterminate, cancelable) |
| Cancellation | Yes | `pi.onCanceled` → terminates workers, resolves promises |
| Data | No | No external data loading |

### 1.3. Adapters

| Adapter | Used | Implementation |
|---|---|---|
| UI adapter | Yes | Datagrok inputs (`ui.input.float` × 7), `ui.bigButton`, `ui.iconFA`, `ui.icons.help`, `ui.dialog`, `ui.input.choice` |
| Display adapter | Yes | Datagrok viewers (line chart, scatter plot, grid), custom `comp_equilibrium` (section 3.5) |
| Worker adapter | Yes | Web Worker wrappers for `task_optimize` (`workers/optimize-worker.ts`) and `task_equilibrium` (`workers/equilibrium-worker.ts`) |
| Progress adapter | Yes | `DG.TaskBarProgressIndicator` (indeterminate, cancelable) |
| Data adapter | No | No external data loading |

### 1.4. Coordinator

The coordinator is the `lotkaVolterraApp()` function in `app.ts`:

- **Input listening:** each `ui.input.float` calls `debouncedRun()` from `onValueChanged`.
- **Reactivity management:** all 7 controls trigger the same primary pipeline; no cascading dependencies between controls.
- **Validation trigger:** `runPrimary()` calls `validate(inputs)` before computation; validators are also attached to each input via `addValidator()` for inline hints.
- **Control state during computations:** `computationsBlocked` flag prevents `runPrimary()` during batch updates (format setting, optimization result write, reset).
- **Computation blocking:** `computationsBlocked = true` before batch writes, `= false` after, then single `runPrimary()`.
- **Results to display:** `updateDataFrame()` creates new `DG.DataFrame` and assigns to `view.dataFrame`, `lineChart.dataFrame`, `scatterPlot.dataFrame`; `updateStatsPanel()` updates labels via `textContent`.
- **Resource lifecycle:** `subs[]` collects subscriptions; `activeWorkers[]` tracks workers; both cleaned up in `onViewRemoved` handler.

### 1.5. Independence Principle

Confirmed. Input controls and their reactivity (debounce, validation hints) function independently of the computational core. The core receives a ready `LotkaVolterraParams` object; it has no knowledge of UI elements. The `validate()` function is a pure function in `core.ts`.

---

## 2. Main View

| Field | Value |
|---|---|
| View type | `DG.TableView` |
| Description | Table view with docked line chart, scatter plot, grid, and left-panel form with controls and stats |

---

## 3. Controls (Inputs)

### 3.1. Primary Pipeline Controls

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_alpha` | Prey birth rate α | `ui.input.float` | `number` | `1.0` | `0.1` | `3.0` | `0.00` | No | Rate at which prey reproduce. Higher α → faster prey growth in the absence of predators. | Model Coefficients |
| `ctrl_beta` | Predation rate β | `ui.input.float` | `number` | `0.1` | `0.01` | `0.5` | `0.000` | No | Rate at which predators consume prey. Higher β → more prey eaten per encounter, reducing prey population faster. | Model Coefficients |
| `ctrl_delta` | Predator efficiency δ | `ui.input.float` | `number` | `0.075` | `0.01` | `0.5` | `0.000` | No | Efficiency of converting consumed prey into predator growth. Higher δ → predators grow faster from each prey consumed. | Model Coefficients |
| `ctrl_gamma` | Predator death rate γ | `ui.input.float` | `number` | `1.5` | `0.1` | `3.0` | `0.00` | No | Natural death rate of predators. Higher γ → predators die off faster without sufficient prey. | Model Coefficients |
| `ctrl_x0` | Initial prey x₀ | `ui.input.float` | `number` | `10` | `1` | `200` | `0.0` | No | Starting prey population at time t=0. | Initial Conditions |
| `ctrl_y0` | Initial predators y₀ | `ui.input.float` | `number` | `5` | `1` | `100` | `0.0` | No | Starting predator population at time t=0. | Initial Conditions |
| `ctrl_T` | Simulation time T | `ui.input.float` | `number` | `100` | `10` | `500` | `0.0` | No | Total simulation time. Longer T shows more oscillation cycles. | Initial Conditions |

### 3.2. Secondary Task Triggers

| ID | Label / icon | Launches task | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_optimize` | `ui.bigButton('Optimize')` | `task_optimize` | Run brute-force grid search over all four model coefficients to maximize peak prey population | Always (disabled during optimization and equilibrium exploration) |
| `btn_analyze` | `ui.iconFA('chart-line')` on ribbon | `task_equilibrium` (via dialog) | Explore equilibrium point | Always |

### 3.3. Secondary Task Controls

`task_optimize` has no dedicated UI; it uses the current slider ranges and primary control values as input.

#### Task controls: `task_equilibrium`

UI type: Datagrok dialog (`ui.dialog('Explore Equilibrium Point')`).

| ID | Label | Control type | Data type | Default | Nullable | Tooltip text |
|---|---|---|---|---|---|---|
| `dlg_param` | Parameter to vary | `ui.input.choice` | `ParamName` | `'alpha'` | No | Select a model parameter to vary across its full range while keeping all other parameters fixed. |

Dialog buttons:
- **Cancel** — closes dialog
- **Run** — closes dialog, launches `task_equilibrium`. Tooltip: "Run equilibrium point dependency exploration for the selected input"

### 3.4. Other Buttons and Actions

| ID | Label / icon | Action | Tooltip text | Ribbon group | Availability condition |
|---|---|---|---|---|---|
| `btn_reset` | `ui.iconFA('undo')` | Reset all controls to default values | Reset all parameters to default values | Actions | Always |
| `btn_help` | `ui.icons.help` (`.lotka-volterra-app-help-icon`) | Open Wikipedia article on Lotka-Volterra equations in a new tab | Learn more about the Predator-Prey model. Opens in a separate page. | Help | Always |

### 3.5. Custom UI Components

| Component ID | Brief description | Role |
|---|---|---|
| `comp_equilibrium` | Display panel showing equilibrium and max values. Built with `ui.divV`/`ui.label` components; values updated via `textContent`. Section headers ("Equilibrium", "Max") are bold labels; values ("Prey* = ...", "Predator* = ...") are right-aligned via `.lotka-volterra-app-stats-value`. Styled via `.lotka-volterra-app-stats-panel`. | Display |

---

## 4. Result Display Elements

### 4.1. Primary Pipeline Display Elements

| ID | Type | Associated output data | Docking location |
|---|---|---|---|
| `view_timeseries` | Datagrok viewer `line chart` | `Time`, `Prey`, `Predators` | Right area, top |
| `view_phase` | Datagrok viewer `scatter plot` | `Prey`, `Predators` | Right area, bottom |
| `view_table` | `DG.TableView` grid | `Time`, `Prey`, `Predators` | Center (default position) |
| `comp_equilibrium` | Custom HTMLElement | Equilibrium + stats | Left panel, below controls |

**Details for `view_timeseries`:**

| Property | Value |
|---|---|
| X axis | `Time` |
| Y axis | `Prey` and `Predators`, two series |
| Series | Prey — line, Predators — line |

**Details for `view_phase`:**

| Property | Value |
|---|---|
| X axis | `Prey` |
| Y axis | `Predators` |
| Size | `Predators` |
| Color | `Prey` |

**Details for `view_table`:**

| Property | Value |
|---|---|
| Color coding on `Prey` | Linear: dark red (min) → yellow (mid) → dark green (max), applied to text. `col.meta.colors.setLinear(...)` + `grid.col('Prey')!.isTextColorCoded = true` |
| Color coding on `Predators` | Linear: dark red (min) → yellow (mid) → dark green (max), applied to text. `col.meta.colors.setLinear(...)` + `grid.col('Predators')!.isTextColorCoded = true` |
| Reapplied | On every DataFrame replacement |
| Column header tooltip (`Prey`, `Predators`) | Custom `ui.divV` (`.lotka-volterra-app-col-tooltip`) shown via `view.grid.onCellTooltip`: `ui.h2(colName)` + range line (`.lotka-volterra-app-tooltip-range`) with **min** (`.lotka-volterra-app-tooltip-min`, bold, dark red) `...` **max** (`.lotka-volterra-app-tooltip-max`, bold, dark green) showing current column min/max values. All styles in CSS. |

### 4.2. Secondary Task Display Elements

`task_optimize` writes results back to primary pipeline controls (batch update), which triggers `task_primary` and updates all primary display elements.

#### Task display: `task_equilibrium`

Results are displayed in a **new TableView** (separate from main view).

| ID | Type | Associated output data | Placement |
|---|---|---|---|
| `eq_view` | `DG.TableView` | DataFrame `[paramName, Prey*, Predator*]` | New view |
| `eq_chart` | Datagrok viewer `line chart` | X: `paramName`, Y: `Prey*` and `Predator*`. `lineWidth: 3`. | Docked `DG.DOCK_TYPE.RIGHT` relative to grid, ratio `0.75` |

---

## 5. Layout and UI Element Placement

### 5.1. Control Placement

| Area | Content |
|---|---|
| Left panel | `ui.form` with grouped controls + equilibrium/stats display |
| Ribbon (Actions group) | `btn_reset` |
| Ribbon (Analyze group) | `btn_analyze` |
| Ribbon (Help group) | `btn_help` |

**Structure of left panel form:**

```
ui.h2('Model Coefficients')
  ctrl_alpha
  ctrl_beta
  ctrl_delta
  ctrl_gamma

ui.h2('Initial Conditions')
  ctrl_x0
  ctrl_y0
  ctrl_T

btn_optimize (button with progress)

ui.h2('Equilibrium & Stats')
  comp_equilibrium
```

### 5.2. Display Element Placement

| Element ID | Docking area | Position / ratio |
|---|---|---|
| `form` (left panel) | `DG.DOCK_TYPE.LEFT` relative to root | ratio `0.6` |
| `view_timeseries` | `DG.DOCK_TYPE.RIGHT` relative to root | ratio `0.7` |
| `view_phase` | `DG.DOCK_TYPE.DOWN` relative to line chart node | ratio `0.5` |
| `view_table` (grid) | Default position (center) | — |

### 5.3. Styles

CSS file: `css/lotka-volterra.css`. Import: `import '../../css/lotka-volterra.css'`. All custom classes use the `lotka-volterra-app-` prefix for isolation.

**Static styles:**

| Element | CSS class(es) | Description |
|---|---|---|
| Stats panel | `.lotka-volterra-app-stats-panel` | Background, padding, border-radius, gap |
| Stats section header | `.lotka-volterra-app-stats-panel > label` | Bold, grey color |
| Stats value | `.lotka-volterra-app-stats-value` | Monospace font, right-aligned |
| Help icon | `.lotka-volterra-app-help-icon` | `font-size: large`, `font-weight: bold`, right margin |
| Tooltip container | `.lotka-volterra-app-col-tooltip` | Vertical gap between header and range |
| Tooltip range row | `.lotka-volterra-app-tooltip-range` | Flex row, baseline-aligned |
| Tooltip min label | `.lotka-volterra-app-tooltip-min` | Bold, dark red (`#8B0000`) |
| Tooltip max label | `.lotka-volterra-app-tooltip-max` | Bold, dark green (`#006400`) |

**Dynamic styles:**

| Element | CSS class(es) | Condition | Description |
|---|---|---|---|
| Optimize button | `.lotka-volterra-app-optimize-btn--disabled` | During optimization | Disables pointer events, reduces opacity |

---

## 6. User Feedback

### 6.1. Control Tooltips

| Control ID | Tooltip text | Mechanism |
|---|---|---|
| `ctrl_alpha` | Rate at which prey reproduce. Higher α → faster prey growth in the absence of predators. | `tooltipText` property |
| `ctrl_beta` | Rate at which predators consume prey. Higher β → more prey eaten per encounter, reducing prey population faster. | `tooltipText` property |
| `ctrl_delta` | Efficiency of converting consumed prey into predator growth. Higher δ → predators grow faster from each prey consumed. | `tooltipText` property |
| `ctrl_gamma` | Natural death rate of predators. Higher γ → predators die off faster without sufficient prey. | `tooltipText` property |
| `ctrl_x0` | Starting prey population at time t=0. | `tooltipText` property |
| `ctrl_y0` | Starting predator population at time t=0. | `tooltipText` property |
| `ctrl_T` | Total simulation time. Longer T shows more oscillation cycles. | `tooltipText` property |
| `btn_optimize` | Run brute-force grid search over all four model coefficients to maximize peak prey population | `ui.bigButton` third argument |
| `btn_reset` | Reset all parameters to default values | `ui.iconFA` third argument |
| `btn_analyze` | Explore equilibrium point | `ui.iconFA` third argument |
| `btn_help` | Learn more about the Predator-Prey model. Opens in a separate page. | `ui.icons.help` second argument |

### 6.2. Validators as Feedback

| Input ID | Validation source | Description |
|---|---|---|
| `ctrl_alpha` … `ctrl_T` | `validate()` from core | Each input has `addValidator()` calling `validate(getInputs())` and returning the error for its own ID, or `null`. Displays inline hint on invalid value. |

### 6.3. Progress Bar

| Task | Progress bar | Type | Cancellation support |
|---|---|---|---|
| `task_primary` | No | — | — |
| `task_optimize` | Yes | Indeterminate (`DG.TaskBarProgressIndicator`) | Yes (`pi.onCanceled` → terminates workers) |
| `task_equilibrium` | Yes | Indeterminate (`DG.TaskBarProgressIndicator`) | Yes (`pi.onCanceled` → terminates worker) |

---

## 7. Validation

### 7.1. Primary Pipeline Validation

| Rule ID | Condition (invalid) | Affected inputs (ID) | Error message |
|---|---|---|---|
| `val_01` | `alpha ≤ 0` | `ctrl_alpha` | "Prey birth rate must be positive" |
| `val_02` | `beta ≤ 0` | `ctrl_beta` | "Predation rate must be positive" |
| `val_03` | `delta ≤ 0` | `ctrl_delta` | "Predator efficiency must be positive" |
| `val_04` | `gamma ≤ 0` | `ctrl_gamma` | "Predator death rate must be positive" |
| `val_05` | `x0 ≤ 0` | `ctrl_x0` | "Initial prey population must be positive" |
| `val_06` | `y0 ≤ 0` | `ctrl_y0` | "Initial predator population must be positive" |
| `val_07` | `T ≤ 0` | `ctrl_T` | "Simulation time must be positive" |

Validation order: val_01 through val_07 (all independent, checked in sequence).

Returned map format: `Map<InputId, string>` where `InputId = 'ctrl_alpha' | ... | 'ctrl_T'`.

### 7.2. Secondary Task Validation

N/A — `task_optimize` reuses the primary pipeline validation. If `validate(inputs).size > 0`, optimization is rejected with `grok.shell.error('Invalid parameters. Fix inputs before optimizing.')`.

---

## 8. Main Pipeline

### 8.1. Primary Pipeline

| Step | Description |
|---|---|
| Parameter input | User sets values through 7 `ui.input.float` controls → `getInputs()` converts to `LotkaVolterraParams` |
| Validation | `validate(inputs)` returns `Map<InputId, string>`. On errors: mark invalid inputs with `d4-invalid` class, call `clearResults()` |
| Computation | `solve(inputs)` — synchronous ODE solution via `mrt()` |
| Result display | `updateDataFrame(result)` creates new DataFrame, assigns to view and viewers; `updateStatsPanel(result)` updates labels |

Reactive trigger: `onValueChanged` with debounce 50 ms.

Error behavior on validation failure: clear results (empty DataFrame).

Error behavior on computation failure: clear results + `grok.shell.error(msg)`.

### 8.2. Secondary Pipelines

#### Task pipeline: `task_optimize`

| Step | Description |
|---|---|
| Trigger | `btn_optimize` click |
| Custom UI | N/A — no dialog |
| Validation | `validate(getInputs())` — reuses primary validation |
| Computation | Generate 14641 grid points → distribute across WebWorkers → collect results → find max `maxPrey` |
| Result display | N/A — no dedicated display |
| Feedback to primary | Optimal `alpha`, `beta`, `delta`, `gamma` written to `ctrl_alpha`, `ctrl_beta`, `ctrl_delta`, `ctrl_gamma` via batch update |

#### Task pipeline: `task_equilibrium`

| Step | Description |
|---|---|
| Trigger | `btn_analyze` click → opens dialog |
| Custom UI | Dialog with `ui.input.choice` for parameter selection (section 3.3) |
| Validation | N/A — parameter choice is always valid |
| Computation | Single WebWorker: vary selected parameter across its range (100000 steps), compute equilibrium for each |
| Result display | New `DG.TableView` with DataFrame `[paramName, Prey*, Predator*]` + line chart (section 4.2) |
| Feedback to primary | No |

### 8.3. Common Pipeline Aspects

#### Control Behavior During Computations

| Pipeline / task | Controls blocked | Which controls |
|---|---|---|
| `task_primary` | No | — (computation < 100ms) |
| `task_optimize` | `btn_optimize` only | `btn_optimize` disabled via `.lotka-volterra-app-optimize-btn--disabled` |
| `task_equilibrium` | All inputs + `btn_optimize` | All `ctrl_*` disabled via `input.enabled = false`; `btn_optimize` disabled via CSS class |

#### Computation Error Handling

| Pipeline / task | Strategy | Notification method |
|---|---|---|
| `task_primary` | Clear results | `grok.shell.error` |
| `task_optimize` (partial failures) | Skip failed points, use valid results | `grok.shell.warning` |
| `task_optimize` (all fail) | Abort optimization | `grok.shell.error` |
| `task_equilibrium` | Abort, no results view | `grok.shell.error` |

### 8.4. Computation Blocking and Batch Input Updates

| Scenario | Source (task) | Target controls (ID) | Blocked pipelines |
|---|---|---|---|
| Optimization result write | `task_optimize` | `ctrl_alpha`, `ctrl_beta`, `ctrl_delta`, `ctrl_gamma` | Primary (blocked) |
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
| Any `ctrl_*` | `comp_equilibrium` | Display update | Recalculate equilibrium and stats |
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

Update strategy: **DataFrame replacement** — on each recomputation a new `DG.DataFrame` is created and assigned to `view.dataFrame` and all viewers. This avoids in-place mutation and ensures viewers always reflect a consistent snapshot.

```
1. Application initialization
   → solve with DEFAULTS → new DG.DataFrame with columns [Time, Prey, Predators]
   → DataFrame added to TableView

2. task_primary completion
   → new DataFrame created from solution arrays
   → assigned to view.dataFrame, lineChart.dataFrame, scatterPlot.dataFrame
   → All viewers redrawn reactively

3. task_optimize completion
   → Optimal values written to sliders (batch update)
   → task_primary runs once → new DataFrame created and assigned

4. btn_reset press
   → All controls reset to defaults
   → task_primary runs → new DataFrame created and assigned
```

---

## 11. Error Handling Beyond Computations

| Error type | Strategy | Notification method |
|---|---|---|
| ODE solver error | Clear results, show message | `grok.shell.error` |
| Worker creation error | Abort optimization | `grok.shell.error` |
| Partial worker errors | Skip failed points, use valid results | `grok.shell.warning` |
| All workers fail | Abort optimization | `grok.shell.error` |

---

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

| Subscription | Event | Cleanup mechanism |
|---|---|---|
| View removal | `grok.events.onViewRemoved` | `sub.unsubscribe()` in cleanup handler |
| Optimization cancel | `pi.onCanceled` | `cancelSub.unsubscribe()` after `Promise.allSettled` |

All subscriptions are collected in `subs[]` and unsubscribed when the view is removed.

### 12.2. Worker Termination

| Worker pool | Created in | Termination mechanism |
|---|---|---|
| Optimization workers | `runOptimization()` | `w.terminate()` for each worker in `terminateWorkers()`, called on completion, cancellation, and view close |
| Equilibrium worker | `runEquilibriumExploration()` | `worker.terminate()` on completion, error, or cancellation via `pi.onCanceled` |

---

## 13. Application Closure

On view close (`grok.events.onViewRemoved`), the coordinator performs:

- [x] All event subscriptions unsubscribed (`subs[]` iteration)
- [x] All web workers terminated (`terminateWorkers()`)
- [x] Pending debounce timer cleared (`clearTimeout(debounceTimer)`)
- [x] No open secondary task dialogs (optimization has no dialog)

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
| `src/tests/lotka-volterra-math-tests.ts` | Math: LV func, Math: LV Equilibrium, Math: LV Solve properties, Math: MRT solver | 16 | ODE RHS verification, equilibrium, solve properties, MRT solver reference problems |
| `src/tests/lotka-volterra-api-tests.ts` | API: Validation | 18 | Validation rules for all 7 inputs, defaults, boundaries, multiple errors |

Tests are run via `grok test` (entry point: `src/package-test.ts`).

### 15.2. Inputs

| Category | Coverage | Description |
|---|---|---|
| Boundary values | `v_bnd_alpha` (`alpha=0.1`), `v_bnd_x0` (`x0=1`) | Lower allowed bounds |
| Invalid values | `v_01a`–`v_07b` (14 tests) | Zero and negative for each of 7 parameters |
| Multiple simultaneous errors | `v_multi` | `alpha=0, beta=0, x0=0, T=0` → ≥ 4 errors |
| Valid defaults | `v_def` | All defaults pass validation |

### 15.3. Mathematical Verification

#### Level 1 Verification

**Formula/equation verification:**

| Test category | Test count | What is verified | Reference source |
|---|---|---|---|
| Math: LV func | 5 | ODE RHS: concrete (α, β, δ, γ, x, y) → expected (dx/dt, dy/dt) | Manual calculation |
| Math: LV Equilibrium | 2 | Equilibrium point: (α, β, δ, γ) → expected (x*, y*) | Manual calculation |

**Output property verification:**

| Test category | Test count | Properties verified |
|---|---|---|
| Math: LV Solve properties | 7 | Non-empty arrays, non-negativity, initial conditions, t starts at 0, equilibrium stability, summary stats, custom initial conditions |

#### Level 2 Verification

**Numerical method verification:**

| Test category | Test count | Reference problems | Tolerance | Source |
|---|---|---|---|---|
| Math: MRT solver | 2 | Non-stiff 1D (`dy/dt = 4·exp(0.8t) − 0.5y`), Stiff 1D (`dy/dt = −1000y + 3000 − 2000·exp(−t)`) | Max error < 0.1 | Chapra & Canale textbook |

**Convergence verification:**

| Status | Description |
|---|---|
| Not yet covered | Could verify that reducing tolerance decreases discrepancy between solutions |

**Asymptotic/equilibrium behavior:**

| Status | Description |
|---|---|
| Implemented | `solve_05`: start at (x*, y*), verify `x(T) ≈ x*`, `y(T) ≈ y*` after T=50 |

#### Validation Tests

| Test ID | Rule | Input data | Expected result |
|---|---|---|---|
| `v_01a` | val_01 | `alpha=0` | Error on `ctrl_alpha` |
| `v_01b` | val_01 | `alpha=-1` | Error on `ctrl_alpha` |
| `v_02a` | val_02 | `beta=0` | Error on `ctrl_beta` |
| `v_02b` | val_02 | `beta=-0.1` | Error on `ctrl_beta` |
| `v_03a` | val_03 | `delta=0` | Error on `ctrl_delta` |
| `v_03b` | val_03 | `delta=-0.05` | Error on `ctrl_delta` |
| `v_04a` | val_04 | `gamma=0` | Error on `ctrl_gamma` |
| `v_04b` | val_04 | `gamma=-1` | Error on `ctrl_gamma` |
| `v_05a` | val_05 | `x0=0` | Error on `ctrl_x0` |
| `v_05b` | val_05 | `x0=-5` | Error on `ctrl_x0` |
| `v_06a` | val_06 | `y0=0` | Error on `ctrl_y0` |
| `v_06b` | val_06 | `y0=-5` | Error on `ctrl_y0` |
| `v_07a` | val_07 | `T=0` | Error on `ctrl_T` |
| `v_07b` | val_07 | `T=-10` | Error on `ctrl_T` |
| `v_def` | all defaults | DEFAULTS | `errors.size = 0` |
| `v_bnd_alpha` | valid boundary | `alpha=0.1` | No error on `ctrl_alpha` |
| `v_bnd_x0` | valid boundary | `x0=1` | No error on `ctrl_x0` |
| `v_multi` | multiple | `alpha=0, beta=0, x0=0, T=0` | ≥ 4 errors |
