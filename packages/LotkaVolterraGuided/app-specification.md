# Application Specification: Lotka-Volterra Predator-Prey Simulation

## 1. General Information

| Field | Value |
|---|---|
| Application name | Lotka-Volterra Predator-Prey Simulation |
| Package | LotkaVolterraGuided |
| Entry function | `lotkaVolterraApp()` |
| Brief description | Interactive ODE simulation of the Lotka-Volterra predator-prey model with real-time visualization, phase portrait, and brute-force parameter optimization via web workers. |
| Main view | `DG.TableView` |

---

## 2. Computational Tasks (Core)

### 2.1. Task List

| Task ID | Name | Pipeline type | Trigger |
|---|---|---|---|
| `task_primary` | Lotka-Volterra RK4 solution | Primary (reactive) | Any input change |
| `task_optimize` | Optimize max prey (grid search) | Secondary (on demand) | Button `btn_optimize` |

### 2.2. Task Description: `task_primary`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Lotka-Volterra RK4 solution |
| Description | Numerical solution of the 2D ODE system `dx/dt = αx − βxy`, `dy/dt = δxy − γy` using RK4, returning trajectories x(t) and y(t), equilibrium points, and summary statistics |
| Synchronicity | Synchronous |
| Execution environment | Main thread |
| Parallelization | No |
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
| 1 | `α=1.0, β=0.1, δ=0.075, γ=1.5, x₀=10, y₀=5` | `dx/dt = 1.0·10 − 0.1·10·5 = 5.0`, `dy/dt = 0.075·10·5 − 1.5·5 = −3.75` | Manual calculation |
| 2 | At equilibrium `x*=γ/δ=20, y*=α/β=10` | `dx/dt = 0`, `dy/dt = 0` | Manual calculation |
| 3 | `α=1.0, β=0.1, δ=0.075, γ=1.5, x=30, y=4` | `dx/dt = 1.0·30 − 0.1·30·4 = 18`, `dy/dt = 0.075·30·4 − 1.5·4 = 3.0` | Manual calculation |

**Level 2 — full formalization:**

- Complete mathematical formulation: classical Lotka-Volterra predator-prey ODE system with constant coefficients. No boundary conditions (IVP).
- Analytical properties: non-trivial equilibrium at (γ/δ, α/β); solutions are periodic orbits around the equilibrium; conserved quantity H = δx + βy − γ ln(x) − α ln(y).
- Numerical method: RK4 (4th-order Runge-Kutta), explicit method, O(h⁴) local truncation error, suitable for non-stiff problems.

**Input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `alpha` | `number` | 1/time | `> 0` | Prey birth rate |
| `beta` | `number` | 1/(ind·time) | `> 0` | Predation rate |
| `delta` | `number` | 1/(ind·time) | `> 0` | Predator growth efficiency |
| `gamma` | `number` | 1/time | `> 0` | Predator death rate |
| `x0` | `number` | individuals | `> 0` | Initial prey population |
| `y0` | `number` | individuals | `> 0` | Initial predator population |
| `T` | `number` | time | `> 0` | Total simulation time |

**Output data:**

| Parameter | Type | Description |
|---|---|---|
| `t` | `Float64Array` | Time values |
| `x` | `Float64Array` | Prey population values |
| `y` | `Float64Array` | Predator population values |
| `xStar` | `number` | Equilibrium prey: γ/δ |
| `yStar` | `number` | Equilibrium predator: α/β |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | External library | `diff-grok` v1.2.0, function `rk4(task: ODEs)`. RK4 is a 4th-order explicit method suitable for non-stiff ODEs. | [diff-grok README](https://github.com/datagrok-ai/diff-grok) |
| 2 | Custom method | Equilibrium: `xStar = gamma / delta`, `yStar = alpha / beta` | — |
| 3 | Custom method | Summary stats: `maxPrey = max(x)`, `maxPredators = max(y)`, `stepCount = t.length` | — |

**Library call:**

```typescript
import { ODEs, rk4 } from 'diff-grok';

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
const solution = rk4(task);
```

**ODE right-hand side reference examples:**

| # | α | β | δ | γ | x | y | Expected dx/dt | Expected dy/dt | Derivation |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.0 | 0.1 | 0.075 | 1.5 | 10 | 5 | 5.0 | −3.75 | `1.0·10−0.1·10·5=5`, `0.075·10·5−1.5·5=−3.75` |
| 2 | 1.0 | 0.1 | 0.075 | 1.5 | 20 | 10 | 0.0 | 0.0 | At equilibrium (γ/δ=20, α/β=10) |
| 3 | 1.0 | 0.1 | 0.075 | 1.5 | 30 | 4 | 18.0 | 3.0 | `1.0·30−0.1·30·4=18`, `0.075·30·4−1.5·4=3` |
| 4 | 2.0 | 0.5 | 0.3 | 1.0 | 5 | 2 | 5.0 | 1.0 | `2.0·5−0.5·5·2=5`, `0.3·5·2−1.0·2=1` |
| 5 | 0.5 | 0.02 | 0.01 | 0.3 | 50 | 10 | 15.0 | −2.0 | `0.5·50−0.02·50·10=15`, `0.01·50·10−0.3·10=2` |

---

### 2.3. Task Description: `task_optimize`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Optimize Max Prey (grid search) |
| Description | Brute-force grid search over 4 coefficients (α, β, δ, γ) with 10% step size. For each combination, solves the ODE and finds the peak prey population. Returns the combination that maximizes peak prey. |
| Synchronicity | Asynchronous |
| Execution environment | WebWorkers (parallel) |
| Parallelization | Yes — grid points distributed across worker pool of size `Math.max(1, navigator.hardwareConcurrency - 2)` |
| Dependency on other tasks | Uses current x₀, y₀, T from primary pipeline controls |

**Input parameters:**

| Parameter | Type | Description |
|---|---|---|
| `alpha_min, alpha_max` | `number` | Current slider range for α |
| `beta_min, beta_max` | `number` | Current slider range for β |
| `delta_min, delta_max` | `number` | Current slider range for δ |
| `gamma_min, gamma_max` | `number` | Current slider range for γ |
| `x0, y0, T` | `number` | Taken from current state of main UI controls (snapshot) |

**Output data:**

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
| 2 | External library in workers | `diff-grok` v1.2.0, `rk4(task)` — for each grid point in a WebWorker |
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
  → as each task completes: progressBar += 1/totalPoints
  → after all complete: find max(maxPrey) → optimal params
```

### 2.4. Dependencies Between Tasks

```
task_primary  — independent
task_optimize — does not depend on task_primary results;
                after completion, triggers a single run of task_primary
```

---

## 3. Controls

### 3.1. Primary Pipeline Controls

| ID | Name (label) | Control type | Data type | Default | Min | Max | Step | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_alpha` | Prey birth rate α | `ui.input.float` | `number` | `1.0` | `0.1` | `3.0` | — | `0.00` | No | Rate at which prey reproduce. Higher α → faster prey growth in the absence of predators. | Model Coefficients |
| `ctrl_beta` | Predation rate β | `ui.input.float` | `number` | `0.1` | `0.01` | `0.5` | — | `0.000` | No | Rate at which predators consume prey. Higher β → more prey eaten per encounter, reducing prey population faster. | Model Coefficients |
| `ctrl_delta` | Predator efficiency δ | `ui.input.float` | `number` | `0.075` | `0.01` | `0.5` | — | `0.000` | No | Efficiency of converting consumed prey into predator growth. Higher δ → predators grow faster from each prey consumed. | Model Coefficients |
| `ctrl_gamma` | Predator death rate γ | `ui.input.float` | `number` | `1.5` | `0.1` | `3.0` | — | `0.00` | No | Natural death rate of predators. Higher γ → predators die off faster without sufficient prey. | Model Coefficients |
| `ctrl_x0` | Initial prey x₀ | `ui.input.float` | `number` | `10` | `1` | `200` | — | `0.0` | No | Starting prey population at time t=0. | Initial Conditions |
| `ctrl_y0` | Initial predators y₀ | `ui.input.float` | `number` | `5` | `1` | `100` | — | `0.0` | No | Starting predator population at time t=0. | Initial Conditions |
| `ctrl_T` | Simulation time T | `ui.input.float` | `number` | `100` | `10` | `500` | — | `0.0` | No | Total simulation time. Longer T shows more oscillation cycles. | Initial Conditions |

### 3.2. Secondary Task Triggers

| ID | Name / icon | Triggers task | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_optimize` | `Optimize Max Prey` button | `task_optimize` | Run brute-force grid search over all four model coefficients to maximize peak prey population | Always (disabled during optimization) |

### 3.3. Other Buttons and Actions

| ID | Name / icon | Action | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_reset` | `ui.iconFA('undo')` | Reset all controls to default values | Reset all parameters to default values | Always |

### 3.4. Custom UI Components

| Component ID | Brief description | Role |
|---|---|---|
| `comp_equilibrium` | Display panel showing equilibrium (x*, y*) and summary stats | Display |
| `comp_start_marker` | Scatter marker showing the starting point on the phase portrait | Display (via viewer options) |

---

## 4. Validation

### 4.1. Primary Pipeline Validation

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

---

## 5. Reactivity and Input Dependencies

### 5.1. Dependency Graph

| Source (input ID) | Target | Reaction type | Logic |
|---|---|---|---|
| Any `ctrl_*` | `comp_equilibrium` | Display update | Recalculate equilibrium and stats |
| Any `ctrl_*` | All viewers | Reactive update | Rerun `task_primary`, update charts and table |

### 5.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| All `ctrl_*` | debounce | 50 |

---

## 6. Behavior During Computations

### 6.1. Primary Pipeline (`task_primary`)

- Controls blocked: No (computation < 100ms)
- Progress bar: No
- Error behavior: Clear results + `grok.shell.error(msg)`

### 6.2. Secondary Pipeline (`task_optimize`)

- `btn_optimize` blocked: Yes (shows progress bar percentage)
- All other controls: No
- Progress bar: Yes — determinate (0–100%), with percentage shown on button
- Cancellation support: No
- Error behavior: Keep previous results + `grok.shell.error(msg)`

---

## 7. Computation Blocking and Batch Update

| Source (task) | Target controls (ID) | Locked pipelines |
|---|---|---|
| `task_optimize` | `ctrl_alpha`, `ctrl_beta`, `ctrl_delta`, `ctrl_gamma` | Primary (blocked) |

Reactivity mode during batch update: Primary pipeline is paused during writing, then runs once with the new values.

---

## 8. Result Display

### 8.1. Primary Pipeline Display Elements

| ID | Type | Associated output data | Docking location |
|---|---|---|---|
| `view_timeseries` | Datagrok viewer `line chart` | `task_primary.t`, `task_primary.x`, `task_primary.y` | Center area, top |
| `view_phase` | Datagrok viewer `scatter plot` | `task_primary.x`, `task_primary.y` | Center area, bottom |
| `view_table` | `DG.TableView` grid | `task_primary.t`, `task_primary.x`, `task_primary.y` | Right area |
| `comp_equilibrium` | Custom HTMLElement | Equilibrium + stats | Left panel, below controls |

**Details for `view_timeseries`:**

| Property | Value |
|---|---|
| X axis | `t`, label "Time" |
| Y axis | `x(t)` and `y(t)`, two series |
| Series | Prey (x) — line, Predators (y) — line |

**Details for `view_phase`:**

| Property | Value |
|---|---|
| X axis | `x` (prey) |
| Y axis | `y` (predators) |
| Start point | Marked with distinct marker at (x₀, y₀) |

---

## 9. Layout

### 9.1. Control Placement

| Area | Content |
|---|---|
| Left panel | `ui.form` with grouped controls + equilibrium/stats display |
| Center area | Time-series chart (top) + Phase portrait (bottom) |
| Right area | Data table (grid) |

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

---

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via `ui.form` controls.

Initial state: all controls initialized with defaults. `task_primary` runs automatically on initialization.

### 10.2. Results Table Lifecycle

```
1. Application initialization
   → solve with DEFAULTS → DG.DataFrame with columns [t, x, y]
   → DataFrame added to TableView

2. task_primary completion
   → DataFrame updated with new arrays
   → All viewers redrawn reactively

3. task_optimize completion
   → Optimal values written to sliders (batch update)
   → task_primary runs once → DataFrame updated

4. btn_reset press
   → All controls reset to defaults
   → task_primary runs → DataFrame updated
```

---

## 11. Error Handling

| Error type | Strategy | Notification method |
|---|---|---|
| ODE solver error | Clear results, show message | `grok.shell.error` |
| Worker creation error | Abort optimization | `grok.shell.error` |
| Partial worker errors | Skip failed points, use valid results | `grok.shell.warning` |
| All workers fail | Abort optimization | `grok.shell.error` |

---

## 12. Testing

### 12.1. Mathematical Verification

#### ODE Right-Hand Side Verification

| Test ID | Input `(α, β, δ, γ, x, y)` | Expected `(dx/dt, dy/dt)` |
|---|---|---|
| `func_01` | `(1.0, 0.1, 0.075, 1.5, 10, 5)` | `(5.0, −3.75)` |
| `func_02` | `(1.0, 0.1, 0.075, 1.5, 20, 10)` | `(0.0, 0.0)` |
| `func_03` | `(1.0, 0.1, 0.075, 1.5, 30, 4)` | `(18.0, 3.0)` |
| `func_04` | `(2.0, 0.5, 0.3, 1.0, 5, 2)` | `(5.0, 1.0)` |
| `func_05` | `(0.5, 0.02, 0.01, 0.3, 50, 10)` | `(15.0, −2.0)` |

#### Equilibrium Verification

| Test ID | Input `(α, β, δ, γ)` | Expected `(x*, y*)` |
|---|---|---|
| `eq_01` | `(1.0, 0.1, 0.075, 1.5)` | `(20, 10)` |
| `eq_02` | `(2.0, 0.5, 0.3, 1.0)` | `(10/3, 4)` |

#### Output Property Verification (solve)

| Test ID | Description | Verified property |
|---|---|---|
| `solve_01` | Default parameters | Non-empty arrays, equal length |
| `solve_02` | x,y values ≥ 0 | All `x[i] ≥ 0`, `y[i] ≥ 0` |
| `solve_03` | Initial conditions | `x[0] = x₀`, `y[0] = y₀` |
| `solve_04` | Equilibrium: dx/dt ≈ 0, dy/dt ≈ 0 at (x*, y*) | `dx/dt < ε`, `dy/dt < ε` |

### 12.2. Validation Tests

| Test ID | Rule | Input data | Expected result |
|---|---|---|---|
| `v_01` | val_01 | `alpha=0` | Error on `ctrl_alpha` |
| `v_02` | val_02 | `beta=-0.1` | Error on `ctrl_beta` |
| `v_03` | val_03 | `delta=0` | Error on `ctrl_delta` |
| `v_04` | val_04 | `gamma=-1` | Error on `ctrl_gamma` |
| `v_05` | val_05 | `x0=0` | Error on `ctrl_x0` |
| `v_06` | val_06 | `y0=-5` | Error on `ctrl_y0` |
| `v_07` | val_07 | `T=0` | Error on `ctrl_T` |
| `v_def` | all defaults | DEFAULTS | `errors.size = 0` |
| `v_multi` | multiple | `alpha=0, beta=0, x0=0, T=0` | ≥ 4 errors |
