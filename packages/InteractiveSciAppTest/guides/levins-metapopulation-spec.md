# Application Specification: Levins Metapopulation Model

## 1. General Information

| Field | Value |
|---|---|
| Application name | Levins Metapopulation Model |
| Package | MetapopulationDynamics |
| Entry function | `levinsMetapopulationApp()` |
| Brief description | Interactive ODE solution for the Levins model: simulation of occupied patch fraction dynamics p(t) with support for the basic model and the extended model (rescue effect). |
| Main view | `DG.TableView` |

---

## 2. Computational Tasks (Core)

### 2.1. Task List

| Task ID | Name | Pipeline type | Trigger |
|---|---|---|---|
| `task_primary` | Levins model ODE solution | Primary (reactive) | Any input change |
| `task_optimize` | Find optimal m by p(t_end) | Secondary (on demand) | Button `btn_optimize` |

### 2.2. Task Description: `task_primary`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Levins model ODE solution |
| Description | Numerical solution of ODE `dp/dt = m·p·(1−p) − e(p)·p`, returning the trajectory p(t) and the equilibrium value p* |
| Synchronicity | Synchronous |
| Execution environment | Main thread |
| Parallelization | No |
| Dependency on other tasks | No |

**Input parameters:**

| Parameter | Type | Description |
|---|---|---|
| `p0` | `number` | Initial fraction of occupied patches |
| `m` | `number` | Colonization rate |
| `e0` | `number` | Baseline extinction rate |
| `rescueEffect` | `boolean` | Enable rescue effect: `e(p) = e0·(1−p)` |
| `t_start` | `number` | Start of integration interval |
| `t_end` | `number` | End of interval |
| `t_step` | `number` | Grid step |
| `tolerance` | `number` | Numerical tolerance of the MRT method |

**Output data:**

| Parameter | Type | Description |
|---|---|---|
| `t` | `Float64Array` | Values of argument t |
| `p` | `Float64Array` | Values of p(t) |
| `p_star` | `number` | Equilibrium value: `1 − e0/m` (for the basic model) |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | External library | `diff-grok` v1.2.0, function `mrt(task: ODEs)` | [diff-grok README](https://github.com/datagrok-ai/diff-grok) |
| 2 | Custom method | Computation of `p_star = 1 − e0/m`; for rescue effect, analytical equilibrium is not used | See model specification (separate document) |

**Library call:**

```typescript
import { ODEs, mrt } from 'diff-grok';

const task: ODEs = {
  name: 'Levins',
  arg: { name: 't', start: t_start, finish: t_end, step: t_step },
  initial: [p0],
  func: (t, y, out) => {
    const e = rescueEffect ? e0 * (1 - y[0]) : e0;
    out[0] = m * y[0] * (1 - y[0]) - e * y[0];
  },
  tolerance: tolerance,
  solutionColNames: ['p(t)'],
};
const solution = mrt(task);
```

---

### 2.2. Task Description: `task_optimize`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Find optimal m by p(t_end) |
| Description | For 10000 uniformly distributed values of m in [min, max], solves the ODE, computes p(t_end), returns m with the maximum p(t_end) |
| Synchronicity | Asynchronous |
| Execution environment | WebWorkers (parallel) |
| Parallelization | Yes — 10000 independent tasks, worker pool of size `Math.max(1, navigator.hardwareConcurrency - 2)` |
| Dependency on other tasks | Uses current values of all `task_primary` inputs except `m` |

**Input parameters:**

| Parameter | Type | Description |
|---|---|---|
| `m_min` | `number` | Lower bound for m search |
| `m_max` | `number` | Upper bound for m search |
| `p0, e0, rescueEffect` | `number / boolean` | Taken from current state of main UI controls (snapshot) |
| `t_start, t_end, t_step, tolerance` | `number` | Taken from current state of main UI controls (snapshot) |

**Output data:**

| Parameter | Type | Description |
|---|---|---|
| `m_optimal` | `number` | Value of m at which p(t_end) is maximized |
| `p_at_t_end_max` | `number` | Achieved maximum value of p(t_end) |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | Custom method | Generate 10000 points: `m_i = m_min + i · (m_max − m_min) / 9999`, i = 0..9999 | — |
| 2 | External library in workers | `diff-grok` v1.2.0, `mrt(task)` — for each point `m_i` in a WebWorker. Passed through pipeline API: `getIvp2WebWorker(ivp)` | [diff-grok pipeline](https://github.com/datagrok-ai/diff-grok) |
| 3 | Custom method | After receiving all results: `m_optimal = m_i` where `p_end` is maximal | — |
| 4 | Datagrok API | Write `m_optimal` to `ctrl_m` via batch update (section 7) | [Datagrok JS API](https://datagrok.ai/api/js/) |

**Parallelization strategy:**

```
Number of workers = Math.max(1, navigator.hardwareConcurrency - 2)

10000 values of m_i
  → worker pool
  → tasks are distributed to workers as they become available (queue)
  → each worker receives: { m_i, p0, e0, rescueEffect,
                             t_start, t_end, t_step, tolerance }
  → each worker returns: { m_i, p_end }
  → as each task completes: progressBar += 1/10000
  → after all 10000: find max(p_end) → m_optimal
```

### 2.3. Dependencies Between Tasks

```
task_primary  — independent
task_optimize — does not depend on task_primary results;
                after completion, triggers a single run of task_primary
```

---

## 3. Controls

### 3.1. Primary Pipeline Controls

| ID | Name (label) | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_p0` | Initial patch fraction p₀ | `ui.input.float` | `number` | `0.5` | `0.001` | `1` | `0.000` | No | Fraction of patches occupied at time t=0. At p₀=0 the population will not recover — computation is not performed. | Initial condition |
| `ctrl_m` | Colonization rate m | `ui.input.float` | `number` | `0.5` | `0.001` | `100` | `0.000` | No | How quickly free patches are colonized from occupied ones. Units: 1/time. | Parameters |
| `ctrl_e0` | Extinction rate e₀ | `ui.input.float` | `number` | `0.2` | `0.001` | `100` | `0.000` | No | Baseline rate of local subpopulation extinction in a patch. With rescue effect — decreases as p grows. Units: 1/time. | Parameters |
| `ctrl_rescue` | Rescue effect | `ui.input.toggle` | `boolean` | `false` | — | — | — | No | If enabled, the extinction rate depends on p: e(p) = e₀·(1−p). The more occupied patches — the lower the local extinction. | Parameters |
| `ctrl_t_start` | Start t₀ | `ui.input.float` | `number` | `0` | `0` | `10000` | `0.0` | No | Simulation start time. Usually 0. | Argument |
| `ctrl_t_end` | End t_end | `ui.input.float` | `number` | `50` | `0.1` | `10000` | `0.0` | No | End time. Recommended ≥ 5/e₀ so the system reaches equilibrium. | Argument |
| `ctrl_t_step` | Step Δt | `ui.input.float` | `number` | `0.1` | `0.001` | `1000` | `0.000` | No | Numerical solution grid step. Affects chart detail but not stability (MRT is an implicit method). | Argument |
| `ctrl_tolerance` | Tolerance | `ui.input.float` | `number` | `1e-7` | `1e-12` | `1e-2` | `0.##E+0` | No | Numerical tolerance of the MRT method. Smaller — more precise, but slower. Recommended 1e-6 … 1e-9. | Solver |

### 3.2. Secondary Task Triggers

| ID | Name / icon | Triggers task | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_optimize` | `ui.iconFA('search')` | `task_optimize` | Find the value of m that maximizes the fraction of occupied patches at time t_end | Always |

### 3.3. Controls for `task_optimize`

UI type: Datagrok dialog window.

| ID | Name (label) | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text |
|---|---|---|---|---|---|---|---|---|---|
| `dlg_m_min` | Minimum m | `ui.input.float` | `number` | `0.1` | `0.001` | `100` | `0.000` | No | Lower bound of the m search range. Must be less than the maximum value. |
| `dlg_m_max` | Maximum m | `ui.input.float` | `number` | `1.0` | `0.001` | `100` | `0.000` | No | Upper bound of the m search range. Must be greater than the minimum value. |

**Dialog buttons:**

| ID | Name | Action | Availability condition |
|---|---|---|---|
| `dlg_btn_ok` | OK | Close dialog → launch `task_optimize` | Only if there are no validation errors |
| `dlg_btn_cancel` | Cancel | Close dialog, task is not launched | Always |

### 3.4. Other Buttons and Actions

| ID | Name / icon | Action | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_reset` | `ui.iconFA('undo')` | Reset all controls to default values | Reset all parameters to default values | Always |

### 3.5. Custom UI Components

| Component ID | Brief description | Role | UI component specification |
|---|---|---|---|
| `comp_rho_badge` | Indicator for ρ = e₀/m with color status (green / red) | Display | — |

> **Note on `comp_rho_badge`:** displays the current value of ρ = e₀/m in real time. Green — the metapopulation persists (ρ < 1), red — heading toward extinction (ρ ≥ 1). Updated reactively when `ctrl_m` or `ctrl_e0` changes.
>
> **Tooltip:** "Ratio of extinction rate to colonization rate. ρ < 1 — metapopulation persists, ρ ≥ 1 — extinction."

---

## 4. Validation

### 4.1. Primary Pipeline Validation (`task_primary`)

#### Complex Validation Rules

| Rule ID | Condition (invalid) | Affected inputs (ID) | Error message |
|---|---|---|---|
| `val_01` | `p0 ≤ 0` | `ctrl_p0` | "Initial patch fraction must be greater than 0" |
| `val_02` | `p0 > 1` | `ctrl_p0` | "Initial patch fraction cannot exceed 1" |
| `val_03` | `m ≤ 0` | `ctrl_m` | "Colonization rate must be positive" |
| `val_04` | `e0 ≤ 0` | `ctrl_e0` | "Extinction rate must be positive" |
| `val_05` | `m ≤ e0` (when `rescueEffect = false`) | `ctrl_m`, `ctrl_e0` | "Colonization rate must exceed extinction rate (m > e₀). At current values, the metapopulation is heading toward extinction" |
| `val_06` | `t_end ≤ t_start` | `ctrl_t_end`, `ctrl_t_start` | "End of interval must be greater than start" |
| `val_07` | `t_step ≤ 0` | `ctrl_t_step` | "Step must be positive" |
| `val_08` | `t_step ≥ t_end − t_start` | `ctrl_t_step` | "Step must be less than the interval length" |
| `val_09` | `tolerance ≤ 0` | `ctrl_tolerance` | "Tolerance must be positive" |

#### Validation Order

```
1. val_01, val_02  (single, ctrl_p0)
2. val_03          (single, ctrl_m)
3. val_04          (single, ctrl_e0)
4. val_05          (combinatorial ctrl_m + ctrl_e0 — only if val_03 and val_04 passed)
5. val_06          (combinatorial ctrl_t_start + ctrl_t_end)
6. val_07          (single, ctrl_t_step)
7. val_08          (combinatorial ctrl_t_step + ctx — only if val_06 and val_07 passed)
8. val_09          (single, ctrl_tolerance)
```

#### Return Map Format

```
Map<inputId, errorMessage>
```

### 4.2. Validation for `task_optimize`

#### Complex Validation Rules

| Rule ID | Condition (invalid) | Affected inputs (ID) | Error message |
|---|---|---|---|
| `opt_val_01` | `m_min ≤ 0` | `dlg_m_min` | "Colonization rate must be positive" |
| `opt_val_02` | `m_max ≤ 0` | `dlg_m_max` | "Colonization rate must be positive" |
| `opt_val_03` | `m_min ≥ m_max` | `dlg_m_min`, `dlg_m_max` | "Minimum value must be less than maximum" |
| `opt_val_04` | `m_max ≤ e0` (when `rescueEffect = false`) | `dlg_m_max` | "At the current e₀, the entire m range leads to extinction (m ≤ e₀). Increase the maximum or decrease e₀" |

> `opt_val_04` — warning, does not block OK.

#### Validation Order

```
1. opt_val_01  (single)
2. opt_val_02  (single)
3. opt_val_03  (combinatorial — only if 01 and 02 passed)
4. opt_val_04  (combinatorial with external parameter e₀ — only if 01-03 passed)
```

---

## 5. Reactivity and Input Dependencies

### 5.1. Dependency Graph

| Source (input ID) | Target (input IDs / components) | Reaction type | Logic |
|---|---|---|---|
| `ctrl_m` | `comp_rho_badge` | Display update | Recalculate ρ = e₀ / m, update badge value and color |
| `ctrl_e0` | `comp_rho_badge` | Display update | Same |
| `ctrl_m` | `comp_rho_badge` | Availability update | If m ≤ e₀ → badge is red, otherwise green |
| `ctrl_e0` | `comp_rho_badge` | Availability update | Same |
| `ctrl_rescue` | `ctrl_e0` (label + tooltip) | Label update | If `rescue = true` → label changes to "Baseline extinction rate e₀", tooltip adds: "Effective rate: e(p) = e₀·(1−p)" |
| `ctrl_t_start` | `ctrl_t_end` | Range update | Minimum allowed value of `ctrl_t_end` = `t_start + t_step` |
| `ctrl_t_start` | `ctrl_t_step` | Range update | Maximum allowed value of `ctrl_t_step` = `t_end − t_start` |
| `ctrl_t_end` | `ctrl_t_step` | Range update | Same: `ctrl_t_step` ≤ `t_end − t_start` |
| `dlg_m_min` | `dlg_m_max` | Range update | Minimum allowed value of `dlg_m_max` = `m_min + ε` |
| `dlg_m_max` | `dlg_m_min` | Range update | Maximum allowed value of `dlg_m_min` = `m_max − ε` |
| `dlg_m_min` | `dlg_btn_ok` | Availability update | OK button is available only if there are no validation errors `opt_val_01–03` |
| `dlg_m_max` | `dlg_btn_ok` | Availability update | Same |

### 5.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| `ctrl_p0` | debounce | 50 |
| `ctrl_m` | debounce | 50 |
| `ctrl_e0` | debounce | 50 |
| `ctrl_t_start` | debounce | 50 |
| `ctrl_t_end` | debounce | 50 |
| `ctrl_t_step` | debounce | 50 |
| `ctrl_tolerance` | debounce | 50 |
| `ctrl_rescue` | none | — |
| `dlg_m_min` | debounce | 50 |
| `dlg_m_max` | debounce | 50 |

---

## 6. Behavior During Computations

### 6.1. Primary Pipeline (`task_primary`)

#### Control Locking

| Control ID | Locked | Note |
|---|---|---|
| `ctrl_p0` | No | Computation takes < 100 ms — locking is unnecessary |
| `ctrl_m` | No | Same |
| `ctrl_e0` | No | Same |
| `ctrl_rescue` | No | Same |
| `ctrl_t_start` | No | Same |
| `ctrl_t_end` | No | Same |
| `ctrl_t_step` | No | Same |
| `ctrl_tolerance` | No | Same |
| `btn_optimize` | No | — |
| `btn_reset` | No | — |

#### Progress Bar

| Field | Value |
|---|---|
| Display | No |
| Type | — |
| Cancellation support | No |

#### Error Behavior

Selected strategy: **Reset results + message**.

| Strategy | Description |
|---|---|
| Reset results | Clear the p(t) chart and the p* value |
| Message | Show an inline message below the chart with the error text |

### 6.2. Secondary Pipeline (`task_optimize`)

#### Control Locking

| Control ID | Locked | Note |
|---|---|---|
| `btn_optimize` | Yes | Re-launch is not possible until completion or cancellation |
| All other controls | No | The main UI remains fully accessible |

> Changing inputs while `task_optimize` is running triggers a `task_primary` recalculation in normal mode but does not affect the already running search — it uses the parameter snapshot from the moment OK was pressed.

#### Progress Bar

| Field | Value |
|---|---|
| Display | Yes |
| Type | Determinate (0–100%, +1/10000 for each completed worker) |
| Cancellation support | Yes — Cancel button terminates all active workers |

#### Error Behavior

Selected strategy: **Last valid + message**.

| Strategy | Description |
|---|---|
| Last valid | The value of `ctrl_m` is not changed |
| Message | `grok.shell.error` with error text |

---

## 7. Computation Locking and Batch Update

### 7.1. Batch Update Scenarios

| Source (task) | Target controls (ID) | Locked pipelines |
|---|---|---|
| `task_optimize` | `ctrl_m` | Primary |

### 7.2. Reactivity Mode During Batch Update

| Scenario | Reactivity mode |
|---|---|
| Writing `m_optimal` → `ctrl_m` | Primary pipeline is paused during writing, then runs once with the new value |

---

## 8. Result Display

### 8.1. Primary Pipeline Display Elements

| ID | Type | Associated output data | Docking location |
|---|---|---|---|
| `view_p_t` | Datagrok viewer `line chart` | `task_primary.t`, `task_primary.p` | Main area |
| `view_p_star` | Custom HTMLElement `comp_rho_badge` | `task_primary.p_star` | Left panel, below "Parameters" group controls |

**Details for `view_p_t`:**

| Property | Value |
|---|---|
| X axis | `t`, label "Time" |
| Y axis | `p(t)`, range `[0, 1]`, label "Fraction of occupied patches" |
| Series | `p(t)` — main trajectory, solid line |
| Update | Reactive — redrawn on each `task_primary` completion |

**Color coding:**

The `p` column of the results table receives `colorCoding` by value:

| Range | Color | Meaning |
|---|---|---|
| `p < e₀/m` | Red | Extinction threat zone |
| `p ≥ e₀/m` | Green | Persistence zone |

> The threshold value `e₀/m` is recalculated and updated in `colorCoding` on each change of `ctrl_m` or `ctrl_e0`.

**Tooltip for `p` column header:**

The `p` column header in the grid receives a tooltip explaining the color coding:

```
Fraction of occupied patches p(t).
Color: green — persistence zone (p ≥ e₀/m),
red — extinction threat zone (p < e₀/m).
Threshold: e₀/m = {current ρ value}.
```

> The tooltip text is updated reactively when `ctrl_m` or `ctrl_e0` changes.

### 8.2. Display Elements for `task_optimize`

Where results are displayed: toast notification + writing to the main control.

| ID | Type | Associated output data | Placement |
|---|---|---|---|
| `view_optimize_result` | `grok.shell.info` | `task_optimize.m_optimal`, `task_optimize.p_at_t_end_max` | Datagrok platform toast notification |
| `ctrl_m` | `ui.input.float` (main UI control) | `task_optimize.m_optimal` | Left panel, "Parameters" group |

**Toast notification text:**
```
Optimal m = {m_optimal}
p(t_end) = {p_at_t_end_max}
```

**Behavior after writing the result:**

```
1. m_optimal → ctrl_m  (batch update, section 7)
2. task_primary runs once with the new m
3. view_p_t is redrawn with the new trajectory
4. grok.shell.info is shown
```

---

## 9. Layout

### 9.1. Control Placement

| Area | Content | Docking | Ratio |
|---|---|---|---|
| Left panel | `ui.form` with groups separated by `ui.h2` headers | `DG.DOCK_TYPE.LEFT` | `0.2` |
| Toolbar | `btn_optimize`, `btn_reset` | — | — |
| Main area (grid) | `DG.TableView` (results table) | Default | — |
| Right area | `view_p_t` (line chart) | `DG.DOCK_TYPE.RIGHT` (relative to grid) | `0.5` |

**Structure of `ui.form` in the left panel:**

```
ui.h2('Initial condition')
  ctrl_p0

ui.h2('Parameters')
  ctrl_m
  ctrl_e0
  ctrl_rescue
  comp_rho_badge

ui.h2('Argument')
  ctrl_t_start
  ctrl_t_end
  ctrl_t_step

ui.h2('Solver')
  ctrl_tolerance
```

### 9.2. Display Element Placement

| Element ID | Type | Area | Note |
|---|---|---|---|
| `view_p_t` | Viewer `line chart` | Right area (dock right, ratio `0.5`) | Docked to the right of the grid, splitting space 50/50 |
| `view_optimize_result` | `grok.shell.info` | — | Datagrok platform toast, placed automatically by the platform |

---

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via `ui.form` controls (section 3).

Initial application state: all controls are initialized with default values from section 3.1 at the time `levinsMetapopulationApp()` is called. `task_primary` runs automatically immediately after initialization.

### 10.2. Loading from Resources

External data loading is not supported.

| Trigger | Resource | Format | Mapping to inputs |
|---|---|---|---|
| — | — | — | — |

### 10.3. Results Table Lifecycle

```
1. Application initialization
   → an empty DG.DataFrame is created with columns [t, p]
   → the DataFrame is added to the TableView

2. task_primary completion
   → the DataFrame is updated: columns [t, p] are overwritten with new Float64Array
   → colorCoding for column p is recalculated (threshold e₀/m)
   → view_p_t is redrawn reactively

3. task_optimize completion
   → m_optimal is written to ctrl_m (batch update, section 7)
   → task_primary runs once → DataFrame is updated per step 2
   → grok.shell.info is shown

4. btn_reset press
   → all controls are reset to default values
   → task_primary runs → DataFrame is updated per step 2

5. task_primary error
   → the DataFrame is cleared (columns [t, p] are zeroed out)
   → view_p_t displays an empty chart
   → an inline error message is shown below the chart
```

### 10.4. Data Lifecycle for `task_optimize`

```
1. Pressing OK in the dialog
   → a snapshot of current parameters { p0, e0, rescueEffect,
     t_start, t_end, t_step, tolerance } is captured
   → an array of 10000 m_i values is generated

2. Worker execution
   → each worker receives { m_i, snapshot }
   → each worker returns { m_i, p_end }
   → intermediate results are not stored anywhere

3. All workers complete
   → m_optimal = m_i at max(p_end) is computed
   → the array { m_i, p_end } is freed from memory

4. Cancellation (progress bar Cancel)
   → all active workers are terminated
   → intermediate results are discarded
   → ctrl_m is not changed
```

---

## 11. Error Handling Beyond Computations

| Error type | Strategy | Notification method |
|---|---|---|
| `diff-grok` initialization error (library failed to load) | Lock `btn_optimize` and all controls, show message | `grok.shell.error`: "Failed to load the solver library. Reload the page." |
| WebWorker creation error (browser does not support or limit exceeded) | Abort `task_optimize`, do not change `ctrl_m` | `grok.shell.error`: "Failed to start parallel computations. Try again later." |
| Worker terminated with error (one or more m_i points) | Skip the point, continue remaining workers, consider only valid results | `grok.shell.warning`: "{N} out of 10000 points were not computed. Result obtained from {10000−N} points." |
| All 10000 workers terminated with error | Abort `task_optimize`, do not change `ctrl_m` | `grok.shell.error`: "Failed to compute any points. Check the parameters." |
| Invalid application state (snapshot contains invalid values) | Abort `task_optimize` before creating workers | `grok.shell.error`: "Internal error: invalid task parameters. Check inputs and retry." |
| Error updating `ctrl_m` after `task_optimize` completion | Show result via `grok.shell.info`, do not write to control | `grok.shell.warning`: "Optimal m = {m_optimal}, but the field could not be updated automatically. Enter the value manually." |

---

## 12. UX

### 12.1. Keyboard Shortcuts

| Combination | Action |
|---|---|
| — | — |

### 12.2. Context Menus

| Context (element) | Menu items |
|---|---|
| — | — |

### 12.3. Undo / Redo

Supported: **No**.

---

## 13. Testing

### 13.1. Computation Part (Core)

#### Tests for `getEquilibrium`

| Test ID | Description | Input data | Expected result |
|---|---|---|---|
| `eq_01` | Basic model equilibrium | `m=0.5, e0=0.2, rescue=false` | `p* = 0.6` |
| `eq_02` | p* = 0 when m ≤ e0 | `m=0.2, e0=0.5, rescue=false` | `p* = 0` |
| `eq_03` | p* = 0 when m = e0 | `m=0.5, e0=0.5, rescue=false` | `p* = 0` |
| `eq_04` | NaN with rescue effect | `m=0.5, e0=0.2, rescue=true` | `NaN` |

#### Tests for `solve` (task_primary)

| Test ID | Description | Input data | Verified property |
|---|---|---|---|
| `solve_01` | Default parameters | DEFAULTS | `t.length > 0`, `p.length > 0`, `t.length = p.length` |
| `solve_02` | p values in [0, 1] | DEFAULTS | All `p[i] ∈ [0, 1]` |
| `solve_03` | p(0) = p0 | DEFAULTS | `p[0] ≈ 0.5` |
| `solve_04` | t starts at t_start | DEFAULTS | `t[0] = 0` |
| `solve_05` | Convergence to p* | `m=0.5, e0=0.2, t_end=200` | `p(t_end) ≈ p*` (tolerance 0.01) |
| `solve_06` | Rescue effect: p in [0, 1] | `m=0.3, e0=0.5, rescue=true, t_end=100` | All `p[i] ∈ [0, 1]` |
| `solve_07` | Higher m → higher p(t_end) | `m=0.5` vs `m=1.0`, `e0=0.2` | `p_end(m=1) > p_end(m=0.5)` |
| `solve_08` | Custom p0 | `p0=0.9` | `p[0] ≈ 0.9` |

### 13.2. Validation for task_primary

| Test ID | Rule | Input data | Expected result |
|---|---|---|---|
| `v_01a` | val_01 | `p0=0` | Error on `ctrl_p0` |
| `v_01b` | val_01 | `p0=-1` | Error on `ctrl_p0` |
| `v_02` | val_02 | `p0=1.1` | Error on `ctrl_p0` |
| `v_p0_lo` | boundary | `p0=0.001` | No error on `ctrl_p0` |
| `v_p0_hi` | boundary | `p0=1` | No error on `ctrl_p0` |
| `v_03a` | val_03 | `m=0` | Error on `ctrl_m` |
| `v_03b` | val_03 | `m=-0.5` | Error on `ctrl_m` |
| `v_04a` | val_04 | `e0=0` | Error on `ctrl_e0` |
| `v_04b` | val_04 | `e0=-0.1` | Error on `ctrl_e0` |
| `v_05a` | val_05 | `m=0.5, e0=0.5, rescue=false` | Error on `ctrl_m` |
| `v_05b` | val_05 | `m=0.3, e0=0.5, rescue=false` | Error on `ctrl_m` |
| `v_05c` | val_05 (rescue) | `m=0.3, e0=0.5, rescue=true` | No error |
| `v_05d` | val_05 dependency | `m=0, e0=0.5` | Message val_03, not val_05 |
| `v_06a` | val_06 | `t_start=10, t_end=10` | Error on `ctrl_t_end` and `ctrl_t_start` |
| `v_06b` | val_06 | `t_start=10, t_end=5` | Error on `ctrl_t_end` |
| `v_07a` | val_07 | `t_step=0` | Error on `ctrl_t_step` |
| `v_07b` | val_07 | `t_step=-0.1` | Error on `ctrl_t_step` |
| `v_08a` | val_08 | `t_step=50 (= t_end-t_start)` | Error on `ctrl_t_step` |
| `v_08b` | val_08 | `t_step=100 (> t_end-t_start)` | Error on `ctrl_t_step` |
| `v_08c` | val_08 dependency val_06 | `t_start=10, t_end=5, t_step=100` | No error on `ctrl_t_step` |
| `v_08d` | val_08 dependency val_07 | `t_step=-1` | Message val_07, not val_08 |
| `v_09a` | val_09 | `tolerance=0` | Error on `ctrl_tolerance` |
| `v_09b` | val_09 | `tolerance=-1e-7` | Error on `ctrl_tolerance` |
| `v_def` | all defaults | DEFAULTS | `errors.size = 0` |
| `v_multi` | multiple | `p0=0, m=0, e0=0, t_step=0, tolerance=0` | ≥ 4 errors |

### 13.3. Validation for task_optimize

| Test ID | Rule | Input data | Expected result |
|---|---|---|---|
| `ov_01` | opt_val_01 | `m_min=0` | Error on `dlg_m_min` |
| `ov_02` | opt_val_02 | `m_max=-1` | Error on `dlg_m_max` |
| `ov_03a` | opt_val_03 | `m_min=0.5, m_max=0.5` | Error on `dlg_m_min` |
| `ov_03b` | opt_val_03 | `m_min=1.0, m_max=0.5` | Error on `dlg_m_min` |
| `ov_04a` | opt_val_04 | `m_max=0.1, e0=0.2, rescue=false` | Warning ≠ null |
| `ov_04b` | opt_val_04 (rescue) | `m_max=0.1, e0=0.2, rescue=true` | Warning = null |
| `ov_valid` | valid | `m_min=0.1, m_max=1.0, e0=0.2` | `errors.size = 0`, `warning = null` |
