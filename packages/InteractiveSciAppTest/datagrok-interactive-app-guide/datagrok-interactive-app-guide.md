# Guide: Building Interactive Scientific Web Applications on Datagrok

> **End-to-end example:** application **Levins Metapopulation Model** (see `example/` directory).
> All file references are relative to the `guides/` directory.
> Application specification: `example/levins-metapopulation-spec.md`.

## 1. General Architecture

An application is a Datagrok package function registered with the `//tags: app` comment.

► **Implementation:** `example/code/src/package.ts` — function `levinsMetapopulationModelApp`, registered via `//name: Levins Metapopulation Model` / `//tags: app`. Delegates logic to `example/code/src/levins/app.ts` → `levinsMetapopulationApp()`.

The application architecture follows the "ports and adapters" pattern (hexagonal architecture). The application consists of three layers and a coordinator.

### 1.1. Core

Computation logic. The core knows nothing about the UI and does not depend on how data is obtained or how results are displayed.

The core contains one or more computational tasks. Each task is a self-contained unit of computation with its own inputs and outputs. For each task, the following are defined independently:

- Input parameters and their types.
- Output data and their format.
- Synchronicity: synchronous or asynchronous.
- Execution environment: main thread or web worker.
- Parallelization: whether parallelization is possible and the strategy.
- Complex validation logic for input parameters.

► **Implementation:** `example/code/src/levins/core.ts` — core with two tasks:
- **`task_primary`** (synchronous, main thread) — `solve(inputs)` solves the Levins model ODE and returns the trajectory `p(t)`.
- **`task_optimize`** (asynchronous, web workers) — searches for the optimal `m` across 10,000 points, executed in parallel in a worker pool.

Examples of computational tasks within a single application:

- **Primary task** — solving the ODE based on model parameters. Triggered reactively when inputs change.
- **Secondary task** — model sensitivity analysis. Triggered by an explicit user action (button), has its own set of parameters.

Tasks can be independent or linked (the secondary task uses results from the primary one).

► **Implementation:** `task_optimize` does not depend on the results of `task_primary`, but uses the current input values of the primary pipeline (snapshot). After optimization completes, the result is written to the `ctrl_m` control via batch update, which triggers a single run of `task_primary`.

#### Computation Formulas and Model

Each computational task is based on a defined transformation of inputs into outputs — from a single formula to a complex system of equations. In this guide, the term "model" refers to any such definition: individual formulas, chains of transformations, ODE/PDE systems, optimization problems, or statistical procedures. The model exists independently of its implementation (custom code, external library, or platform API). The model definition is the primary source of verification criteria: what cannot be defined cannot be verified.

The model definition has two levels.

**Level 1 — required minimum (before implementation):**

- **Variables:** name, meaning, units of measurement, valid domain (e.g., `p ∈ (0, 1]`, `m > 0`).
- **Relationships:** equations, recurrences, algorithmic steps that connect inputs to outputs. The notation must be unambiguous — another developer must be able to independently implement the same computation from this description.
- **Output properties:** constraints that must hold on the result — bounds, monotonicity, symmetry, conservation laws, limiting/degenerate cases. These properties directly become verification tests.
- **Reference examples:** for each computational path (each mode, branch, or regime of the model), at least one concrete input → expected output pair with the source (manual calculation, literature, reference implementation).

**Level 2 — full formalization (can be developed incrementally alongside the implementation):**

- **Complete mathematical formulation:** equations, initial/boundary conditions, parameterization. For ODE/PDE — the system in explicit form.
- **Analytical properties:** equilibria, asymptotic behavior, stability conditions, bifurcation points.
- **Numerical method justification:** why this particular method is chosen, its properties (stability, order of accuracy, applicability to stiff/non-stiff problems), reference to literature or documentation.
Both Level 1 and Level 2 content can be placed directly in the main application specification or extracted into a separate model specification document — depending on complexity. For simple models (a few formulas), the application specification is sufficient. For complex models (multi-step pipelines, multiple computational paths, extensive reference data), a separate document avoids cluttering the main specification. The main application specification then includes a brief model description and a link to the model specification.

Level 2 need not be complete before implementation begins, but must be complete before the computational part is considered verified.

► **Implementation:**
- **Level 1** is defined in the application specification (`example/levins-metapopulation-spec.md`): variables `p, m, e₀` with units and domains; ODE `dp/dt = m·p·(1−p) − e(p)·p` with two modes (`e = e₀` and `e = e₀·(1−p)`); output property `p(t) ∈ [0, 1]`; equilibrium `p* = 1 − e₀/m`; reference examples for each mode verified in tests (`Math: Levins func` — 5 tests covering base model and rescue effect at specific `p` values).
- **Level 2:** equilibrium analysis, MRT method justification (A-stable implicit method suitable for stiff ODEs), convergence properties — documented in the specification. MRT solver verified against analytical solutions in `Math: MRT solver` tests (non-stiff and stiff reference problems from Chapra & Canale).

#### Computation Implementation

Computations for each task are implemented using one or a combination of the following approaches:

- **Datagrok API computation methods.** The core can use computation methods provided by the Datagrok API. Available only on the main thread — the Datagrok API is not available in web workers.

- **External libraries.** Computations are performed using third-party libraries. At the specification stage, the following is determined: which specific libraries are used, which versions, which functions/methods are applied, and a link to the library documentation (API reference, README, or guide). A documentation link is mandatory — without it, it is impossible to correctly implement and verify the calls. The order of library usage (which calls, in what sequence, with what parameters) is either defined in the specification or deferred to a separate agreement — if the usage approach is non-trivial or allows for variations. Additionally, for libraries that implement numerical methods (solvers, optimizers, fitting): the specification states which properties of the method are relevant (stability, order of accuracy, applicability class), expected accuracy for the application's use case, and the verification strategy — how the library's results will be validated (reference problems with known solutions, comparison with an alternative implementation, etc.). See section 15.3.

- **Custom methods.** Computations are implemented within the application. Each custom method is described in a separate document — a method specification. The method specification contains: mathematical formulation (formulas, equations), step-by-step algorithm, input and output data, constraints and assumptions, edge cases, and references to literature. The main application specification includes a brief method description and a link to the method specification document. Additionally, the method specification defines expected accuracy and the verification strategy: reference examples with expected outputs and their sources (manual calculation, literature, reference implementation). See section 15.3.

A single task can combine multiple approaches — for example, a custom method for data preprocessing, an external library for numerical solution, and a Datagrok computation method for postprocessing.

► **Implementation:**
- `task_primary` combines an **external library** (`diff-grok`, function `mrt`) and a **custom method** (`getEquilibrium` in `example/code/src/levins/model.ts`).
- `task_optimize` uses an **external library** (`diff-grok`, `mrt`) inside workers (`example/code/src/levins/optimize-worker.ts`) and a **custom method** (finding the maximum `p_end` in `example/code/src/levins/app.ts`, lines 421–425).
- The `mrt` call is wrapped in `createLevinsODE()` (`example/code/src/levins/model.ts`), which allows reusing the ODE specification in both the main thread and workers.

Execution environment constraint: tasks that use Datagrok API computation methods can only run on the main thread. Tasks that use only external libraries and custom methods can run on either the main thread or in a web worker.

► **Implementation:** `task_primary` runs on the main thread (synchronous, < 100 ms). `task_optimize` is distributed across workers — `diff-grok` does not depend on the Datagrok API.

General core properties:

- Does not depend on how data is obtained or how results are displayed.
- Easily testable in isolation — each task is tested separately.

► **Implementation:** `example/code/src/levins/core.ts` does not import `datagrok-api` or `ui` — only `diff-grok` and `./model`. Core tests in `example/code/src/tests/levins-api-tests.ts` test `validate`, `solve`, `validateOptimize`, `getEquilibrium` without UI.

Computation core implementation references: patterns for working with raw data and null handling (`reference/COMPUTATION-PATTERNS.md`), efficient typed array operations (`reference/ARRAY-OPERATIONS.md`).

### 1.2. Ports

Interfaces through which the core communicates with the outside world. Ports contain no implementation — only contracts. Each computational task of the core has its own set of ports:

- **Input port** — describes what parameters and types the task expects.
- **Output port** — describes the format of the task's results.
- **Progress port** — interface for reporting execution progress (percentage, stage).
- **Cancellation port** — interface for checking whether the user has requested cancellation.

Additionally, at the application level:

- **Data port** — interface for loading data from external resources.

► **Implementation:**
- **Input port `task_primary`:** interface `LevinsParams` (`example/code/src/levins/model.ts`, line 6) — `{ p0, m, e0, rescueEffect, t_start, t_end, t_step, tolerance }`.
- **Output port `task_primary`:** interface `LevinsSolution` (`example/code/src/levins/core.ts`, line 14) — `{ t: Float64Array, p: Float64Array, p_star: number }`.
- **Input port `task_optimize`:** interface `WorkerTask` (`example/code/src/levins/core.ts`, line 129) — data passed to the worker via `postMessage`.
- **Output port `task_optimize`:** interface `WorkerResult` (`example/code/src/levins/core.ts`, line 140) — `{ m_i, p_end, error? }`.
- **Progress port:** in `task_optimize` implemented via `DG.TaskBarProgressIndicator` (`example/code/src/levins/app.ts`, line 302). Progress is updated upon each completed worker.
- **Cancellation port:** `canceled` flag (`example/code/src/levins/app.ts`, line 303), checked before sending the next task to a worker.
- **Data port:** not used (the application does not load external data).

### 1.3. Adapters

Concrete implementations of ports for the Datagrok environment.

- **UI adapter** — Datagrok inputs (`ui.input.*`), buttons, custom `HTMLElement`. Converts user input into typed data for the task's input port.
- **Display adapter** — Datagrok viewers, custom `HTMLElement`, docking to the main view. Receives data from the task's output port and visualizes it.
- **Worker adapter** — wrapper for running the core in a web worker. Implements input and output ports via `postMessage`. Implementation references: worker-utils infrastructure and lifecycle (`reference/WORKER-GUIDE.md`), distribution across multiple workers (`reference/PARALLEL-EXECUTION.md`).
- **Progress adapter** — Datagrok progress bar. Implements the progress port.
- **Data adapter** — loading from a specific resource. Which resource and which mechanism — defined by the specification.

► **Implementation in `example/code/src/levins/app.ts`:**
- **UI adapter:** controls `ctrlP0`, `ctrlM`, `ctrlE0`, `ctrlRescue`, `ctrlTStart`, `ctrlTEnd`, `ctrlTStep`, `ctrlTolerance` (lines 43–99). Function `getInputs()` (line 154) converts control state to `LevinsParams`.
- **Display adapter:** `updateDataFrame()` (line 265) updates the `DG.DataFrame` and the `line chart` viewer; `updateColorCoding()` (line 188) sets conditional color coding for the `p` column.
- **Worker adapter:** `example/code/src/levins/optimize-worker.ts` — the worker receives `WorkerTask` via `onmessage`, calls `mrt(createLevinsODE(...))`, returns `WorkerResult`. The worker pool is created in `runOptimization()` (line 292).
- **Progress adapter:** `DG.TaskBarProgressIndicator.create('Optimizing m...')` (line 302), updated via `pi.update(...)` (line 365).
- **Data adapter:** not used.

### 1.4. Coordinator (Application Service)

The coordinator connects adapters and the core. It:

- Listens for input changes via the UI adapter.
- Manages reactivity: cascading dependencies between inputs, range updates, defaults, availability — based on the specification.
- Triggers validation and computations through the corresponding ports.
- Manages control state (enabled/disabled) during computations.
- Manages computation blocking: can suspend reactive pipeline execution during batch input updates (see section 8.4).
- Passes results to the display adapter.
- Manages resource lifecycle (subscriptions, workers).

► **Implementation:** function `levinsMetapopulationApp()` in `example/code/src/levins/app.ts` fulfills all coordinator roles:
- Listens to `onValueChanged` via input callbacks (lines 43–99).
- Reactivity: `updateRhoBadge()` (line 124), `updateRescueLabel()` (line 136).
- Validation and computation: `runPrimary()` (line 226).
- Blocking: `computationsBlocked` flag (line 19), used during input formatting (lines 102–110), during reset (lines 506–518), and after optimization (lines 428–437).
- Lifecycle: `subs[]` (line 21), `activeWorkers[]` (line 22), cleanup in `onViewRemoved` (line 561).

### 1.5. Independence Principle

Input behavior does not depend on the core's computational part. The UI adapter and reactivity between inputs form a standalone layer managed by the coordinator based on the specification. The core receives a ready, validated set of parameters.

► **Implementation:** the function `runPrimary()` first calls `validate(inputs)` from the core, and only if `errors.size === 0` passes the data to `solve(inputs)`. The core (`core.ts`) is unaware of the inputs' existence — it accepts a plain `LevinsParams`.

## 2. Main View

The application has a main view. If a scientific application is based on a table, the main view should be `DG.TableView`.

► **Implementation:** `example/code/src/levins/app.ts`, line 33: `const view = grok.shell.addTableView(df)`.

## 3. Controls (Inputs)

Users set input parameters for computations through controls. Control types:

- **Standard Datagrok inputs** — created via `ui.input.*`.
- **Datagrok buttons** — perform a specified action when clicked.
- **Custom HTMLElement** — for example, a `div` element with specific styles that performs an action when clicked. Each custom element is described in a separate document — a UI component specification. The UI component specification contains: visual description (sketch/mockup), states (normal, hover, disabled, active), events (what it emits on interaction), styles (CSS classes), accessibility (tooltips, aria). The main application specification includes a brief element description and a link to the UI component specification document.

► **Implementation:**
- **Standard inputs:** 7 numeric + 1 toggle in `example/code/src/levins/app.ts`, lines 43–99 (`ctrlP0`, `ctrlM`, `ctrlE0`, `ctrlRescue`, `ctrlTStart`, `ctrlTEnd`, `ctrlTStep`, `ctrlTolerance`).
- **Buttons:** `optimizeBtn = ui.iconFA('search', ...)` (line 503), `resetBtn = ui.iconFA('undo', ...)` (line 505).
- **Custom HTMLElement:** `rhoBadge = ui.div([], 'd4-tag levins-rho-badge')` (line 37) — a rho = e0/m indicator with dynamic color switching via CSS classes.

### 3.1. Input Options

When creating standard Datagrok inputs (`ui.input.*`), the specification defines options for each input. These options are passed during input creation and control its behavior:

- **Min** — minimum allowed value. For numeric inputs, sets the lower bound.
- **Max** — maximum allowed value. For numeric inputs, sets the upper bound.
- **Format** — value display format (e.g., `0.000` for three decimal places, `0.0` for one, `0.##E+0` for scientific notation). The format determines how the value is displayed in the input.

Input options differ from validation (section 7): options set basic control-level constraints (range, format), while validation checks complex conditions across combinations of multiple input values.

► **Implementation:** min/max are set during input creation (e.g., `ctrlP0`: `min: 0.001, max: 1`). Formats are set in a separate block (lines 102–109), wrapped in `computationsBlocked = true/false` — so that format assignment does not trigger a side-effect recomputation.

### 3.2. Control Classification

Controls are classified by ownership:

- **Main view controls** — inputs of the primary pipeline, placed in the main application interface.
- **Secondary task triggers** — buttons or icons that launch secondary pipelines (see section 8.2).
- **Secondary task controls** — inputs placed in the secondary task's own UI (e.g., in a dialog).

► **Implementation:**
- **Main view controls:** `ctrlP0`…`ctrlTolerance` — in the left panel form.
- **Secondary task trigger:** `optimizeBtn` (`ui.iconFA('search')`, line 503) → opens the optimization dialog.
- **Secondary task controls:** `dlgMMin`, `dlgMMax` — dialog inputs in `showOptimizeDialog()` (lines 448–499).

## 4. Result Display Elements

Computation results are displayed using:

- **Standard Datagrok viewers.**
- **Custom HTMLElement.** Each custom display element is described in a separate UI component specification (similar to custom controls, see section 3).

By default, these elements are docked to the main view.

► **Implementation:**
- **Viewer:** `line chart` — `view.addViewer('Line chart', {...})` (line 546), docked to the right of the grid.
- **Custom element:** `rhoBadge` — displays the current rho value with color indication (green/red).
- **Color coding of column `p`:** `updateColorCoding()` (line 188) sets conditional colors (green — persistence zone, red — extinction threat) with dynamic threshold recalculation `e0/m`.
- **Column `p` header tooltip:** `setupGridTooltip()` (line 202) via `view.grid.onCellTooltip`.

## 5. Layout and UI Element Placement

The placement of controls and display elements (panels, ribbon, toolbar, side panels, docking area) is defined by the application specification.

► **Implementation (`example/code/src/levins/app.ts`):**
- **Left panel:** `ui.form` with groups via `ui.h2` (lines 523–540), docked as `DG.DOCK_TYPE.LEFT`, ratio `0.2` (line 543).
- **Toolbar (ribbon):** `view.setRibbonPanels([[optimizeBtn, resetBtn]])` (line 520).
- **Main area:** `DG.TableView` (grid) — default.
- **Right area:** `line chart`, docked `DG.DOCK_TYPE.RIGHT` relative to the grid, ratio `0.5` (lines 552–554).

## 5.1. Styles

All visual styles of the application are placed in a separate CSS file (`css/<app-name>.css`). Inline styles in TypeScript code are not allowed — CSS classes are used instead.

- **Static styles** — element styling that does not change during operation. Set via CSS class when creating the element.
- **Dynamic styles** — styles that depend on application state (e.g., indicator color, button activity). Implemented via CSS class toggling (`classList.toggle`, `classList.add/remove`), not via direct `element.style.*` assignment.

The CSS file is imported via ES import (`import '../css/<app-name>.css'`). Webpack with `style-loader` + `css-loader` injects styles into the DOM when the bundle loads.

► **Implementation:**
- CSS file: `example/code/css/levins.css` — contains classes `.levins-rho-badge`, `.levins-rho-badge--persists`, `.levins-rho-badge--extinct`, `.levins-btn--disabled`.
- Import: `import '../../css/levins.css'` (`example/code/src/levins/app.ts`, line 12).
- **Static styles:** `rhoBadge` is created with classes `'d4-tag levins-rho-badge'` (line 37).
- **Dynamic styles:** `rhoBadge.classList.toggle('levins-rho-badge--persists', persists)` (line 130); `optimizeBtn.classList.toggle('levins-btn--disabled', !enabled)` (line 289).

## 6. User Feedback

### 6.1. Control Tooltips

- For standard Datagrok inputs, tooltips are defined at input creation time via the `tooltipText` property.
- For elements that are not Datagrok inputs, tooltips are bound via `ui.tooltip.bind`.

► **Implementation:**
- All inputs have `tooltipText` (e.g., `ctrlP0`: `tooltipText: 'Fraction of patches occupied at t=0...'`, line 46).
- `rhoBadge`: `ui.tooltip.bind(rhoBadge, '...')` (line 38).
- Buttons `optimizeBtn` and `resetBtn`: tooltip is passed as the third argument of `ui.iconFA` (lines 503, 518).

### 6.2. Validators as Feedback

Validators are added to standard Datagrok inputs — they display inline hints about the validity of the current value.

► **Implementation:** function `addValidators()` (`example/code/src/levins/app.ts`, line 168) adds a validator to each input. The validator calls `validate(getInputs())` from the core and returns an error for the specific `InputId`.

### 6.3. Progress Bar

During long computations, the standard Datagrok progress bar indicator is displayed with cancellation support.

► **Implementation:** `task_optimize` uses `DG.TaskBarProgressIndicator.create('Optimizing m...')` (line 302), updated via `pi.update(percent, label)` upon each completed worker (line 365).

## 7. Validation

Validation is performed through the Datagrok validator mechanism. Each computational task has its own validation rules.

### 7.1. Complex Validation

Validation is complex: the entire set of task inputs is analyzed as a whole. If a combination of values is invalid, a `Map` with error specifications is returned. This specification is then used by the Datagrok validator mechanism to display errors on the corresponding inputs.

► **Implementation:**
- `task_primary`: function `validate(inputs: LevinsParams): ValidationErrors` (`example/code/src/levins/core.ts`, line 40) — returns `Map<InputId, string>`. Rules val_01…val_09 check both individual values and combinations (e.g., val_05: `m <= e0` when `rescueEffect = false`; val_08: `t_step >= t_end - t_start`).
- `task_optimize`: function `validateOptimize(opt, e0, rescueEffect)` (`example/code/src/levins/core.ts`, line 106) — returns `{ errors: Map<OptInputId, string>, warning: string | null }`. Rules opt_val_01…opt_val_04.

### 7.2. Validation Order

The order of determining input set validity is defined by the specification for each task.

► **Implementation:** in `validate()`, combinatorial rules are checked only after basic rules pass:
- val_05 is checked only if val_03 and val_04 passed (line 59: `if (!errors.has('ctrl_m') && !errors.has('ctrl_e0') && ...)`).
- val_08 is checked only if val_06 and val_07 passed (line 73: `if (!errors.has('ctrl_t_end') && !errors.has('ctrl_t_step') && ...)`).

## 8. Main Pipeline

Each computational task of the core has its own pipeline. All pipelines are orchestrated by the coordinator (see section 1.4).

### 8.1. Primary Pipeline

The primary pipeline is bound to the main application controls and operates reactively.

**Parameter input.** The user sets input values through controls. The UI adapter converts the input into typed data for the task's input port.

**Validation.** The input set is validated (see section 7). Validation is performed by the core through the input port.

**Computation.** Validated inputs are passed to the core task. Execution characteristics (synchronicity, environment, parallelization) are defined by the task specification (see section 1.1).

**Result display.** The core returns results through the task's output port. The coordinator passes them to the display adapter (see section 4).

► **Implementation:** function `runPrimary()` (`example/code/src/levins/app.ts`, line 226):
1. Checks `computationsBlocked` (line 227).
2. Collects inputs: `getInputs()` (line 230).
3. Validates: `validate(inputs)` (line 231).
4. On errors — visually marks invalid inputs (`d4-invalid`) and calls `clearResults()` (lines 237–245).
5. On success — `solve(inputs)` (line 248), then `updateDataFrame(result)` + `updateColorCoding()` (lines 249–250).

Reactive trigger: each input has `onValueChanged: () => debouncedRun()` (debounce 50 ms, line 258).

### 8.2. Secondary Pipelines

A secondary pipeline is bound to an explicit user action (button, icon, menu item) and runs on demand.

**Trigger.** The user initiates an action (e.g., clicks the "Sensitivity Analysis" button).

**Custom UI.** A secondary pipeline may have its own parameter input interface — for example, a Datagrok dialog with its own set of inputs, validation, and tooltips. This UI is defined by the secondary task specification.

**Validation.** Secondary task parameters are validated independently — using their own complex validation rules.

**Computation.** The core task is executed. The secondary task may use primary task results as part of its input data.

**Result display.** Secondary task results are displayed with their own elements — these can be additional viewers docked to the main view, dialog content, or a separate window. Defined by the specification.

**Feedback to primary pipeline.** The secondary task may return values that are substituted into primary pipeline controls. In this case, the computation blocking and batch update mechanism is used (see section 8.4).

► **Implementation of `task_optimize`:**
1. **Trigger:** `optimizeBtn` → `showOptimizeDialog()` (line 503).
2. **Custom UI:** dialog `ui.dialog('Find optimal m')` with inputs `dlgMMin`, `dlgMMax` (lines 448–499), with its own validators and tooltips.
3. **Validation:** `validateDialog()` (line 463) calls `validateOptimize()` from the core.
4. **Computation:** `runOptimization(mMin, mMax)` (line 292) — creates a worker pool, distributes 10,000 tasks, collects results.
5. **Display:** `grok.shell.info(...)` with the found optimum (line 433).
6. **Feedback:** `ctrlM.value = best.m_i` via batch update (lines 428–433).

### 8.3. Common Pipeline Aspects

The following aspects are defined by the specification for each pipeline independently:

#### Control Behavior During Computations

Controls may become disabled during computations — which ones and for which pipeline is defined by the specification.

► **Implementation:** `task_primary` does not block controls (< 100 ms). `task_optimize` blocks only `optimizeBtn` via `setOptimizeBtnEnabled(false)` (lines 288–289, 300).

#### Progress and Cancellation

During long computations, the standard Datagrok progress bar is displayed with cancellation support.

► **Implementation:** `task_primary` — no progress bar. `task_optimize` — `DG.TaskBarProgressIndicator` with completion percentage and cancellation via the `canceled` flag (line 303).

#### Computation Error Handling

The course of action upon computation failure is defined by the specification for each task.

► **Implementation:**
- `task_primary`: `catch` → `clearResults()` + `grok.shell.error(msg)` (lines 251–255).
- `task_optimize`: partial errors (some workers failed) → `grok.shell.warning(...)` (line 413); complete failure → `grok.shell.error(...)` (line 416); error writing to control → `grok.shell.warning(...)` with instructions to enter manually (line 436).

### 8.4. Computation Blocking and Batch Input Updates

In certain scenarios, the result of a secondary task must be substituted into primary pipeline controls. In this case, each individual input change should not trigger a reactive recomputation — the computation should happen once, after all values have been substituted.

The coordinator supports a computation blocking mode:

**Block request.** Before batch update, the coordinator suspends reactive execution of specified pipelines. Which pipelines are blocked is defined by the specification.

**Batch update.** Values are substituted into controls. Cascading dependencies between inputs (reactivity, section 9) can be handled in one of two modes — defined by the specification:

- Reactivity between inputs works as usual, but computations are not triggered.
- Reactivity between inputs is also suspended until the batch update completes.

**Unblock.** After all values are substituted, the coordinator removes the block. Full input set validation occurs, then computation runs, then results are displayed — the standard pipeline (section 8.1).

Example: the user launches parameter optimization (secondary task). Upon optimization completion, the coordinator blocks the primary pipeline, substitutes the found parameter values into all controls, removes the block — the primary pipeline runs once for the complete set of optimal parameters.

► **Implementation:** the `computationsBlocked` flag (`example/code/src/levins/app.ts`, line 19) is used in three scenarios:
1. **Format initialization** (lines 102–110): blocking prevents side-effect recomputations during `format` assignment.
2. **Reset** (lines 506–518): `computationsBlocked = true` → reset all values → `computationsBlocked = false` → `runPrimary()`.
3. **Optimization result** (lines 428–437): `computationsBlocked = true` → `ctrlM.value = best.m_i` → `computationsBlocked = false` → `runPrimary()`.

Blocking check: `runPrimary()` first checks `if (computationsBlocked) return` (line 227).

## 9. Reactivity and Dependencies Between Inputs

Dependencies between inputs (cascading updates of ranges, defaults, availability) are defined by the application specification. Reactivity is managed by the coordinator (see section 1.4) and operates entirely at the UI adapter level — independent of the core's computational part.

Reactivity can be temporarily suspended by the coordinator in batch input update mode (see section 8.4).

► **Implementation:**
- **`ctrl_m` / `ctrl_e0` → `rhoBadge`:** `updateRhoBadge()` (line 124) is called from `onValueChanged` of both inputs. Recalculates `rho = e0/m`, updates text and CSS class.
- **`ctrl_rescue` → `ctrl_e0` (label + tooltip):** `updateRescueLabel()` (line 136) switches caption and tooltip based on the toggle state.
- **`ctrl_t_start` / `ctrl_t_end` → ranges:** `updateArgRanges()` (line 149) — scaffold for range updates; actual checking via complex validation (val_06, val_07, val_08).
- **Debounce:** numeric inputs use `debouncedRun()` (debounce 50 ms, line 258), toggle `ctrl_rescue` calls `runPrimary()` directly (without debounce).

## 10. Data Lifecycle

### 10.1. Data Input

The primary approach is manual entry through application controls.

► **Implementation:** all model parameters are entered by the user through the form in the left panel. The initial state uses default values from `DEFAULTS` (`example/code/src/levins/core.ts`, line 27). `task_primary` runs automatically on initialization: `solve(DEFAULTS)` (line 26).

### 10.2. Loading from a Resource

Via buttons or icons — loading data from an external resource. Which specific resource and loading mechanism is defined by the application specification.

► **Implementation:** not used in this application.

## 11. Error Handling Beyond Computations

The strategy for handling data loading errors, network errors, invalid input files, and incorrect application state is defined by the application specification.

► **Implementation:** see specification (`example/levins-metapopulation-spec.md`, section 11). Examples:
- Worker creation error: `reject(new Error('Failed to start parallel computations...'))` (line 380).
- Partial worker errors: `errorCount` counting and `grok.shell.warning(...)` (lines 305, 412–413).
- Error writing to control: fallback with `grok.shell.warning(...)` (line 436).

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

All Datagrok event subscriptions (`onValueChanged`, `onAfterDraw`, etc.) must be collected and unsubscribed when the application closes via `sub.unsubscribe()`.

► **Implementation:** array `subs` (`example/code/src/levins/app.ts`, line 21) collects subscriptions. On close — `for (const sub of subs) sub.unsubscribe()` (line 564).

### 12.2. Worker Termination

When the application closes, all web workers must be properly terminated.

► **Implementation:** array `activeWorkers` (line 22), function `terminateWorkers()` (line 440) calls `w.terminate()` for each worker. Called on close (line 563) and after optimization completes (line 405).

## 13. Application Closure

When the application closes, the coordinator performs:

- All event subscriptions are unsubscribed — for both primary and secondary pipelines (see section 12.1).
- All web workers are terminated — including secondary task workers (see section 12.2).
- All associated UI elements are closed (including open secondary task dialogs).
- Pending requests are cancelled.

► **Implementation:** handler `grok.events.onViewRemoved.subscribe(...)` (`example/code/src/levins/app.ts`, lines 561–568):
1. `terminateWorkers()` — terminates all active workers.
2. `for (const sub of subs) sub.unsubscribe()` — unsubscribes subscriptions.
3. `clearTimeout(debounceTimer)` — cancels pending debounce.

## 14. Accessibility and UX

Keyboard shortcuts, context menus, undo/redo, and other UX elements are defined by the application specification.

► **Implementation:** in the current version of the application, keyboard shortcuts, context menus, and undo/redo are not implemented (see specification, section 12).

## 15. Testing

### 15.1. Computational Part (Core)

Core correctness verification: unit tests for each computational task separately. The core is tested in isolation — without UI and adapters.

► **Implementation:** tests are split across two files:

**`example/code/src/tests/levins-api-tests.ts`** — 2 categories:
- **API: Validation** — 24 tests: val_01…val_09 + boundary values + dependency order + multiple errors + defaults.
- **API: Optimization Validation** — 7 tests: opt_val_01…opt_val_04 + valid input data.

**`example/code/src/tests/levins-math-tests.ts`** — 4 categories:
- **Math: MRT solver** — 3 tests: non-stiff 1D, stiff 1D, stiff 2D (van der Pol) — verifying the `mrt` solver against analytical/reference solutions.
- **Math: Levins func** — 5 tests: ODE right-hand side correctness for the base model and rescue effect at specific `p` values.
- **Math: Equilibrium** — 4 tests: `getEquilibrium` for the base model and rescue effect.
- **Math: Solve output properties** — 8 tests: output invariants of `solve()` — bounds `p ∈ [0, 1]`, initial conditions, convergence to `p*`, monotonicity in `m`.

Tests are run via `grok test` (entry point: `example/code/src/package-test.ts`).

### 15.2. Inputs

Input verification for each task: all cases including edge cases (boundary values, invalid combinations, empty values, extreme values).

► **Implementation:** validation tests cover:
- Boundary values: `p0 = 0.001` (lower allowed bound), `p0 = 1` (upper).
- Invalid combinations: `m <= e0` without rescue, `t_step >= t_end - t_start`.
- Dependencies between rules: val_05 is skipped when val_03 fails, val_08 is skipped when val_06/val_07 fail.
- Multiple errors: simultaneously invalid `p0`, `m`, `e0`, `t_step`, `tolerance`.

### 15.3. Mathematical Verification

Verification that the implemented computation matches the model definition (see section 1.1, "Computation Formulas and Model"). Test categories correspond to the two levels of the model definition. Verification criteria — reference examples, expected accuracies, output property constraints, reference problems for the numerical method — are defined by the model specification, not invented during test writing. Tests implement what is specified; the specification is the source of truth.

#### Level 1 verification (required)

**Formula/equation verification.** The implemented transformation is checked at control points with manually computed expected values. For each computational path (mode, branch, regime), at least one test substitutes concrete inputs and compares the output against a hand-calculated result.

**Output property verification.** Constraints declared in the model definition (bounds, monotonicity, conservation laws, limiting cases) are checked on actual computation results. These tests do not compare against a specific expected value — they verify that the result satisfies a declared invariant.

► **Implementation:**
- **Formula verification:** `Math: Levins func` — 5 tests. Each test substitutes specific `(m, e₀, p)` into the ODE right-hand side and compares `dp/dt` against a hand-calculated value (e.g., `func_01`: `dp/dt = 0.5·0.5·0.5 − 0.2·0.5 = 0.025`). Both computational paths are covered: base model (3 tests) and rescue effect (2 tests).
- **Equilibrium verification:** `Math: Equilibrium` — 4 tests. `getEquilibrium` is checked against the analytical formula `p* = 1 − e₀/m` for the base model, and `NaN` for the rescue effect (no closed-form equilibrium).
- **Output properties:** `Math: Solve output properties` — 8 tests. Verifies invariants declared in the specification: `p(t) ∈ [0, 1]` (solve_02, solve_06), `p(0) = p0` (solve_03, solve_08), convergence to `p*` (solve_05), monotonicity in `m` (solve_07), non-empty output arrays (solve_01), `t[0] = t_start` (solve_04).

#### Level 2 verification (for full formalization)

**Numerical method verification.** The solver (or library) is tested on reference problems with known analytical solutions. The test verifies that the numerical error stays within the expected tolerance. Reference problems should cover the solver's applicability range (e.g., non-stiff and stiff problems for an ODE solver). Each reference problem must cite its source (textbook, paper, test suite).

**Convergence verification.** Solving the same problem with decreasing step size or tolerance produces solutions that converge. The test compares solutions at two different precision levels and verifies that the discrepancy decreases.

**Asymptotic/equilibrium behavior.** The numerical solution on a sufficiently long interval approaches the analytically predicted equilibrium or asymptote.

► **Implementation:**
- **Numerical method:** `Math: MRT solver` — 3 tests. Non-stiff 1D and stiff 1D problems verified against analytical solutions (Chapra & Canale, pp. 736, 767). Stiff 2D van der Pol (µ=1000) verified for solver stability (reference: VDPOL test set). Tolerance threshold: max absolute error < 0.1.
- **Convergence:** not yet covered. Candidate: solve the Levins ODE with `tolerance = 1e-5` and `tolerance = 1e-9`, verify that the discrepancy between solutions decreases.
- **Asymptotic behavior:** not yet covered. Candidate: verify that `p(t_end)` approaches `p* = 1 − e₀/m` for sufficiently large `t_end`.
