# Guide: Building Interactive Scientific Web Applications on Datagrok

> **End-to-end example:** application **Lotka-Volterra Predator-Prey Simulation** (see `./examples/lotka-volterra-spec/` directory).
> Application specification: `./examples/lotka-volterra-spec/lotka-volterra-spec.md`.

## 1. General Architecture

An application is a Datagrok package function registered with the `//tags: app` comment.

► **Implementation:** `./examples/lotka-volterra-spec/code/src/package.ts` — function `lotkaVolterraSimulation`, registered via `//name: Lotka-Volterra Simulation` / `//tags: app`. Delegates logic to `./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts` → `lotkaVolterraApp()`.

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

► **Implementation:** `./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts` — core with three tasks:
- **`task_primary`** (synchronous, main thread) — `solve(inputs)` solves the Lotka-Volterra 2D ODE system and returns trajectories `x(t)`, `y(t)`, equilibrium point, and summary statistics.
- **`task_optimize`** (asynchronous, web workers) — brute-force grid search over 4 coefficients (α, β, δ, γ) with 11 steps each (14,641 points total), executed in parallel in a worker pool.
- **`task_equilibrium`** (asynchronous, single web worker) — varies a selected parameter across its range (100,000 steps), computing the equilibrium point for each value. Opens results in a new TableView.

Examples of computational tasks within a single application:

- **Primary task** — solving the ODE based on model parameters. Triggered reactively when inputs change.
- **Secondary task** — model sensitivity analysis. Triggered by an explicit user action (button), has its own set of parameters.

Tasks can be independent or linked (the secondary task uses results from the primary one).

► **Implementation:** `task_optimize` does not depend on the results of `task_primary`, but uses the current input values of the primary pipeline (snapshot of `x0`, `y0`, `T`). After optimization completes, the result is written to `ctrl_alpha`, `ctrl_beta`, `ctrl_delta`, `ctrl_gamma` controls via batch update, which triggers a single run of `task_primary`. `task_equilibrium` is fully independent — it opens results in a new view.

#### Computation Formulas and Model

Each computational task is based on a defined transformation of inputs into outputs — from a single formula to a complex system of equations. In this guide, the term "model" refers to any such definition: individual formulas, chains of transformations, ODE/PDE systems, optimization problems, or statistical procedures. The model exists independently of its implementation (custom code, external library, or platform API). The model definition is the primary source of verification criteria: what cannot be defined cannot be verified.

The model definition has two levels.

**Level 1 — required minimum (before implementation):**

- **Variables:** name, meaning, units of measurement, valid domain (e.g., `x ≥ 0`, `α > 0`).
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
- **Level 1** is defined in the application specification (`./examples/lotka-volterra-spec/lotka-volterra-spec.md`): variables `x, y, α, β, δ, γ` with units and domains; ODE system `dx/dt = αx − βxy`, `dy/dt = δxy − γy`; output properties `x(t) ≥ 0`, `y(t) ≥ 0`, `x(0) = x₀`, `y(0) = y₀`; equilibrium `(x*, y*) = (γ/δ, α/β)`; 5 reference examples verified in tests (`Math: LV func` — 5 tests covering different parameter sets and equilibrium point).
- **Level 2:** equilibrium analysis (periodic orbits around non-trivial equilibrium, conserved Hamiltonian `H = δx + βy − γ ln(x) − α ln(y)`), MRT method justification — documented in the specification. MRT solver verified against analytical solutions in `Math: MRT solver` tests (non-stiff and stiff reference problems from Chapra & Canale).

#### Computation Implementation

Computations for each task are implemented using one or a combination of the following approaches:

- **Datagrok API computation methods.** The core can use computation methods provided by the Datagrok API. Available only on the main thread — the Datagrok API is not available in web workers.

- **External libraries.** Computations are performed using third-party libraries. At the specification stage, the following is determined: which specific libraries are used, which versions, which functions/methods are applied, and a link to the library documentation (API reference, README, or guide). A documentation link is mandatory — without it, it is impossible to correctly implement and verify the calls. The order of library usage (which calls, in what sequence, with what parameters) is either defined in the specification or deferred to a separate agreement — if the usage approach is non-trivial or allows for variations. Additionally, for libraries that implement numerical methods (solvers, optimizers, fitting): the specification states which properties of the method are relevant (stability, order of accuracy, applicability class), expected accuracy for the application's use case, and the verification strategy — how the library's results will be validated (reference problems with known solutions, comparison with an alternative implementation, etc.). See section 15.3.

- **Custom methods.** Computations are implemented within the application. Each custom method is described in a separate document — a method specification. The method specification contains: mathematical formulation (formulas, equations), step-by-step algorithm, input and output data, constraints and assumptions, edge cases, and references to literature. The main application specification includes a brief method description and a link to the method specification document. Additionally, the method specification defines expected accuracy and the verification strategy: reference examples with expected outputs and their sources (manual calculation, literature, reference implementation). See section 15.3.

A single task can combine multiple approaches — for example, a custom method for data preprocessing, an external library for numerical solution, and a Datagrok computation method for postprocessing.

► **Implementation:**
- `task_primary` combines an **external library** (`diff-grok`, function `mrt`) and a **custom method** (`getEquilibrium` in `./examples/lotka-volterra-spec/code/src/lotka-volterra/model.ts`, plus summary stats computation in `solve()` in `core.ts`).
- `task_optimize` uses an **external library** (`diff-grok`, `mrt`) inside workers (`./examples/lotka-volterra-spec/code/src/lotka-volterra/workers/optimize-worker.ts`) and a **custom method** (grid generation via `linspace()` and finding the maximum `maxPrey` in `./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`).
- `task_equilibrium` uses a **custom method** in a worker (`./examples/lotka-volterra-spec/code/src/lotka-volterra/workers/equilibrium-worker.ts`) — computing `x* = γ/δ`, `y* = α/β` for each parameter value.
- The `mrt` call is wrapped in `createLotkaVolterraODE()` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/model.ts`), which allows reusing the ODE specification in both the main thread and workers.

Execution environment constraint: tasks that use Datagrok API computation methods can only run on the main thread. Tasks that use only external libraries and custom methods can run on either the main thread or in a web worker.

► **Implementation:** `task_primary` runs on the main thread (synchronous, < 100 ms). `task_optimize` is distributed across workers — `diff-grok` does not depend on the Datagrok API. `task_equilibrium` runs in a single worker.

General core properties:

- Does not depend on how data is obtained or how results are displayed.
- Easily testable in isolation — each task is tested separately.

► **Implementation:** `./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts` does not import `datagrok-api` or `ui` — only `diff-grok` and `./model`. Core tests in `./examples/lotka-volterra-spec/code/src/tests/lotka-volterra-api-tests.ts` test `validate`, `solve`, `getEquilibrium` without UI.

Computation core implementation references: patterns for working with raw data and null handling (`./reference/COMPUTATION-PATTERNS.md`), efficient typed array operations (`./reference/ARRAY-OPERATIONS.md`).

### 1.2. Ports

Interfaces through which the core communicates with the outside world. Ports contain no implementation — only contracts. Each computational task of the core has its own set of ports:

- **Input port** — describes what parameters and types the task expects.
- **Output port** — describes the format of the task's results.
- **Progress port** — interface for reporting execution progress (percentage, stage).
- **Cancellation port** — interface for checking whether the user has requested cancellation.

Additionally, at the application level:

- **Data port** — interface for loading data from external resources.

► **Implementation:**
- **Input port `task_primary`:** interface `LotkaVolterraParams` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/model.ts`, line 6) — `{ alpha, beta, delta, gamma, x0, y0, T }`.
- **Output port `task_primary`:** interface `LotkaVolterraSolution` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts`, line 14) — `{ t, x, y, xStar, yStar, maxPrey, maxPredators, stepCount }`.
- **Input port `task_optimize`:** interface `WorkerTask` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts`, line 115) — `{ alpha, beta, delta, gamma, x0, y0, T }` passed to worker via `postMessage`.
- **Output port `task_optimize`:** interface `WorkerResult` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts`, line 125) — `{ alpha, beta, delta, gamma, maxPrey, error? }`.
- **Input port `task_equilibrium`:** interface `EquilibriumTask` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts`, line 140) — `{ paramName, paramMin, paramMax, steps, baseParams }`.
- **Output port `task_equilibrium`:** interface `EquilibriumResult` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts`, line 148) — `{ paramValues[], xStar[], yStar[], error? }`.
- **Progress port:** `DG.TaskBarProgressIndicator` (indeterminate, cancelable) for both `task_optimize` and `task_equilibrium`.
- **Cancellation port:** `pi.onCanceled` → terminates workers and resolves promises.
- **Data port:** not used (the application does not load external data).

### 1.3. Adapters

Concrete implementations of ports for the Datagrok environment.

- **UI adapter** — Datagrok inputs (`ui.input.*`), buttons, custom `HTMLElement`. Converts user input into typed data for the task's input port.
- **Display adapter** — Datagrok viewers, custom `HTMLElement`, docking to the main view. Receives data from the task's output port and visualizes it.
- **Worker adapter** — wrapper for running the core in a web worker. Implements input and output ports via `postMessage`. Worker source files must be placed in a dedicated `workers/` subdirectory (e.g., `src/<app-name>/workers/optimize-worker.ts`) to separate them from the main application code. Implementation references: worker-utils infrastructure and lifecycle (`./reference/WORKER-GUIDE.md`), distribution across multiple workers (`./reference/PARALLEL-EXECUTION.md`).
- **Progress adapter** — Datagrok progress bar. Implements the progress port.
- **Data adapter** — loading from a specific resource. Which resource and which mechanism — defined by the specification.

► **Implementation in `./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`:**
- **UI adapter:** controls `ctrlAlpha`, `ctrlBeta`, `ctrlDelta`, `ctrlGamma`, `ctrlX0`, `ctrlY0`, `ctrlT` (lines 76–124). Function `getInputs()` (line 149) converts control state to `LotkaVolterraParams`.
- **Display adapter:** `updateDataFrame()` (line 226) updates the `DG.DataFrame` and all viewers (line chart, scatter plot, grid); `updateStatsPanel()` (line 65) updates equilibrium and max labels; `applyColorCoding()` (line 221) sets linear color coding on `Prey` and `Predators` columns.
- **Worker adapter:** `./examples/lotka-volterra-spec/code/src/lotka-volterra/workers/optimize-worker.ts` — receives `WorkerTask[]` via `onmessage`, calls `mrt(createLotkaVolterraODE(...))`, returns `WorkerResult[]`. `./examples/lotka-volterra-spec/code/src/lotka-volterra/workers/equilibrium-worker.ts` — receives `EquilibriumTask`, computes equilibrium for each parameter value, returns `EquilibriumResult`. The worker pools are created in `runOptimization()` (line 261) and `runEquilibriumExploration()` (line 410).
- **Progress adapter:** `DG.TaskBarProgressIndicator.create(...)` — indeterminate, cancelable — used for both `task_optimize` (line 297) and `task_equilibrium` (line 423).
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

► **Implementation:** function `lotkaVolterraApp()` in `./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts` fulfills all coordinator roles:
- Listens to `onValueChanged` via input callbacks (lines 76–124).
- Reactivity: all 7 controls trigger `debouncedRun()` — no cascading dependencies between controls in this application.
- Validation and computation: `runPrimary()` (line 182).
- Blocking: `computationsBlocked` flag (line 22), used during input formatting (lines 127–135), during reset (lines 505–516), and after optimization (lines 376–387).
- Lifecycle: `subs[]` (line 24), `activeWorkers[]` (line 25), cleanup in `onViewRemoved` (line 605).

### 1.5. Independence Principle

Input behavior does not depend on the core's computational part. The UI adapter and reactivity between inputs form a standalone layer managed by the coordinator based on the specification. The core receives a ready, validated set of parameters.

► **Implementation:** the function `runPrimary()` first calls `validate(inputs)` from the core, and only if `errors.size === 0` passes the data to `solve(inputs)`. The core (`core.ts`) is unaware of the inputs' existence — it accepts a plain `LotkaVolterraParams`.

## 2. Main View

The application has a main view. If a scientific application is based on a table, the main view should be `DG.TableView`.

► **Implementation:** `./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, line 40: `const view = grok.shell.addTableView(df)`.

## 3. Controls (Inputs)

Users set input parameters for computations through controls. Control types:

- **Standard Datagrok inputs** — created via `ui.input.*`.
- **Datagrok buttons** — perform a specified action when clicked.
- **Custom HTMLElement** — for example, a `div` element with specific styles that performs an action when clicked. Each custom element is described in a separate document — a UI component specification. The UI component specification contains: visual description (sketch/mockup), states (normal, hover, disabled, active), events (what it emits on interaction), styles (CSS classes), accessibility (tooltips, aria). The main application specification includes a brief element description and a link to the UI component specification document.

► **Implementation:**
- **Standard inputs:** 7 numeric inputs in `./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, lines 76–124 (`ctrlAlpha`, `ctrlBeta`, `ctrlDelta`, `ctrlGamma`, `ctrlX0`, `ctrlY0`, `ctrlT`).
- **Buttons:** `optimizeBtn = ui.bigButton('Optimize', ...)` (line 254), `resetBtn = ui.iconFA('undo', ...)` (line 505), `analyzeBtn = ui.iconFA('chart-line', ...)` (line 518), `helpBtn = ui.icons.help(...)` (line 521).
- **Custom HTMLElement:** `statsPanel` — a `ui.divV` panel with `ui.label` components (lines 46–63) showing equilibrium values and max population statistics.

### 3.1. Input Options

When creating standard Datagrok inputs (`ui.input.*`), the specification defines options for each input. These options are passed during input creation and control its behavior:

- **Min** — minimum allowed value. For numeric inputs, sets the lower bound.
- **Max** — maximum allowed value. For numeric inputs, sets the upper bound.
- **Format** — value display format (e.g., `0.000` for three decimal places, `0.0` for one, `0.##E+0` for scientific notation). The format determines how the value is displayed in the input.

Input options differ from validation (section 7): options set basic control-level constraints (range, format), while validation checks complex conditions across combinations of multiple input values.

► **Implementation:** min/max are set during input creation from `RANGES` (e.g., `ctrlAlpha`: `min: RANGES.alpha.min, max: RANGES.alpha.max`). Formats are set in a separate block (lines 127–134), wrapped in `computationsBlocked = true/false` — so that format assignment does not trigger a side-effect recomputation.

### 3.2. Control Classification

Controls are classified by ownership:

- **Main view controls** — inputs of the primary pipeline, placed in the main application interface.
- **Secondary task triggers** — buttons or icons that launch secondary pipelines (see section 8.2).
- **Secondary task controls** — inputs placed in the secondary task's own UI (e.g., in a dialog).

► **Implementation:**
- **Main view controls:** `ctrlAlpha`…`ctrlT` — in the left panel form.
- **Secondary task triggers:** `optimizeBtn` (`ui.bigButton('Optimize')`, line 254) → launches optimization directly; `analyzeBtn` (`ui.iconFA('chart-line')`, line 518) → opens the equilibrium exploration dialog.
- **Secondary task controls:** `choiceInput` (`ui.input.choice`) — dialog input in `showEquilibriumDialog()` (line 488) for selecting the parameter to vary.

## 4. Result Display Elements

Computation results are displayed using:

- **Standard Datagrok viewers.**
- **Custom HTMLElement.** Each custom display element is described in a separate UI component specification (similar to custom controls, see section 3).

By default, these elements are docked to the main view.

► **Implementation:**
- **Viewers:** `line chart` — `view.addViewer('Line chart', {...})` (line 552), showing Time vs Prey and Predators; `scatter plot` — `view.addViewer('Scatter plot', {...})` (line 559), showing the phase portrait (Prey vs Predators).
- **Custom element:** `statsPanel` — displays equilibrium point (`Prey* = γ/δ`, `Predator* = α/β`) and maximum population values, updated via `updateStatsPanel()` (line 65).
- **Color coding:** `applyColorCoding()` (line 221) sets linear color coding on `Prey` and `Predators` columns (dark red → yellow → dark green). Grid columns have `isTextColorCoded = true` (lines 42–43, 237–238).
- **Column header tooltips:** `columnTooltip()` (line 575) via `view.grid.onCellTooltip` (line 592) — shows column name, min (dark red), and max (dark green) for `Prey` and `Predators` columns.

## 5. Layout and UI Element Placement

The placement of controls and display elements (panels, ribbon, toolbar, side panels, docking area) is defined by the application specification.

► **Implementation (`./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`):**
- **Left panel:** `ui.form` with groups via `ui.h2` (lines 530–547), docked as `DG.DOCK_TYPE.LEFT`, ratio `0.6` (line 549).
- **Ribbon:** three groups — `[resetBtn]` (Actions), `[analyzeBtn]` (Analyze), `[helpBtn]` (Help) — set via `view.setRibbonPanels([[resetBtn], [analyzeBtn], [helpBtn]])` (line 527).
- **Main area:** `DG.TableView` (grid) — default.
- **Right area:** `line chart`, docked `DG.DOCK_TYPE.RIGHT`, ratio `0.7` (line 569); `scatter plot`, docked `DG.DOCK_TYPE.DOWN` relative to line chart, ratio `0.5` (line 572).

## 5.1. Styles

All visual styles of the application are placed in a separate CSS file (`css/<app-name>.css`). Inline styles in TypeScript code are not allowed — CSS classes are used instead.

- **Style isolation** — all custom CSS classes must use a unique application-specific prefix (`<app-name>-`) to prevent collisions with platform styles and other packages. For example, `lotka-volterra-app-stats-panel`, `lotka-volterra-app-tooltip-min`. Do not use generic class names without a prefix.
- **Static styles** — element styling that does not change during operation. Set via CSS class when creating the element.
- **Dynamic styles** — styles that depend on application state (e.g., indicator color, button activity). Implemented via CSS class toggling (`classList.toggle`, `classList.add/remove`), not via direct `element.style.*` assignment.

The CSS file is imported via ES import (`import '../css/<app-name>.css'`). Webpack with `style-loader` + `css-loader` injects styles into the DOM when the bundle loads.

► **Implementation:**
- CSS file: `./examples/lotka-volterra-spec/code/css/lotka-volterra.css` — contains classes `.lotka-volterra-app-stats-panel`, `.lotka-volterra-app-stats-value`, `.lotka-volterra-app-optimize-btn--disabled`, `.lotka-volterra-app-help-icon`, `.lotka-volterra-app-col-tooltip`, `.lotka-volterra-app-tooltip-range`, `.lotka-volterra-app-tooltip-min`, `.lotka-volterra-app-tooltip-max`.
- Import: `import '../../css/lotka-volterra.css'` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, line 13).
- **Static styles:** `statsPanel` is created with class `'lotka-volterra-app-stats-panel'` (line 63); stats value labels with `'lotka-volterra-app-stats-value'` (lines 51–54); `helpBtn` with `'lotka-volterra-app-help-icon'` (line 525).
- **Dynamic styles:** `optimizeBtn.classList.toggle('lotka-volterra-app-optimize-btn--disabled', !enabled)` (line 258).

## 6. User Feedback

### 6.1. Control Tooltips

- For standard Datagrok inputs, tooltips are defined at input creation time via the `tooltipText` property.
- For elements that are not Datagrok inputs, tooltips are bound via `ui.tooltip.bind`.

► **Implementation:**
- All inputs have `tooltipText` (e.g., `ctrlAlpha`: `tooltipText: 'Rate at which prey reproduce...'`, line 79).
- Buttons `optimizeBtn`, `resetBtn`, `analyzeBtn`: tooltip is passed as the third argument of `ui.bigButton` / `ui.iconFA` (lines 254, 505, 518).
- `helpBtn`: tooltip is passed as the second argument of `ui.icons.help` (line 522).

### 6.2. Validators as Feedback

Validators are added to standard Datagrok inputs — they display inline hints about the validity of the current value.

► **Implementation:** function `addValidators()` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, line 162) adds a validator to each input. The validator calls `validate(getInputs())` from the core and returns an error for the specific `InputId`.

### 6.3. Progress Bar

During long computations, the standard Datagrok progress bar indicator is displayed with cancellation support.

► **Implementation:** `task_optimize` uses `DG.TaskBarProgressIndicator.create('Optimizing Max Prey...')` (line 297) — indeterminate, cancelable. `task_equilibrium` uses `DG.TaskBarProgressIndicator.create('Exploring equilibrium...')` (line 423) — indeterminate, cancelable.

## 7. Validation

Validation is performed through the Datagrok validator mechanism. Each computational task has its own validation rules.

### 7.1. Complex Validation

Validation is complex: the entire set of task inputs is analyzed as a whole. If a combination of values is invalid, a `Map` with error specifications is returned. This specification is then used by the Datagrok validator mechanism to display errors on the corresponding inputs.

► **Implementation:**
- `task_primary`: function `validate(inputs: LotkaVolterraParams): ValidationErrors` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts`, line 56) — returns `Map<InputId, string>`. Rules val_01…val_07 check that each of the 7 parameters is positive. All rules are independent — no combinatorial dependencies.
- `task_optimize`: reuses `validate()` from the primary pipeline. If validation fails, optimization is rejected with `grok.shell.error('Invalid parameters. Fix inputs before optimizing.')`.
- `task_equilibrium`: no validation needed — the parameter choice is always valid.

### 7.2. Validation Order

The order of determining input set validity is defined by the specification for each task.

► **Implementation:** in `validate()`, all 7 rules are independent and checked in sequence (val_01 through val_07). There are no combinatorial rules that depend on other rules passing first.

## 8. Main Pipeline

Each computational task of the core has its own pipeline. All pipelines are orchestrated by the coordinator (see section 1.4).

### 8.1. Primary Pipeline

The primary pipeline is bound to the main application controls and operates reactively.

**Parameter input.** The user sets input values through controls. The UI adapter converts the input into typed data for the task's input port.

**Validation.** The input set is validated (see section 7). Validation is performed by the core through the input port.

**Computation.** Validated inputs are passed to the core task. Execution characteristics (synchronicity, environment, parallelization) are defined by the task specification (see section 1.1).

**Result display.** The core returns results through the task's output port. The coordinator passes them to the display adapter (see section 4).

► **Implementation:** function `runPrimary()` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, line 182):
1. Checks `computationsBlocked` (line 183).
2. Collects inputs: `getInputs()` (line 186).
3. Validates: `validate(inputs)` (line 187).
4. On errors — visually marks invalid inputs (`d4-invalid`) and calls `clearResults()` (lines 193–200).
5. On success — `solve(inputs)` (line 204), then `updateDataFrame(result)` + `updateStatsPanel(result)` (lines 205–206).

Reactive trigger: each input has `onValueChanged: () => debouncedRun()` (debounce 50 ms, line 214).

### 8.2. Secondary Pipelines

A secondary pipeline is bound to an explicit user action (button, icon, menu item) and runs on demand.

**Trigger.** The user initiates an action (e.g., clicks the "Optimize" button).

**Custom UI.** A secondary pipeline may have its own parameter input interface — for example, a Datagrok dialog with its own set of inputs, validation, and tooltips. This UI is defined by the secondary task specification.

**Validation.** Secondary task parameters are validated independently — using their own complex validation rules.

**Computation.** The core task is executed. The secondary task may use primary task results as part of its input data.

**Result display.** Secondary task results are displayed with their own elements — these can be additional viewers docked to the main view, dialog content, or a separate window. Defined by the specification.

**Feedback to primary pipeline.** The secondary task may return values that are substituted into primary pipeline controls. In this case, the computation blocking and batch update mechanism is used (see section 8.4).

► **Implementation of `task_optimize`:**
1. **Trigger:** `optimizeBtn` click (line 254).
2. **Custom UI:** none — uses current slider ranges and primary control values.
3. **Validation:** `validate(getInputs())` — reuses primary validation (line 263).
4. **Computation:** `runOptimization()` (line 261) — generates 14,641 grid points (11⁴), distributes across worker pool, collects results.
5. **Display:** no dedicated display — result fed back to primary pipeline.
6. **Feedback:** `ctrlAlpha.value = best.alpha`, `ctrlBeta.value = best.beta`, etc. via batch update (lines 376–387).

► **Implementation of `task_equilibrium`:**
1. **Trigger:** `analyzeBtn` → `showEquilibriumDialog()` (line 518).
2. **Custom UI:** dialog `ui.dialog('Explore Equilibrium Point')` with `ui.input.choice` for parameter selection (line 488).
3. **Validation:** none — parameter choice is always valid.
4. **Computation:** `runEquilibriumExploration(paramName)` (line 410) — single worker, 100,000 steps.
5. **Display:** new `DG.TableView` with DataFrame `[paramName, Prey*, Predator*]` and a line chart docked right (lines 466–485).
6. **Feedback:** none — results displayed in a separate view.

### 8.3. Common Pipeline Aspects

The following aspects are defined by the specification for each pipeline independently:

#### Control Behavior During Computations

Controls may become disabled during computations — which ones and for which pipeline is defined by the specification.

► **Implementation:** `task_primary` does not block controls (< 100 ms). `task_optimize` blocks only `optimizeBtn` via CSS class `lotka-volterra-app-optimize-btn--disabled` (line 258). `task_equilibrium` blocks all inputs via `input.enabled = false` and `optimizeBtn` via CSS class (function `setControlsEnabled()`, line 404).

#### Progress and Cancellation

During long computations, the standard Datagrok progress bar is displayed with cancellation support.

► **Implementation:** `task_primary` — no progress bar. `task_optimize` — `DG.TaskBarProgressIndicator` (indeterminate, cancelable) with cancellation via `pi.onCanceled` → `terminateWorkers()` (line 328). `task_equilibrium` — `DG.TaskBarProgressIndicator` (indeterminate, cancelable) with cancellation via `pi.onCanceled` → `worker.terminate()` (line 446).

#### Computation Error Handling

The course of action upon computation failure is defined by the specification for each task.

► **Implementation:**
- `task_primary`: `catch` → `clearResults()` + `grok.shell.error(msg)` (lines 207–211).
- `task_optimize`: partial errors (some workers failed) → `grok.shell.warning(...)` (line 361); complete failure → `grok.shell.error(...)` (line 364); error writing to controls → `grok.shell.warning(...)` (line 386).
- `task_equilibrium`: worker error → `grok.shell.error(...)` (line 462).

### 8.4. Computation Blocking and Batch Input Updates

In certain scenarios, the result of a secondary task must be substituted into primary pipeline controls. In this case, each individual input change should not trigger a reactive recomputation — the computation should happen once, after all values have been substituted.

The coordinator supports a computation blocking mode:

**Block request.** Before batch update, the coordinator suspends reactive execution of specified pipelines. Which pipelines are blocked is defined by the specification.

**Batch update.** Values are substituted into controls. Cascading dependencies between inputs (reactivity, section 9) can be handled in one of two modes — defined by the specification:

- Reactivity between inputs works as usual, but computations are not triggered.
- Reactivity between inputs is also suspended until the batch update completes.

**Unblock.** After all values are substituted, the coordinator removes the block. Full input set validation occurs, then computation runs, then results are displayed — the standard pipeline (section 8.1).

Example: the user launches parameter optimization (secondary task). Upon optimization completion, the coordinator blocks the primary pipeline, substitutes the found parameter values into all controls, removes the block — the primary pipeline runs once for the complete set of optimal parameters.

► **Implementation:** the `computationsBlocked` flag (`./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, line 22) is used in three scenarios:
1. **Format initialization** (lines 127–135): blocking prevents side-effect recomputations during `format` assignment.
2. **Reset** (lines 506–516): `computationsBlocked = true` → reset all 7 values → `computationsBlocked = false` → `runPrimary()`.
3. **Optimization result** (lines 376–387): `computationsBlocked = true` → set `ctrlAlpha.value`, `ctrlBeta.value`, `ctrlDelta.value`, `ctrlGamma.value` → `computationsBlocked = false` → `runPrimary()`.

Blocking check: `runPrimary()` first checks `if (computationsBlocked) return` (line 183).

## 9. Reactivity and Dependencies Between Inputs

Dependencies between inputs (cascading updates of ranges, defaults, availability) are defined by the application specification. Reactivity is managed by the coordinator (see section 1.4) and operates entirely at the UI adapter level — independent of the core's computational part.

Reactivity can be temporarily suspended by the coordinator in batch input update mode (see section 8.4).

► **Implementation:** in this application, there are no cascading dependencies between inputs — all 7 controls independently trigger `debouncedRun()` (debounce 50 ms, line 214). Any input change causes validation and recomputation via the primary pipeline. The `statsPanel` (equilibrium and max values) is updated as part of the primary pipeline result display, not as a separate reactive dependency.

## 10. Data Lifecycle

### 10.1. Data Input

The primary approach is manual entry through application controls.

► **Implementation:** all model parameters are entered by the user through the form in the left panel. The initial state uses default values from `DEFAULTS` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts`, line 32). `task_primary` runs automatically on initialization: `solve(DEFAULTS)` (line 30).

### 10.2. Loading from a Resource

Via buttons or icons — loading data from an external resource. Which specific resource and loading mechanism is defined by the application specification.

► **Implementation:** not used in this application.

## 11. Error Handling Beyond Computations

The strategy for handling data loading errors, network errors, invalid input files, and incorrect application state is defined by the application specification.

► **Implementation:** see specification (`./examples/lotka-volterra-spec/lotka-volterra-spec.md`, section 11). Examples:
- Worker creation error: `reject(new Error('Failed to start parallel computations...'))` (line 307).
- Partial worker errors: `errorCount` counting and `grok.shell.warning(...)` (lines 345, 361).
- Equilibrium worker error: `grok.shell.error(...)` (line 462).

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

All Datagrok event subscriptions (`onValueChanged`, `onAfterDraw`, etc.) must be collected and unsubscribed when the application closes via `sub.unsubscribe()`.

► **Implementation:** array `subs` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, line 24) collects subscriptions. On close — `for (const sub of subs) sub.unsubscribe()` (line 609).

### 12.2. Worker Termination

When the application closes, all web workers must be properly terminated.

► **Implementation:** array `activeWorkers` (line 25), function `terminateWorkers()` (line 390) calls `w.terminate()` for each worker. Called on close (line 607) and after optimization completes (line 338). Equilibrium worker is terminated individually on completion, error, or cancellation (lines 437, 442, 447).

## 13. Application Closure

When the application closes, the coordinator performs:

- All event subscriptions are unsubscribed — for both primary and secondary pipelines (see section 12.1).
- All web workers are terminated — including secondary task workers (see section 12.2).
- All associated UI elements are closed (including open secondary task dialogs).
- Pending requests are cancelled.

► **Implementation:** handler `grok.events.onViewRemoved.subscribe(...)` (`./examples/lotka-volterra-spec/code/src/lotka-volterra/app.ts`, lines 605–613):
1. `terminateWorkers()` — terminates all active optimization workers.
2. `for (const sub of subs) sub.unsubscribe()` — unsubscribes subscriptions.
3. `clearTimeout(debounceTimer)` — cancels pending debounce.

## 14. Accessibility and UX

Keyboard shortcuts, context menus, undo/redo, and other UX elements are defined by the application specification.

► **Implementation:** in the current version of the application, keyboard shortcuts, context menus, and undo/redo are not implemented (see specification, section 14). The `btn_reset` button resets all controls to default values as the only rollback mechanism. The `btn_help` button opens the Wikipedia article on Lotka-Volterra equations in a new browser tab.

## 15. Testing

### 15.1. Computational Part (Core)

Core correctness verification: unit tests for each computational task separately. The core is tested in isolation — without UI and adapters.

► **Implementation:** tests are split across two files:

**`./examples/lotka-volterra-spec/code/src/tests/lotka-volterra-api-tests.ts`** — 1 category:
- **API: Validation** — 18 tests: val_01…val_07 (zero and negative values for each of 7 parameters) + boundary values + valid defaults + multiple simultaneous errors.

**`./examples/lotka-volterra-spec/code/src/tests/lotka-volterra-math-tests.ts`** — 4 categories:
- **Math: LV func** — 5 tests: ODE right-hand side correctness at specific `(α, β, δ, γ, x, y)` points, including equilibrium.
- **Math: LV Equilibrium** — 2 tests: `getEquilibrium` verified against analytical formulas `x* = γ/δ`, `y* = α/β`.
- **Math: LV Solve properties** — 7 tests: output invariants of `solve()` — non-negativity, initial conditions, equilibrium stability, summary stats.
- **Math: MRT solver** — 2 tests: non-stiff 1D and stiff 1D reference problems verified against analytical solutions (Chapra & Canale).

Tests are run via `grok test` (entry point: `./examples/lotka-volterra-spec/code/src/package-test.ts`).

### 15.2. Inputs

Input verification for each task: all cases including edge cases (boundary values, invalid combinations, empty values, extreme values).

► **Implementation:** validation tests cover:
- Boundary values: `alpha = 0.1` (lower allowed bound), `x0 = 1` (lower allowed bound).
- Invalid values: zero and negative for each of the 7 parameters (14 tests).
- Multiple errors: simultaneously invalid `alpha=0`, `beta=0`, `x0=0`, `T=0` → ≥ 4 errors.
- Valid defaults: all default values pass validation.

### 15.3. Mathematical Verification

Verification that the implemented computation matches the model definition (see section 1.1, "Computation Formulas and Model"). Test categories correspond to the two levels of the model definition. Verification criteria — reference examples, expected accuracies, output property constraints, reference problems for the numerical method — are defined by the model specification, not invented during test writing. Tests implement what is specified; the specification is the source of truth.

#### Level 1 verification (required)

**Formula/equation verification.** The implemented transformation is checked at control points with manually computed expected values. For each computational path (mode, branch, regime), at least one test substitutes concrete inputs and compares the output against a hand-calculated result.

**Output property verification.** Constraints declared in the model definition (bounds, monotonicity, conservation laws, limiting cases) are checked on actual computation results. These tests do not compare against a specific expected value — they verify that the result satisfies a declared invariant.

► **Implementation:**
- **Formula verification:** `Math: LV func` — 5 tests. Each test substitutes specific `(α, β, δ, γ, x, y)` into the ODE right-hand side and compares `(dx/dt, dy/dt)` against hand-calculated values (e.g., `func_01`: `dx/dt = 1.0·10 − 0.1·10·5 = 5.0`, `dy/dt = 0.075·10·5 − 1.5·5 = −3.75`). Multiple parameter sets are covered (3 tests with default coefficients, 2 tests with different coefficients).
- **Equilibrium verification:** `Math: LV Equilibrium` — 2 tests. `getEquilibrium` is checked against analytical formulas `x* = γ/δ`, `y* = α/β` for two different parameter sets.
- **Output properties:** `Math: LV Solve properties` — 7 tests. Verifies invariants declared in the specification: `x(t) ≥ 0, y(t) ≥ 0` (solve_02), `x(0) = x₀, y(0) = y₀` (solve_03, solve_07), equilibrium stability (solve_05), `t[0] = 0` (solve_04), non-empty arrays (solve_01), summary stats correctness (solve_06).

#### Level 2 verification (for full formalization)

**Numerical method verification.** The solver (or library) is tested on reference problems with known analytical solutions. The test verifies that the numerical error stays within the expected tolerance. Reference problems should cover the solver's applicability range (e.g., non-stiff and stiff problems for an ODE solver). Each reference problem must cite its source (textbook, paper, test suite).

**Convergence verification.** Solving the same problem with decreasing step size or tolerance produces solutions that converge. The test compares solutions at two different precision levels and verifies that the discrepancy decreases.

**Asymptotic/equilibrium behavior.** The numerical solution on a sufficiently long interval approaches the analytically predicted equilibrium or asymptote.

► **Implementation:**
- **Numerical method:** `Math: MRT solver` — 2 tests. Non-stiff 1D and stiff 1D problems verified against analytical solutions (Chapra & Canale, pp. 736, 767). Tolerance threshold: max absolute error < 0.1.
- **Convergence:** not yet covered. Candidate: solve the Lotka-Volterra system with `tolerance = 1e-4` and `tolerance = 1e-8`, verify that the discrepancy between solutions decreases.
- **Asymptotic behavior:** partially covered in `solve_05`: starting at equilibrium `(x*, y*)`, verify that the solution stays near equilibrium after T=50.
