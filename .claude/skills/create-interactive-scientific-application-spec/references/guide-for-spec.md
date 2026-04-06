# Architecture Guide for Specification Writing

> This is a condensed version of the full implementation guide. It contains
> only the architectural concepts, section expectations, and naming conventions
> needed to write a correct application specification. For implementation
> details, see the full guide used by the `implement-from-spec` skill.

---

## 1. General Architecture

An application is a Datagrok package function registered with `//tags: app`.
The architecture follows the **ports and adapters** pattern (hexagonal architecture)
with three layers and a coordinator.

### 1.1. Core

Computation logic. The core knows nothing about the UI and does not depend on
how data is obtained or how results are displayed.

The core contains one or more **computational tasks**. Each task is a
self-contained unit of computation. For each task, the specification must
independently define:

- Input parameters and their types.
- Output data and their format.
- Synchronicity: synchronous or asynchronous.
- Execution environment: main thread or web worker.
- Parallelization: whether parallelization is possible and the strategy.
- Complex validation logic for input parameters.

**Task taxonomy:**

- **Primary task** — triggered reactively when inputs change. Typically
  the main computation (e.g., solving an ODE, running a simulation).
- **Secondary task** — triggered by an explicit user action (button).
  Has its own parameters and may have its own UI (e.g., a dialog).

Tasks can be independent or linked (a secondary task may use results from
the primary one, or feed values back into primary controls).

**Core constraints:**

- Must NOT import `datagrok-api` or `ui`.
- Must NOT contain UI-prefixed identifiers. Use domain names (`alpha`,
  `beta`) — the coordinator maps them to control IDs (`ctrl_alpha`).
- Mathematical domains (e.g., `> 0`) belong to the core. Slider ranges
  (min/max for UI) belong to the UI adapter.
- Must be testable in isolation — each task tested separately.

**Execution environment rule:** tasks that use Datagrok API computation
methods can only run on the main thread. Tasks using only external libraries
and custom methods can run on either the main thread or in a web worker.

#### Computation Formulas and Model

The model definition has two levels — both should be addressed in the
specification (Section 1.1).

**Level 1 — required minimum (before implementation):**

- **Variables:** name, meaning, units, valid domain (e.g., `x ≥ 0`, `α > 0`).
- **Relationships:** equations, recurrences, algorithmic steps connecting
  inputs to outputs. Must be unambiguous — another developer must be able
  to independently implement the same computation.
- **Output properties:** invariants that must hold — bounds, monotonicity,
  symmetry, conservation laws, limiting/degenerate cases. These become
  verification tests.
- **Reference examples:** for each computational path, at least one concrete
  input → expected output pair with source (manual calculation, literature).

**Level 2 — full formalization (can be developed incrementally):**

- Complete mathematical formulation (equations, initial/boundary conditions).
- Analytical properties (equilibria, stability, bifurcation points).
- Numerical method justification (why this method, its properties, reference).

Level 2 need not be complete before implementation begins, but must be
complete before the computational part is considered verified.

#### Computation Implementation Approaches

For each task, specify which approach (or combination) is used:

- **Datagrok API methods** — available only on the main thread.
- **External libraries** — specify: which library, version, which
  functions/methods, link to documentation. For numerical libraries:
  which method properties matter, expected accuracy, verification strategy.
- **Custom methods** — each described in a separate method specification
  (formulas, algorithm, edge cases, references). The main spec includes
  a brief description and a link.

### 1.2. Ports

Interfaces through which the core communicates with the outside world.
Ports contain no implementation — only contracts. Each task has:

- **Input port** — what parameters and types the task expects.
- **Output port** — the format of the task's results.
- **Progress port** — reporting execution progress (percentage, stage).
- **Cancellation port** — checking whether the user requested cancellation.

At the application level:

- **Data port** — loading data from external resources.

The specification must define ports for every task. For simple applications
where tasks are synchronous and fast (< 100 ms), Progress and Cancellation
ports can be marked N/A.

### 1.3. Adapters

Concrete implementations of ports for the Datagrok environment. The
specification defines **what** each adapter does, not how it is coded.

- **UI adapter** — Datagrok inputs, buttons, custom elements. Converts
  user input into typed data for the task's input port.
- **Display adapter** — Datagrok viewers, custom elements, docking.
  Receives data from the output port and visualizes it.
- **Worker adapter** — runs the core in a web worker. Implements input
  and output ports via `postMessage`. Worker files go in
  `src/<app-name>/workers/`.
- **Progress adapter** — Datagrok progress bar with cancellation support.
- **Data adapter** — loading from a specific resource (defined by spec).

### 1.4. Coordinator

The coordinator connects adapters and the core. The specification must
describe the coordinator's responsibilities:

- Listening for input changes via the UI adapter.
- Managing reactivity: cascading dependencies between inputs (range updates,
  defaults, availability).
- Triggering validation and computations.
- Managing control state (enabled/disabled) during computations.
- Managing computation blocking for batch input updates (see Section 8.4).
- Passing results to the display adapter.
- Managing resource lifecycle (subscriptions, workers).

### 1.5. Independence Principle

Input behavior does not depend on the core's computational part. The UI
adapter and reactivity between inputs form a standalone layer managed by the
coordinator. The core receives a ready, validated set of parameters.

This means: when specifying controls (Section 3) and reactivity (Section 9),
describe them independently from the computation logic. The core should not
know about controls, their names, or their behavior.

---

## 2. Main View

The specification must define the main view type. If the application is based
on a table, use `DG.TableView`. Otherwise, describe the main view structure.

---

## 3. Controls

The specification must define every control with full detail.

**Control types:**

- **Standard Datagrok inputs** — created via `ui.input.*`.
- **Datagrok buttons** — action buttons.
- **Custom HTMLElement** — described in a separate UI component specification
  (visual description, states, events, styles, accessibility).

### 3.1. Input Options

For each standard Datagrok input, the specification must define:

- **Min** — minimum allowed value.
- **Max** — maximum allowed value.
- **Format** — display format (e.g., `0.000`, `0.0`, `0.##E+0`).

Input options differ from validation (Section 7): options set basic
control-level constraints, while validation checks complex conditions
across multiple input values.

### 3.2. Control Classification

Controls are classified by ownership:

- **Main view controls** — inputs of the primary pipeline, placed in
  the main application interface.
- **Secondary task triggers** — buttons/icons that launch secondary
  pipelines (see Section 8.2).
- **Secondary task controls** — inputs in the secondary task's own UI
  (e.g., in a dialog).

### Naming Convention

Use a consistent prefix:

- `ctrl_` for inputs: `ctrl_alpha`, `ctrl_x0`, `ctrl_drag_coeff`
- `btn_` for buttons: `btn_optimize`, `btn_reset`, `btn_help`
- `view_` for viewers: `view_line_chart`, `view_phase_plot`

---

## 4. Result Display Elements

The specification must define every display element:

- **Standard Datagrok viewers** — type, associated data columns, options.
- **Custom HTMLElement** — described in a separate UI component specification.

Display elements are docked to the main view by default. Specify docking
position and ratios.

---

## 5. Layout and UI Element Placement

The specification must define the placement of all controls and display
elements: panels, ribbon, toolbar, side panels, docking area. Include
docking types and ratios.

### 5.1. Styles

All styles go in a CSS file (`css/<app-name>.css`). The specification should
list the visual appearance requirements that need custom styling. Key rules:

- **Style isolation** — all CSS classes use `<app-name>-` prefix.
- **Static styles** — styling that does not change during operation.
- **Dynamic styles** — styles that depend on application state. Specify
  which states trigger which visual changes.

---

## 6. User Feedback

### 6.1. Control Tooltips

Every control and action button must have a tooltip. The specification must
define the tooltip text for each control — not a generic restatement of the
label, but a domain-aware explanation.

### 6.2. Validators as Feedback

Validators display inline hints about the validity of the current value
on standard Datagrok inputs. The specification defines which inputs have
validators and what they check (see Section 7).

### 6.3. Progress Bar

For long-running tasks, the specification must define:

- Whether a progress bar is shown.
- Whether it is determinate or indeterminate.
- Whether cancellation is supported.

---

## 7. Validation

Validation uses the Datagrok validator mechanism. Each task has its own
validation rules.

### 7.1. Complex Validation

The specification must define validation as analyzing the entire input set
as a whole. The result is a `Map<InputId, string>` with error messages per
input. Define:

- Every validation rule with a unique ID (e.g., `val_01`, `val_02`).
- The condition that triggers the error.
- The exact error message.
- Dependencies between rules (if any).

### 7.2. Validation Order

The specification must define the order in which rules are checked and whether
any rules depend on others passing first.

---

## 8. Main Pipeline

### 8.1. Primary Pipeline

The primary pipeline is reactive — bound to main controls. The specification
must describe the flow:

1. **Parameter input** — controls → UI adapter → typed data.
2. **Validation** — input set validated by the core.
3. **Computation** — validated inputs → core task → results.
4. **Result display** — results → display adapter.

Also specify: debounce strategy (if any) and what happens on validation
failure (clear results, show errors, etc.).

### 8.2. Secondary Pipelines

Each secondary pipeline is bound to an explicit user action. The
specification must define for each:

1. **Trigger** — what initiates the pipeline (button, icon, menu item).
2. **Custom UI** — any dedicated parameter interface (dialog, panel).
3. **Validation** — independent validation rules for the secondary task.
4. **Computation** — the core task, its execution characteristics.
5. **Result display** — how and where results are shown (docked viewer,
   new TableView, dialog content).
6. **Feedback to primary pipeline** — whether results are substituted
   into primary controls (triggers computation blocking, see 8.4).

### 8.3. Common Pipeline Aspects

For each pipeline, the specification must define:

- **Control behavior during computations** — which controls are disabled
  and for which pipeline.
- **Progress and cancellation** — see Section 6.3.
- **Error handling** — what happens on computation failure (error message,
  clear results, partial results).

### 8.4. Computation Blocking and Batch Input Updates

When a secondary task feeds values back into primary controls, each
individual change should not trigger a recomputation. The specification
must define:

- **Which pipelines are blocked** during batch update.
- **Reactivity mode during batch update:** (a) reactivity between inputs
  works but computations don't trigger, or (b) reactivity is also suspended.
- **Unblock behavior:** validation → computation → display after all
  values are substituted.

For simple applications without secondary-to-primary feedback, this section
is N/A.

---

## 9. Reactivity and Dependencies Between Inputs

The specification must define all cascading dependencies between inputs:

- When input A changes, does it update the range/default/availability
  of input B?
- Are there groups of inputs that change together?
- Is reactivity suspended during batch updates (see Section 8.4)?

If there are no dependencies (all inputs independently trigger computation),
state this explicitly.

---

## 10. Data Lifecycle

### 10.1. Data Input

The specification must define the primary data input method:

- Manual entry through controls (most common for simulations).
- File upload.
- External data source.

### 10.2. Loading from a Resource

If the application loads data from external resources, the specification
must define: which resource, loading mechanism, error handling, and
how loaded data integrates with the application state.

---

## 11. Error Handling Beyond Computations

The specification must define error handling strategy for:

- Data loading errors.
- Network errors.
- Invalid input files.
- Incorrect application state.
- Worker creation failures.
- Partial failures (some workers failed, others succeeded).

---

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

The specification must list all event subscriptions the application will
use (e.g., `onValueChanged`, `onAfterDraw`, `onViewRemoved`). All
subscriptions must be collected and unsubscribed on close.

### 12.2. Worker Termination

If the application uses web workers, the specification must define:

- When workers are created and terminated.
- How workers are cleaned up on application close.
- How workers are terminated on cancellation.

---

## 13. Application Closure

The specification must define what happens when the application closes:

- All event subscriptions are unsubscribed.
- All web workers are terminated.
- All associated UI elements are closed.
- Pending requests are cancelled.

---

## 14. Accessibility and UX

The specification must define (or explicitly mark as not implemented):

- Keyboard shortcuts.
- Context menus.
- Undo/redo mechanism.
- Any other UX elements (e.g., help button, reset button).

---

## 15. Testing

The specification must define test categories and expectations.

### 15.1. Computational Part (Core)

Unit tests for each computational task separately, tested in isolation
without UI and adapters. Define test categories and expected coverage.

### 15.2. Inputs

Input verification for each task: boundary values, invalid values,
invalid combinations, empty values, extreme values. Define specific
test cases.

### 15.3. Mathematical Verification

Test categories corresponding to the two levels of the model definition:

**Level 1 (required):**

- **Formula/equation verification** — control points with manually
  computed expected values. At least one test per computational path.
- **Output property verification** — invariants declared in the model
  definition (bounds, monotonicity, conservation laws) checked on
  actual results.

**Level 2 (for full formalization):**

- **Numerical method verification** — reference problems with known
  analytical solutions. Specify: which problems, sources, tolerance.
- **Convergence verification** — decreasing step size/tolerance produces
  converging solutions.
- **Asymptotic/equilibrium behavior** — numerical solution approaches
  predicted equilibrium on sufficiently long intervals.

For each test category, the specification defines: what is tested,
expected results, and their sources. Tests implement what is specified —
the specification is the source of truth.
