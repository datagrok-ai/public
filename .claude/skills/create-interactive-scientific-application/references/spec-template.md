# Datagrok Interactive Scientific Application Specification Template

> This template follows the structure of `guide.md`.
> Top-level section numbers (1–15) match the guide. Subsection numbering
> within sections 3, 4, 5, and 9 is template-specific (finer-grained than the guide).
> Fill in every section. Mark sections that do not apply as "N/A" with a brief explanation.

## 1. General Architecture

### 1.0. General Information

| Field | Value |
|---|---|
| Application name | |
| Package | |
| Entry function | |
| Brief description | What the application does, what scientific problem it solves |
| Main view type | `DG.TableView` / other |

### 1.1. Core

The core contains one or more computation tasks. Each task is described separately.

#### Task List

| Task ID | Name | Pipeline type | Trigger | Synchronicity | Execution environment | Parallelization |
|---|---|---|---|---|---|---|
| task_primary | ... | Primary (reactive) | Input change | Sync / async | Main thread / web worker | No / yes — strategy |
| task_secondary_1 | ... | Secondary (on demand) | Button / icon / menu item | ... | ... | ... |

#### Dependencies Between Tasks

```
Example:
  task_primary → result is used by task_secondary_1
  task_secondary_2 is independent of task_primary
```

#### Description of Each Task

A separate block is filled in for each task.

---

##### Task: `task_id`

**General characteristics:**

| Field | Value |
|---|---|
| Name | |
| Description | What the task computes |
| Dependency on other tasks | No / uses results of task `task_id` |

**Computation Formulas and Model**

> See guide section 1.1 "Computation Formulas and Model" for the two-level definition.

**Level 1 — required minimum (before implementation):**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| ... | ... | ... | e.g., `p ∈ (0, 1]` |

Relationships (equations, recurrences, algorithmic steps connecting inputs to outputs):

```
...
```

Output properties (invariants that must hold on the result):

| Property | Description |
|---|---|
| ... | e.g., bounds, monotonicity, conservation laws, limiting cases |

Reference examples (at least one per computational path / mode / branch):

| # | Inputs | Expected output | Computational path | Source |
|---|---|---|---|---|
| 1 | ... | ... | e.g., base model | Manual calculation / literature / reference implementation |

**Level 2 — full formalization (can be developed incrementally):**

- Complete mathematical formulation: ___(equations, initial/boundary conditions, parameterization)___
- Analytical properties: ___(equilibria, asymptotic behavior, stability, bifurcation points)___
- Numerical method justification: ___(why this method, stability, order of accuracy, applicability, literature reference)___

Level 2 document location: in this specification / separate document: ___

> Level 2 need not be complete before implementation begins, but must be complete before the computational part is considered verified.

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| ... | ... | ... | ... | ... |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| ... | ... | ... |

**Computation implementation:**

For each computation step:

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | Datagrok API | Method: ... (main thread only) | [Datagrok JS API](https://datagrok.ai/api/js/) |
| 2 | External library | Library: ..., version: ..., functions: ... | Link to API reference / README (required) |
| 3 | Custom method | Brief description: ... | Link to method specification (separate document, required) |

For **external libraries** implementing numerical methods: state which method properties are relevant (stability, order, applicability class), expected accuracy, and the verification strategy (reference problems, comparison with alternatives). See section 15.3.

For **custom methods**: the method specification (separate document) must contain: mathematical formulation, step-by-step algorithm, input/output data, constraints, edge cases, literature references, expected accuracy, and the verification strategy (reference examples with sources). See section 15.3.

**Execution environment constraint:** if the task uses Datagrok API — main thread only.

---

### 1.2. Ports

For each task, define input/output ports. Additionally, define application-level ports.

#### Task ports: `task_id`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | ... | Interface name / type | What parameters the task expects |
| Output | ... | Interface name / type | What the task returns |

#### Application-level ports

| Port | Used | Description |
|---|---|---|
| Progress | Yes / No | Interface for reporting execution progress (percentage, stage) |
| Cancellation | Yes / No | Interface for checking cancellation requests |
| Data | Yes / No | Interface for loading data from external resources |

### 1.3. Adapters

| Adapter | Used | Implementation |
|---|---|---|
| UI adapter | Yes / No | Datagrok inputs (`ui.input.*`), buttons, custom HTMLElement |
| Display adapter | Yes / No | Datagrok viewers, custom HTMLElement, docking |
| Worker adapter | Yes / No | Web worker wrapper for core tasks. Worker files in `src/<app-name>/workers/`. |
| Progress adapter | Yes / No | Datagrok progress bar |
| Data adapter | Yes / No | Loading from a specific resource |

For custom HTMLElements used as adapters, specify the component ID from section 3 (Custom UI Components).

### 1.4. Coordinator

The coordinator connects adapters and the core. Describe:

- How input changes are listened to (UI adapter).
- Reactivity management strategy (cascading dependencies, see section 9).
- How validation and computations are triggered.
- How control state is managed during computations.
- How computation blocking works (see section 8.4).
- How results are passed to the display adapter.
- Resource lifecycle management (subscriptions, workers — see section 12).

### 1.5. Independence Principle

Confirm that input behavior does not depend on the core's computational part. The UI adapter and reactivity form a standalone layer. The core receives a ready, validated set of parameters.

## 2. Main View

| Field | Value |
|---|---|
| View type | `DG.TableView` / custom |
| Description | ... |

## 3. Controls (Inputs)

### 3.1. Primary Pipeline Controls

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|
| ... | ... | `ui.input.int` / ... | `number` / ... | ... | ... | ... | e.g., `0.000` | Yes / No | ... | ... |

### 3.2. Secondary Task Triggers

Buttons, icons, menu items that launch secondary pipelines:

| ID | Label / icon | Launches task | Tooltip text | Availability condition |
|---|---|---|---|---|
| ... | ... | `task_id` | ... | ... |

### 3.3. Secondary Task Controls

For each secondary task that has its own UI:

#### Task controls: `task_id`

UI type: Datagrok dialog / other.

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text |
|---|---|---|---|---|---|---|---|---|---|
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

### 3.4. Other Buttons and Actions

Buttons that are not secondary task triggers (data loading, export, etc.):

| ID | Label / icon | Action | Tooltip text | Availability condition |
|---|---|---|---|---|
| ... | ... | ... | ... | ... |

### 3.5. Custom UI Components

If the application uses custom HTMLElements (controls or display elements):

- **Simple components** (static display panels built from `ui.divV`/`ui.label`
  with no complex interaction) — describe inline in this table.
- **Complex components** (custom interactive elements with multiple states,
  events, or non-trivial behavior) — describe in a separate UI component
  specification document and link from this table.

| Component ID | Brief description | Role (control / display / trigger) | UI component specification |
|---|---|---|---|
| ... | ... | ... | Inline / Link to separate document |

**UI component specification** (separate document, for complex components) must contain: visual description (sketch/mockup), states (normal, hover, disabled, active), events (what it emits on interaction), styles (CSS classes), accessibility (tooltips, aria).

## 4. Result Display Elements

### 4.1. Primary Pipeline Display Elements

| ID | Type | Associated output data (task.parameter) | Docking location |
|---|---|---|---|
| ... | Datagrok viewer (scatter plot / line chart / grid / ...) / custom HTMLElement | ... | ... |

For custom HTMLElements, specify the component ID from section 3.5.

### 4.2. Secondary Task Display Elements

For each secondary task:

#### Task display: `task_id`

Where results are displayed: additional viewers in main view / dialog content / separate window.

| ID | Type | Associated output data | Placement |
|---|---|---|---|
| ... | ... | ... | ... |

## 5. Layout and UI Element Placement

### 5.1. Control Placement

| Area | Content (control IDs) |
|---|---|
| Left panel | ... |
| Right panel | ... |
| Top panel / Ribbon | ... |
| Toolbar | ... |
| Main area | ... |

### 5.2. Display Element Placement

| Element ID | Docking area | Position / ratio |
|---|---|---|
| ... | ... | ... |

### 5.3. Styles

CSS file: `css/<app-name>.css`

All custom CSS classes must use a unique application-specific prefix (`<app-name>-`) to prevent collisions with platform styles and other packages (e.g., `lotka-volterra-app-stats-panel`, `lotka-volterra-app-tooltip-min`).

**Static styles** (do not change during operation):

| Element / class | CSS class(es) | Description |
|---|---|---|
| ... | ... | ... |

**Dynamic styles** (depend on application state — implemented via `classList.toggle/add/remove`):

| Element / class | CSS class(es) | Condition | Description |
|---|---|---|---|
| ... | ... | ... | ... |

CSS import: `import '../css/<app-name>.css'` in the main application file.

## 6. User Feedback

### 6.1. Control Tooltips

| Control ID | Tooltip text | Mechanism |
|---|---|---|
| ... | ... | `tooltipText` property / `ui.tooltip.bind` / `ui.iconFA` third argument |

### 6.2. Validators as Feedback

Validators are added to standard Datagrok inputs — they display inline hints about the validity of the current value.

| Input ID | Validation source | Description |
|---|---|---|
| ... | `validate()` from core / custom | ... |

### 6.3. Progress Bar

| Task | Progress bar | Type | Cancellation support |
|---|---|---|---|
| ... | Yes / No | Determinate / indeterminate | Yes / No |

## 7. Validation

### 7.1. Primary Pipeline Validation

#### Complex Validation Rules

| Rule ID | Condition (invalid combination) | Affected inputs (ID) | Error message |
|---|---|---|---|
| ... | ... | ... | ... |

#### Validation Order

```
1. Rule_A
2. Rule_B (checked only if Rule_A passed)
3. Rule_C
```

#### Returned Map Format

```
Map<inputId, errorMessage>
```

### 7.2. Secondary Task Validation

For each secondary task — a similar block.

#### Task validation: `task_id`

| Rule ID | Condition | Affected inputs (ID) | Error message |
|---|---|---|---|
| ... | ... | ... | ... |

Validation order:

```
1. ...
```

Return format: `{ errors: Map<inputId, errorMessage>, warning: string | null }`

## 8. Main Pipeline

### 8.1. Primary Pipeline

| Step | Description |
|---|---|
| Parameter input | User sets values through controls → UI adapter converts to typed data |
| Validation | Input set validated by core (section 7) |
| Computation | Validated inputs passed to core task |
| Result display | Results passed to display adapter (section 4) |

Reactive trigger: ___(e.g., `onValueChanged` with debounce ___ ms)___

Error behavior on validation failure: ___(clear results / keep previous / show message)___

### 8.2. Secondary Pipelines

For each secondary task:

#### Task pipeline: `task_id`

| Step | Description |
|---|---|
| Trigger | Button / icon / menu item: `control_id` |
| Custom UI | Dialog / panel with own inputs (section 3.3) |
| Validation | Independent validation rules (section 7.2) |
| Computation | Core task execution |
| Result display | Where and how results are shown (section 4.2) |
| Feedback to primary | Values substituted into primary controls: ___ / No |

### 8.3. Common Pipeline Aspects

Defined for each pipeline independently:

#### Control Behavior During Computations

| Pipeline / task | Controls blocked | Which controls |
|---|---|---|
| task_primary | Yes / No | ... |
| task_secondary_1 | Yes / No | ... |

#### Computation Error Handling

| Pipeline / task | Strategy | Notification method |
|---|---|---|
| task_primary | Reset results / keep previous / message | `grok.shell.error` / inline / ... |
| task_secondary_1 | ... | ... |

### 8.4. Computation Blocking and Batch Input Updates

| Scenario | Source (task) | Target controls (ID) | Blocked pipelines |
|---|---|---|---|
| ... | `task_id` | ... | Primary / ... |

Reactivity mode during batch update:

| Scenario | Reactivity mode |
|---|---|
| ... | Active but computations not triggered / Fully suspended |

## 9. Reactivity and Dependencies Between Inputs

### 9.1. Dependency Graph

| Source (input ID) | Target (input IDs) | Reaction type | Logic |
|---|---|---|---|
| ... | ... | Range / default / availability / option list / label update | ... |

### 9.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| ... | debounce / throttle / none | ... |

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via controls (section 3).

### 10.2. Loading from Resources

| Trigger (button ID) | Resource | Format | Mapping to inputs (ID) |
|---|---|---|---|
| ... | File / URL / DB / API | ... | ... |

### 10.3. Results Table Lifecycle

Update strategy: **DataFrame replacement** / **in-place mutation** / other.

Describe how the results DataFrame is created, updated, and replaced
throughout the application lifecycle:

```
1. Application initialization → ...
2. Primary task completion → ...
3. Secondary task completion → ...
4. Reset → ...
```

## 11. Error Handling Beyond Computations

| Error type | Strategy | Notification method |
|---|---|---|
| Data loading error | ... | `grok.shell.warning` / `grok.shell.error` / inline |
| Network error | ... | ... |
| Invalid input file | ... | ... |
| Worker creation error | ... | ... |
| Partial worker errors | ... | ... |

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

| Subscription | Event | Cleanup mechanism |
|---|---|---|
| ... | `onValueChanged` / `onAfterDraw` / ... | `sub.unsubscribe()` in cleanup handler |

All subscriptions must be collected and unsubscribed when the application closes.

### 12.2. Worker Termination

| Worker pool | Created in | Termination mechanism |
|---|---|---|
| ... | task / function name | `w.terminate()` in cleanup handler |

## 13. Application Closure

On view close, the coordinator performs:

- [ ] All event subscriptions unsubscribed (section 12.1)
- [ ] All web workers terminated (section 12.2)
- [ ] All associated UI elements closed (including open secondary task dialogs)
- [ ] Pending requests cancelled (debounce timers, in-flight operations)

Closure handler: ___(e.g., `grok.events.onViewRemoved.subscribe(...)`)___

## 14. Accessibility and UX

### 14.1. Keyboard Shortcuts

| Combination | Action |
|---|---|
| ... | ... |

### 14.2. Context Menus

| Context (element) | Menu items |
|---|---|
| ... | ... |

### 14.3. Undo / Redo

Supported: Yes / No.

If yes — which actions support rollback: ___

## 15. Testing

### 15.1. Computational Part (Core)

Core correctness verification: unit tests for each computational task separately. The core is tested in isolation — without UI and adapters.

Test files:

| File | Categories | Test count | Description |
|---|---|---|---|
| ... | ... | ... | ... |

Tests are run via `grok test` (entry point: `src/package-test.ts`).

### 15.2. Inputs

Input verification for each task: all cases including edge cases.

| Category | Coverage | Description |
|---|---|---|
| Boundary values | ... | e.g., lower/upper allowed bounds |
| Invalid combinations | ... | e.g., cross-parameter constraints |
| Dependencies between rules | ... | e.g., rule B skipped when rule A fails |
| Multiple simultaneous errors | ... | ... |

### 15.3. Mathematical Verification

> Verification criteria are defined by the model specification (section 1.1), not invented during test writing.
> Tests implement what is specified; the specification is the source of truth.

#### Level 1 Verification (required)

**Formula/equation verification:**

| Test category | Test count | What is verified | Reference source |
|---|---|---|---|
| ... | ... | Concrete inputs → expected output for each computational path | Manual calculation / literature |

**Output property verification:**

| Test category | Test count | Properties verified |
|---|---|---|
| ... | ... | e.g., bounds, initial conditions, convergence, monotonicity |

#### Level 2 Verification (for full formalization)

**Numerical method verification:**

| Test category | Test count | Reference problems | Tolerance | Source |
|---|---|---|---|---|
| ... | ... | e.g., non-stiff 1D, stiff 1D, stiff 2D | ... | Textbook / paper / test suite |

**Convergence verification:**

| Status | Description |
|---|---|
| Implemented / Not yet covered | e.g., solve with two tolerance levels, verify discrepancy decreases |

**Asymptotic/equilibrium behavior:**

| Status | Description |
|---|---|
| Implemented / Not yet covered | e.g., p(t_end) → p* for large t_end |
