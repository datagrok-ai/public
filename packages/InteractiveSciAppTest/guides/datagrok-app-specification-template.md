# Datagrok Interactive Scientific Application Specification Template

## 1. General Information

| Field | Value |
|---|---|
| Application name | |
| Package | |
| Entry function | |
| Brief description | What the application does, what scientific problem it solves |
| Main view | `DG.TableView` / other |

## 2. Computation Tasks (Core)

The core contains one or more computation tasks. Each task is described separately.

### 2.1. Task List

| Task ID | Name | Pipeline type | Trigger |
|---|---|---|---|
| task_primary | ... | Primary (reactive) | Input change |
| task_secondary_1 | ... | Secondary (on demand) | Button / icon / menu item |
| ... | ... | ... | ... |

### 2.2. Description of Each Task

A separate block is filled in for each task.

---

#### Task: `task_id`

**General characteristics:**

| Field | Value |
|---|---|
| Name | |
| Description | What the task computes |
| Synchronicity | Synchronous / asynchronous |
| Execution environment | Main thread / web worker |
| Parallelization | No / yes — strategy |
| Dependency on other tasks | No / uses results of task `task_id` |

**Task input parameters:**

| Parameter | Type | Description |
|---|---|---|
| ... | ... | ... |

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

**Method specification** (separate document) must contain: mathematical formulation (formulas, equations), step-by-step algorithm, input and output data, constraints and assumptions, edge cases, literature references.

If the usage of an external library is non-trivial: link to a separate agreement document or detailed description.

**Execution environment constraint:** if the task uses Datagrok API — main thread only.

---

### 2.3. Dependencies Between Tasks

```
Example:
  task_primary → result is used by task_secondary_1
  task_secondary_2 is independent of task_primary
```

## 3. Controls

### 3.1. Primary Pipeline Controls

Main view inputs:

| ID | Label | Control type | Data type | Default | Range / allowed values | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|
| ... | ... | `ui.input.int` / ... | `number` / ... | ... | ... | Yes / No | ... | ... |

For custom HTMLElements, specify the component ID from section 3.5 as the control type.

### 3.2. Secondary Task Triggers

Buttons, icons, menu items that launch secondary pipelines:

| ID | Label / icon | Launches task | Tooltip text | Availability condition |
|---|---|---|---|---|
| ... | ... | `task_id` | ... | ... |

### 3.3. Secondary Task Controls

For each secondary task that has its own UI:

#### Task controls: `task_id`

UI type: Datagrok dialog / other.

| ID | Label | Control type | Data type | Default | Range / allowed values | Nullable | Tooltip text |
|---|---|---|---|---|---|---|---|
| ... | ... | ... | ... | ... | ... | ... | ... |

For custom HTMLElements, specify the component ID from section 3.5 as the control type.

### 3.4. Other Buttons and Actions

Buttons that are not secondary task triggers (data loading, export, etc.):

| ID | Label / icon | Action | Tooltip text | Availability condition |
|---|---|---|---|---|
| ... | ... | ... | ... | ... |

### 3.5. Custom UI Components

If the application uses custom HTMLElements (controls or display elements), each is described in a separate document — a UI component specification.

| Component ID | Brief description | Role (control / display / trigger) | UI component specification |
|---|---|---|---|
| ... | ... | ... | Link to separate document (required) |

**UI component specification** (separate document) must contain: visual description (sketch/mockup), states (normal, hover, disabled, active), events (what it emits on interaction), styles (CSS classes), accessibility (tooltips, aria).

## 4. Validation

Each task has its own validation rules.

### 4.1. Primary Pipeline Validation

#### Complex Validation Rules

| Rule ID | Condition (invalid combination) | Affected inputs (ID) | Error message |
|---|---|---|---|
| ... | ... | ... | ... |

#### Validation Order

```
1. Rule_A
2. Rule_B
3. Rule_C
```

#### Returned Map Format

```
Map<inputId, errorMessage>
```

### 4.2. Secondary Task Validation

For each secondary task — a similar block.

#### Task validation: `task_id`

| Rule ID | Condition | Affected inputs (ID) | Error message |
|---|---|---|---|
| ... | ... | ... | ... |

Validation order:

```
1. ...
```

## 5. Reactivity and Dependencies Between Inputs

### 5.1. Dependency Graph

| Source (input ID) | Target (input IDs) | Reaction type | Logic |
|---|---|---|---|
| ... | ... | Range / default / availability / option list update | ... |

### 5.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| ... | debounce / throttle / none | ... |

## 6. Behavior During Computations

Defined independently for each pipeline.

### 6.1. Primary Pipeline

#### Control Blocking

| Control ID | Blocked | Note |
|---|---|---|
| ... | Yes / No | ... |

#### Progress Bar

| Field | Value |
|---|---|
| Display | Yes / No |
| Type | Determinate (with percentage) / indeterminate |
| Cancellation support | Yes / No |

#### Error Behavior

| Strategy | Description |
|---|---|
| Reset results | Clear display area |
| Last valid | Keep previous results |
| Message | Show error text to user |

Selected strategy: ___

### 6.2. Secondary Pipelines

For each secondary task — a similar block.

#### Task pipeline: `task_id`

Control blocking:

| Control ID | Blocked | Note |
|---|---|---|
| ... | Yes / No | ... |

Progress bar:

| Field | Value |
|---|---|
| Display | Yes / No |
| Type | Determinate / indeterminate |
| Cancellation support | Yes / No |

Error behavior: ___

## 7. Computation Blocking and Batch Update

Fill in if a secondary task returns values to primary pipeline controls.

### 7.1. Batch Update Scenarios

| Source (task) | Target controls (ID) | Blocked pipelines |
|---|---|---|
| `task_id` | ... | Primary / ... |

### 7.2. Reactivity Mode During Batch Update

| Scenario | Reactivity mode |
|---|---|
| ... | Active but computations are not triggered / Fully suspended |

## 8. Results Display

### 8.1. Primary Pipeline Display Elements

| ID | Type | Associated output data (task.parameter) | Docking location |
|---|---|---|---|
| ... | Datagrok viewer (scatter plot / line chart / grid / ...) / custom HTMLElement | ... | ... |

For custom HTMLElements, specify the component ID from section 3.5 (UI component specification).

### 8.2. Secondary Task Display Elements

For each secondary task:

#### Task display: `task_id`

Where results are displayed: additional viewers in main view / dialog content / separate window.

| ID | Type | Associated output data | Placement |
|---|---|---|---|
| ... | ... | ... | ... |

For custom HTMLElements, specify the component ID from section 3.5.

## 9. Layout

### 9.1. Control Placement

| Area | Content (control IDs) |
|---|---|
| Left panel | ... |
| Right panel | ... |
| Top panel / Ribbon | ... |
| Toolbar | ... |
| Main area | ... |

### 9.2. Display Element Placement

| Element ID | Docking area | Position |
|---|---|---|
| ... | ... | ... |

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via controls (section 3).

### 10.2. Loading from Resources

| Trigger (button ID) | Resource | Format | Mapping to inputs (ID) |
|---|---|---|---|
| ... | File / URL / DB / API | ... | ... |

## 11. Error Handling Outside of Computations

| Error type | Strategy | Notification method |
|---|---|---|
| Data loading error | ... | `grok.shell.warning` / `grok.shell.error` / inline |
| Network error | ... | ... |
| Invalid input file | ... | ... |
| Incorrect state | ... | ... |

## 12. UX

### 12.1. Keyboard Shortcuts

| Combination | Action |
|---|---|
| ... | ... |

### 12.2. Context Menus

| Context (element) | Menu items |
|---|---|
| ... | ... |

### 12.3. Undo / Redo

Supported: Yes / No.

If yes — which actions support rollback: ___
