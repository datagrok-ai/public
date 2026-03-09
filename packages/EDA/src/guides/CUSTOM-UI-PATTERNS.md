# Custom UI Patterns (Ports and Adapters)

This guide describes the UI architecture for analytical methods in the EDA package.
For computation patterns see `COMPUTATION-PATTERNS.md`.
For worker patterns see `WORKER-GUIDE.md`.

## Architecture

The UI is built using the **Ports and Adapters** pattern (hexagonal architecture).
Three layers isolated through typed contracts:

```
┌───────────────┐      ┌─────────────┐      ┌────────────────┐
│ Input Adapter │─────>│    Core      │─────>│ Output Adapter │
│ (UI controls) │  P₁  │ (computation)│  P₂  │ (viewers)      │
└───────────────┘      └─────────────┘      └────────────────┘
```

- **Port P₁ (input)** — the `ComputeInputs` interface, describes *what the core needs*.
- **Port P₂ (output)** — the `ComputeOutputs` interface, describes *what the core returns*.
- **Core** — a pure function `compute(inputs): outputs`. Has no knowledge of the UI.
- **Input Adapter** — collects data from `ui.input.*` controls into `ComputeInputs`.
- **Output Adapter** — takes `ComputeOutputs` and distributes them across viewers.

### Why

- The core can be tested without the UI.
- The input adapter is replaceable (dialog, panel form, programmatic call).
- The output adapter is replaceable (dock panels, dialog, file export).
- Each layer changes independently.

---

## Step 1. Ports (interfaces)

Ports are data types that describe the contract between layers.
They depend on neither the UI nor any specific computation implementation.

### Input port

Describes all parameters required by the core:

```typescript
/** Computation input data */
type ComputeInputs = {
  table: DG.DataFrame;
  descriptors: DG.Column[];
  desirability: DG.Column;
  pValue: number;
  r2: number;
  qCutoff: number;
  useSigmoid: boolean;
};
```

### Output port

Describes all computation results:

```typescript
/** Computation results */
type ComputeOutputs = {
  prediction: DG.Column;
  statsTable: DG.DataFrame;
  roc: { fpr: Float32Array; tpr: Float32Array; auc: number; threshold: number };
  selectedDescriptors: string[];
};
```

### Validation port

Describes the result of input data validation:

```typescript
/** Validation result */
type ValidationResult = {
  valid: boolean;
  errors: Map<string, string>;  // inputId → error message
};
```

### Rules for defining ports

- Ports are `type` or `interface`, not classes.
- Ports use Datagrok types (`DG.DataFrame`, `DG.Column`) since these are core platform data types, not UI.
- Ports must not contain references to `ui.*`, `DG.Viewer`, or `HTMLElement`.
- One method — one set of ports. Do not attempt to create a universal port for all methods.

---

## Step 2. Core (computation)

The core is a pure function (or a set of functions) that takes `ComputeInputs` and returns `ComputeOutputs`.

### Rules

- The core **does not import** `ui` and does not create DOM elements.
- The core **does not show** warnings (`grok.shell.warning`). Instead, it throws typed errors.
- The core **has no knowledge** of viewers, forms, or dialogs.
- The core may import `grok` and `DG` (for data manipulation).

### Example

```typescript
import * as DG from 'datagrok-api/dg';

/** Computation error — thrown by the core, caught by the wiring layer */
class ComputeError extends Error {
  constructor(message: string) { super(message); }
}

/** Input data validation (no side effects) */
function validate(inputs: ComputeInputs): ValidationResult {
  const errors = new Map<string, string>();

  if (inputs.descriptors.length < 1)
    errors.set('descriptors', 'Select at least one column.');

  if (inputs.pValue < 0.001 || inputs.pValue > 1)
    errors.set('pValue', 'p-value must be between 0.001 and 1.');

  const zeroVarCols = inputs.descriptors
    .filter((col) => col.stats.stdev === 0)
    .map((col) => col.name);

  if (zeroVarCols.length > 0)
    errors.set('descriptors', `Zero variance: ${zeroVarCols.join(', ')}`);

  return { valid: errors.size === 0, errors };
}

/** Main computation */
function compute(inputs: ComputeInputs): ComputeOutputs {
  const validation = validate(inputs);
  if (!validation.valid)
    throw new ComputeError([...validation.errors.values()].join('; '));

  // ... computations ...

  return { prediction, statsTable, roc, selectedDescriptors };
}
```

---

## Step 3. Input Adapter (UI controls)

The input adapter creates controls and collects `ComputeInputs` from them.
It knows about `ui.input.*`, but **has no knowledge** of viewers and does not invoke computations directly.

### Input adapter contract

```typescript
/** The input adapter returns a form and a data collection function */
type InputAdapter = {
  form: HTMLElement;                        // root DOM element of the form
  collectInputs: () => ComputeInputs;       // collect current values
  onChanged: (cb: () => void) => void;      // subscribe to any change
  validate: () => ValidationResult;         // validate all controls
};
```

### Datagrok input catalog

All inputs are created via `ui.input.*`. Each is an instance of `DG.InputBase`:
- `.value` — current value
- `.addValidator(fn)` — register a validator (returns an error string or `null`)
- `.validate()` — run validation (returns `boolean`)
- `.enabled` / `.root.hidden` — availability / visibility
- `.setTooltip(text)` — tooltip

#### Column selection

```typescript
// Single column
ui.input.column('Desirability', {
  table: df,
  nullable: false,
  value: initialCol,
  filter: (col) => col.isNumerical,
  tooltipText: 'Column with target values.',
  onValueChanged: (value) => { /* ... */ },
});

// Multiple columns
ui.input.columns('Descriptors', {
  table: df,
  nullable: false,
  available: numericColNames,
  checked: preselectedNames,
  tooltipText: 'Descriptor columns.',
  onValueChanged: (value) => { /* ... */ },
});
```

#### Numeric parameters

```typescript
ui.input.float('p-value', {
  nullable: false,
  min: 0.001, max: 1.0, step: 0.001,
  value: 0.05,
  format: '0.000',
  tooltipText: 'p-value threshold.',
  onValueChanged: (value) => { /* ... */ },
});
```

#### Dropdown selection

```typescript
ui.input.choice('Condition', {
  value: '>=',
  items: ['>', '<', '>=', '<='],
  nullable: false,
  tooltipText: 'Comparison sign.',
  onValueChanged: (value) => { /* ... */ },
});
```

#### Boolean parameters

```typescript
ui.input.bool('Auto-tuning', {
  value: false,
  tooltipText: 'Automatic parameter tuning.',
  onValueChanged: (value) => { /* ... */ },
});
```

#### Multi-choice categories

```typescript
ui.input.multiChoice('Preferred', {
  value: preselected,
  items: allCategories,
  nullable: false,
  tooltipText: 'Select desired categories.',
  onValueChanged: (value) => { /* ... */ },
});
```

#### Buttons

```typescript
ui.button('Save', async () => { /* ... */ }, 'Save model.');
ui.bigButton('Run', async () => { /* ... */ });
```

#### Custom elements

```typescript
const el = ui.div('Click me', 'my-css-class');
el.addEventListener('click', () => { /* ... */ });
```

`ui.div()`, `ui.divV()`, `ui.divH()` — regular `HTMLDivElement`, any DOM events.

#### Toggle icons

```typescript
ui.icons.settings(
  () => { metricsDiv.hidden = !metricsDiv.hidden; },
  'Show/hide metrics settings',
);
```

### Input adapter example

```typescript
function createInputAdapter(df: DG.DataFrame): InputAdapter {
  // --- Create controls ---

  const descrInput = ui.input.columns('Descriptors', {
    table: df,
    available: getNumericColNames(df),
    onValueChanged: () => notifyChanged(),
  });

  const desInput = ui.input.column('Desirability', {
    table: df,
    filter: (col) => col.isNumerical || col.type === DG.COLUMN_TYPE.BOOL,
    onValueChanged: (value) => {
      if (value != null) {
        updateAuxInputs(value);
        notifyChanged();
      }
    },
  });

  const pInput = ui.input.float('p-value', {
    nullable: false, min: 0.001, max: 1.0, step: 0.001, value: 0.05,
    onValueChanged: () => notifyChanged(),
  });

  const rInput = ui.input.float('R²', {
    nullable: false, min: 0.0, max: 1.0, step: 0.01, value: 0.5,
    onValueChanged: () => notifyChanged(),
  });

  const qInput = ui.input.float('q-cutoff', {
    nullable: false, min: 0.0, max: 1.0, step: 0.01, value: 0.1,
    onValueChanged: () => notifyChanged(),
  });

  const sigmaInput = ui.input.bool('σ correction', {
    value: true,
    onValueChanged: () => notifyChanged(),
  });

  // --- Register validators ---

  descrInput.addValidator(() => {
    if (descrInput.value == null || descrInput.value.length < 1)
      return 'Select at least one column.';
    return null;
  });

  pInput.addValidator(() => {
    if (pInput.value == null || pInput.value < 0.001 || pInput.value > 1)
      return 'p-value: from 0.001 to 1.';
    return null;
  });

  // --- Form layout ---

  const form = ui.form([]);
  form.append(ui.h2('Training data'));
  form.append(descrInput.root);
  form.append(desInput.root);
  form.append(ui.h2('Settings'));
  form.append(sigmaInput.root);
  form.append(pInput.root);
  form.append(rInput.root);
  form.append(qInput.root);

  // --- Conditional visibility ---

  const auxInputsDiv = ui.divV([/* dependent controls */]);
  const updateAuxInputs = (col: DG.Column) => {
    auxInputsDiv.hidden = col.type === DG.COLUMN_TYPE.BOOL;
  };

  // --- Change subscription ---

  let changeCallback: (() => void) | null = null;
  const notifyChanged = () => { if (changeCallback) changeCallback(); };

  // --- Contract ---

  return {
    form,
    collectInputs: () => ({
      table: df,
      descriptors: descrInput.value ?? [],
      desirability: desInput.value!,
      pValue: pInput.value!,
      r2: rInput.value!,
      qCutoff: qInput.value!,
      useSigmoid: sigmaInput.value,
    }),
    onChanged: (cb) => { changeCallback = cb; },
    validate: () => {
      const errors = new Map<string, string>();
      if (!descrInput.validate()) errors.set('descriptors', 'Invalid');
      if (!pInput.validate()) errors.set('pValue', 'Invalid');
      if (!rInput.validate()) errors.set('r2', 'Invalid');
      if (!qInput.validate()) errors.set('qCutoff', 'Invalid');
      return { valid: errors.size === 0, errors };
    },
  };
}
```

### Advanced input adapter patterns

#### Grouping auxiliary inputs

```typescript
const auxDiv = ui.divV([signInput.root, thresholdInput.root]);
form.append(auxDiv);
auxDiv.hidden = someCondition;  // show/hide the entire group
```

#### Dynamic input addition/removal

```typescript
let dynamicInput: DG.InputBase<string[] | null> | null = null;

const updateDynamicInput = () => {
  if (dynamicInput != null) {
    dynamicInput.root.remove();
    dynamicInput = null;
  }

  if (needsDynamicInput()) {
    dynamicInput = ui.input.multiChoice('Categories', {
      items: getCategories(),
      onValueChanged: () => notifyChanged(),
    });
    dynamicInput.addValidator(() => /* ... */);
    containerDiv.append(dynamicInput.root);
  }
};
```

#### Blocking cascading updates

When one input programmatically changes another:

```typescript
let blocked = false;

colInput.onValueChanged = (value) => {
  blocked = true;
  thresholdInput.value = Math.round(value.stats.avg * 100) / 100;
  blocked = false;
  notifyChanged();
};

thresholdInput.onValueChanged = () => {
  if (blocked) return;
  notifyChanged();
};
```

#### Conditional availability

```typescript
const setEnability = (toEnable: boolean) => {
  pInput.enabled = toEnable;
  rInput.enabled = toEnable;
  qInput.enabled = toEnable;
};

autoTuneInput.onValueChanged = (value) => {
  setEnability(!value);
};
```

---

## Step 4. Output Adapter (viewers)

The output adapter takes `ComputeOutputs` and updates visual components.
It knows about `DG.Viewer`, but **has no knowledge** of input controls and does not invoke computations.

### Output adapter contract

```typescript
/** Output adapter — viewers and an update function */
type OutputAdapter = {
  viewers: DG.Viewer[];                       // all viewers for placement
  update: (outputs: ComputeOutputs) => void;  // update all viewers
};
```

### Viewer catalog

#### Grid viewer

```typescript
const statGrid = DG.Viewer.grid(statsDataFrame, {
  showTitle: true,
  title: 'Statistics',
});

statGrid.setOptions({ rowHeight: 75 });

const col = statGrid.col('Desirability');
col!.width = 200;
col!.cellType = 'html';
```

#### Line chart

```typescript
const rocViewer = DG.Viewer.lineChart(rocDataFrame, {
  xColumnName: 'FPR',
  yColumnName: 'TPR',
  linesWidth: 5,
  markerType: 'dot',
  title: `ROC Curve (AUC = ${auc.toFixed(3)})`,
});

// Formula line
rocDataFrame.meta.formulaLines.addLine({
  title: 'Baseline',
  formula: '${TPR} = ${FPR}',
  width: 1, style: 'dashed', min: 0, max: 1,
});
```

#### Specialized viewers

```typescript
const confMatrix = DG.Viewer.fromType('Confusion matrix', df, {
  xColumnName: desColName,
  yColumnName: predColName,
  description: `Threshold: ${threshold.toFixed(3)}`,
});
```

### Updating viewers

```typescript
// Assign a new dataFrame + setOptions
rocViewer.dataFrame = newRocDf;
rocViewer.setOptions({ title: `ROC (AUC = ${newAuc.toFixed(3)})` });
```

### Grid cell tooltips

```typescript
grid.onCellTooltip((cell, x, y) => {
  if (cell.isColHeader) {
    const tip = tooltipMap.get(cell.tableColumn!.name);
    if (tip) {
      ui.tooltip.show(ui.divV([ui.p(tip)]), x, y);
      return true;
    }
  }
  return false;
});
```

### HTML in grid cells

```typescript
grid.onCellPrepare((cell) => {
  if (!cell.isTableCell || cell.tableColumn?.name !== 'Profile')
    return;

  const key = grid.cell('Name', cell.gridRow).value;
  const element = htmlElementsMap.get(key);
  if (element != null)
    cell.element = element;
});
```

### Color coding

```typescript
// Categorical
column.colors.setCategorical({
  'Selected': DG.Color.fromHtml('#1a921a'),
  'Excluded': DG.Color.fromHtml('#d03943'),
});

// Linear gradient
column.colors.setLinear(palette, { min: column.stats.min, max: column.stats.max });
```

### Updating the main grid

```typescript
df.columns.remove(predictionName);
df.columns.add(newPrediction);

view.grid.sort([predictionName], [false]);

const gridCol = view.grid.col(predictionName);
gridCol!.format = '0.0000';
gridCol!.isTextColorCoded = true;
```

### Output adapter example

```typescript
function createOutputAdapter(df: DG.DataFrame): OutputAdapter {
  // --- Create viewers ---

  const statsGrid = DG.Viewer.grid(DG.DataFrame.create(), { showTitle: true, title: 'Statistics' });

  const rocViewer = DG.Viewer.lineChart(DG.DataFrame.create(), {
    xColumnName: 'FPR', yColumnName: 'TPR',
    linesWidth: 5, markerType: 'dot',
  });

  const confMatrix = DG.Viewer.fromType('Confusion matrix', df);

  // --- Update ---

  const update = (outputs: ComputeOutputs) => {
    // Statistics grid
    statsGrid.dataFrame = outputs.statsTable;

    // ROC curve
    const rocDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('FPR', outputs.roc.fpr),
      DG.Column.fromFloat32Array('TPR', outputs.roc.tpr),
    ]);
    rocDf.meta.formulaLines.addLine({
      title: 'Baseline', formula: '${TPR} = ${FPR}',
      width: 1, style: 'dashed', min: 0, max: 1,
    });
    rocViewer.dataFrame = rocDf;
    rocViewer.setOptions({ title: `ROC (AUC = ${outputs.roc.auc.toFixed(3)})` });

    // Prediction column
    df.columns.remove(outputs.prediction.name);
    df.columns.add(outputs.prediction);

    // Confusion matrix
    confMatrix.setOptions({
      xColumnName: 'Desirability',
      yColumnName: outputs.prediction.name,
      description: `Threshold: ${outputs.roc.threshold.toFixed(3)}`,
    });
  };

  return {
    viewers: [statsGrid, rocViewer, confMatrix],
    update,
  };
}
```

---

## Step 5. Wiring

The wiring is a thin layer that connects the input adapter, the core, and the output adapter.
This is the only place where all three components meet.

### Basic wiring

```typescript
function wire(inputAdapter: InputAdapter, outputAdapter: OutputAdapter): void {
  const run = () => {
    const validation = inputAdapter.validate();
    if (!validation.valid) return;

    try {
      const inputs = inputAdapter.collectInputs();
      const outputs = compute(inputs);
      outputAdapter.update(outputs);
    } catch (err) {
      if (err instanceof Error)
        grok.shell.warning(err.message);
    }
  };

  inputAdapter.onChanged(() => run());

  // Initial run after UI rendering
  setTimeout(() => run(), 10);
}
```

### Wiring with conditional logic

```typescript
function wireWithAutoTune(
  inputAdapter: InputAdapter,
  outputAdapter: OutputAdapter,
  autoTuneEnabled: () => boolean,
  optimize: (inputs: ComputeInputs) => Promise<Partial<ComputeInputs>>,
): void {
  const run = async () => {
    const validation = inputAdapter.validate();
    if (!validation.valid) return;

    try {
      let inputs = inputAdapter.collectInputs();

      if (autoTuneEnabled())
        inputs = { ...inputs, ...(await optimize(inputs)) };

      const outputs = compute(inputs);
      outputAdapter.update(outputs);
    } catch (err) {
      if (err instanceof Error)
        grok.shell.warning(err.message);
    }
  };

  inputAdapter.onChanged(() => run());
}
```

### Full assembly

```typescript
function runMyMethod(): void {
  const df = grok.shell.t;
  if (df == null) { grok.shell.warning('No table'); return; }

  const view = grok.shell.getTableView(df.name);

  // 1. Create adapters
  const inputAdapter = createInputAdapter(df);
  const outputAdapter = createOutputAdapter(df);

  // 2. Wire them together
  wire(inputAdapter, outputAdapter);

  // 3. Place in UI (this is a separate responsibility — layout)
  view.dockManager.dock(inputAdapter.form, DG.DOCK_TYPE.LEFT, null, undefined, 0.15);
  // ... dock viewers ...
}
```

---

## Reference: real-world examples

| Component | File |
|-----------|------|
| Full UI (form + viewers + state) | `probabilistic-scoring/prob-scoring.ts` |
| Simple dialog with validation | `anova/anova-ui.ts` |
| Conditional visibility, async dialog | `missing-values-imputation/ui.ts` |

> **Note:** Existing code is written in Reactive Form style (direct callbacks without abstractions).
> New code should follow the Ports and Adapters pattern described in this guide.
