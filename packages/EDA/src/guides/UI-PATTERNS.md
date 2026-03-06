# UI Patterns for EDA Methods

Reference for implementing user-facing dialogs and result visualization in the EDA package.
For computation patterns, see `COMPUTATION-PATTERNS.md`.
For worker-based methods, see `WORKER-GUIDE.md`.

## Dialog Pattern

Standard flow: dialog with validated inputs -> computation -> dock results into the current table view.

```typescript
export function runMyMethod(): void {
  const df = grok.shell.t;
  if (df === null) { grok.shell.warning('No dataframe is opened'); return; }

  // 1. Build inputs with validation
  const input = ui.input.column('Feature', {
    table: df,
    filter: (col: DG.Column) => [DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.INT].includes(col.type),
    nullable: false,
  });

  // 2. Create dialog
  const dlg = ui.dialog({title: 'My Method', helpUrl: '/help/...'});
  const view = grok.shell.getTableView(df.name);
  view.root.appendChild(dlg.root);

  dlg.addButton('Run', () => {
    dlg.close();
    try {
      const result = computeMyMethod(/* args */);
      // 3. Dock results (see Docking section below)
    } catch (error) {
      if (error instanceof Error)
        grok.shell.warning(error.message);
    }
  });

  dlg.add(input).show();
}
```

Reference: `anova/anova-ui.ts` (`runOneWayAnova`), `missing-values-imputation/ui.ts` (`runKNNImputer`).

---

## Dynamic Input Validation

Validate inputs on change and enable/disable the run button accordingly:

```typescript
const dlg = ui.dialog({title: 'My Method'});

const onInputChanged = () => {
  const isValid = /* validation logic */;
  dlg.getButton('Run').disabled = !isValid;
};

const input = ui.input.column('Feature', {
  table: df,
  onValueChanged: () => onInputChanged(),
});
```

For complex cases, show/hide additional inputs based on current state:

```typescript
const showWidgets = () => {
  dlg.getButton('Run').disabled = false;
  advancedInput.root.hidden = false;
};

const hideWidgets = () => {
  dlg.getButton('Run').disabled = true;
  advancedInput.root.hidden = true;
};
```

Reference: `missing-values-imputation/ui.ts` (`checkApplicability`, `showWidgets`, `hideWidgets`).

---

## Async Dialog Pattern

Return a `Promise<void>` for callers that need to await dialog completion:

```typescript
export async function runMyMethod(): Promise<void> {
  let resolve: (value: void | PromiseLike<void>) => void;
  let reject: (reason?: any) => void;
  let okClicked = false;

  const promise = new Promise<void>((res, rej) => {
    resolve = res;
    reject = rej;
  });

  const dlg = ui.dialog({title: 'My Method'});

  dlg.addButton('Run', () => {
    okClicked = true;
    dlg.close();
    try {
      // ... computation ...
      resolve();
    } catch (err) {
      reject(err);
    }
  });

  dlg.show().onClose.subscribe(() => !okClicked && resolve());

  return promise;
}
```

Reference: `missing-values-imputation/ui.ts` (lines 251-290).

---

## Docking Results

Dock viewers and reports into the current table view:

```typescript
const view = grok.shell.getTableView(df.name);
grok.shell.v = view;

// Dock a chart
const chart = DG.Viewer.boxPlot(df, {
  categoryColumnNames: [factorName],
  valueColumnName: featureName,
});
const node = view.dockManager.dock(chart, DG.DOCK_TYPE.RIGHT, null, 'Results');

// Dock a report below the chart
const report = ui.tabControl({
  'Analysis': ui.panel([grid.root]),
  'Details': ui.panel([detailsDiv]),
});
view.dockManager.dock(report.root, DG.DOCK_TYPE.DOWN, node, '', 0.25);
```

Reference: `anova/anova-ui.ts` (`addVizualization`).

---

## Result Grid with Tooltips

Display tabular results with column header tooltips:

```typescript
const grid = DG.Viewer.grid(DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Source', ['Between', 'Within', 'Total']),
  DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Value', [1.23, 4.56, 5.79]),
]));

const tooltips = new Map([
  ['Source', 'Description of sources.'],
  ['Value', 'Description of values.'],
]);

grid.onCellTooltip((cell, x, y) => {
  if (cell.isColHeader) {
    ui.tooltip.show(ui.divV([ui.p(tooltips.get(cell.tableColumn!.name)!)]), x, y);
    return true;
  }
});
```

Reference: `anova/anova-ui.ts` (`getAnovaGrid`).
