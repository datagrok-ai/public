---
name: datagrok-exec
description: Use whenever you need the Datagrok browser to actually execute JavaScript — adding viewers, filtering, modifying the view, or returning a result widget to the chat. Plain javascript blocks do NOT execute. Open this skill before emitting a datagrok-exec block.
---

# datagrok-exec

The ONLY way code runs in the Datagrok browser is inside a fenced block tagged
`datagrok-exec`. Regular ` ```javascript ` blocks render as snippets only.

Emit a block to *perform* an action the user asked for (add a viewer, filter,
open a file, run a function). For informational questions — "how do I…", "what
is…", "explain…", "can you…" — answer in plain text; do not execute code.

## Globals inside the block

| Variable | Type            | Available                                         |
|----------|-----------------|---------------------------------------------------|
| `grok`   | module          | Always                                            |
| `ui`     | module          | Always                                            |
| `DG`     | module          | Always                                            |
| `view`   | `DG.ViewBase`   | Always — check `view.type` for specific view type |
| `t`      | `DG.DataFrame`  | Only when `view.type === 'TableView'`             |

The block runs in an async IIFE, so `await` works directly.

## Native top-menu commands

If an operation already exists as a top-menu command (aggregate, join, cluster,
add column, …), invoke it instead of building a custom `ui.dialog`. Walk the menu
with `find()` (chain through `|` groups), matching leaf text exactly:

```datagrok-exec
grok.shell.topMenu.find('Data').find('Aggregate Rows...').click();
grok.shell.topMenu.find('Edit').find('Go To').find('Row...').click();
```

| Menu     | Commands                                                                                                                  |
|----------|--------------------------------------------------------------------------------------------------------------------------|
| `Edit`   | `Undo`, `Redo`, `Add New Column...`, `Add Rows...`, `Column Properties...`, `Find and Replace...`, `Go To \| Row...`, `Go To \| Next Selected`, `Go To \| Prev Selected`, `Remove \| Selected Rows`, `Remove \| Selected Columns`, `Remove \| Selected Rows or Columns` |
| `View`   | `Columns`, `Console`, `Context Panel`, `Schema`, `Tables`, `Variables`, `Toolbox`, `Windows`, `Home`, `Reset Filter`, `Edit Tooltip...`, `Embed...`, `Full Screen`, `Presentation Mode`, `Layout \| Save to Gallery`, `Layout \| Open Gallery`, `Layout \| Clone View`, `Layout \| Download`, `Layout \| Clear` |
| `Select` | `All`, `None`, `Invert`, `All Columns`, `Duplicates...`, `Missing Values...`, `Random Rows...`, `Extract Selected Rows`, `Filter to Column...`, `Selection to Column...`, `Selection to Filter` |
| `Data`   | `Aggregate Rows...` (also pivot), `Join Tables...`, `Link Tables...`, `Compare Tables...`, `Compare Columns...`, `Append Tables...`, `Unpivot...`, `Split Column...`, `Split Column by RegExp...`, `Categorize...`, `Anonymize...`, `Batch Edit...` |
| `ML`     | `Analyze \| PCA...`, `Analyze \| PLS...`, `Analyze \| ANOVA...`, `Analyze \| Multivariate Analysis...`, `Cluster \| DBSCAN...`, `Cluster \| MCL...`, `Models \| Train Model...`, `Models \| Apply Model...`, `Notebooks \| New Notebook...`, `Reduce Dimensionality...`, `Impute Missing Values...`, `Pareto Front...` |

## Per-area skills

For task-specific API surface, open the matching skill before emitting code:

| User intent                                    | Skill                          |
|------------------------------------------------|--------------------------------|
| Find / describe / add / remove / rename columns; set semType/units/format/friendly name; color coding | `datagrok-df-and-columns`   |
| Calculated (formula-driven) columns            | `datagrok-calc-column`         |
| Filter rows (range, categorical, contains, substructure, predicate) | `datagrok-filtering`     |
| Select rows; current row; selection ↔ filter   | `datagrok-selection`           |
| Add / configure / find / close viewers         | `datagrok-viewers`             |
| Grid sort, visibility, widths, pins, color coding, freeze | `datagrok-grid-customization` |
| Cheminformatics: SMILES/MolBlock/InChI/canonicalize | `datagrok-chem-data` / `datagrok-chem-toolkit` |

## Multiple blocks in one response

Each `datagrok-exec` block runs in its own scope (a fresh `new Function(...)`
IIFE). Blocks execute sequentially — block N+1 awaits block N — but JS
variables declared in earlier blocks are NOT visible in later ones.

Persists across blocks:
- The `view` reference and the `t` DataFrame reference (column additions,
  filter changes, and other mutations are visible to later blocks)
- Anything pushed into the platform: viewers added, dialogs opened, columns
  appended, server-side state

Does NOT persist: `const` / `let` / `var` bindings, helper functions, cached
values like `const ic50Col = t.col('IC50')`. Block 2 must re-derive from `t`.

Use multiple blocks when each step is independently meaningful to the user.
Use a single block when steps share complex local state (intermediate arrays,
helper functions, derived columns referenced multiple times).

Before consuming a query/function result, check the wrapper's return type in
the apiRef:

```js
// CORRECT — wrapper returns Promise<string>
const helm = await grok.data.query('Biologics:GetBiologicsPeptideHelmByIdentifier',
  { peptideIdentifier: 'GROKPEP-000002' });

// WRONG — assuming a DataFrame because most queries return one
const df = await grok.data.query('Biologics:GetBiologicsPeptideHelmByIdentifier', {...});
const helm = df.columns.byIndex(0).get(0);
```

## Returning a result to the chat

If the code modifies the view directly (adds a viewer, filters, color-codes,
appends columns), do NOT return anything.

If the code produces a result to *show* the user, `return` an `HTMLElement` —
it will be appended to the chat. Raw scalars, strings, or `DG.DataFrame` will
NOT render; convert via the table below:

| Result type                       | How to render                                  |
|-----------------------------------|------------------------------------------------|
| scalar                            | `ui.divText(String(value))`                    |
| key-value pairs                   | `ui.tableFromMap({key: value})`                |
| list of items                     | `ui.divV(items.map((x) => ui.divText(x)))`     |
| `DG.DataFrame`                    | `DG.Viewer.grid(df).root`                      |
| `DG.Viewer` / `DG.Widget`         | `obj.root`                                     |
| `DG.ViewBase` (incl. apps)        | open via `grok.shell.addView(v)` — see [Launching apps & views](#launching-apps--views) — do NOT return |
| molecule (SMILES / molblock)      | `grok.chem.drawMolecule(smiles, 300, 200)`     |
| macromolecule (HELM)              | see [HELM output](#helm-output) below          |
| graphics                          | see [Graphics output](#graphics-output) below  |

## Launching apps & views

Datagrok **apps** are functions with `meta.role: app`. Calling one via
`grok.functions.call('Pkg:appName')` may return a `DG.ViewBase` (or nothing if
the app opens its own view). A returned view is NOT yet attached to the
workspace — hand it to `grok.shell.addView()`. Do NOT `return` the view from
the block; open it.

```datagrok-exec
const result = await grok.functions.call('Pkg:appName');
if (result instanceof DG.ViewBase)
  grok.shell.addView(result);
// If the app already added its own view (returns null/undefined), do nothing.
```

Never assign to `grok.shell.v` to launch an app — that just swaps the current
view reference without registering it.

## HELM output

For HELM macromolecule notation, render via the async HELM input.

```datagrok-exec
const helmInput = await ui.input.helmAsync('', { editable: false });
helmInput.stringValue = helmString;
return helmInput.root;
```

## Graphics output

Returned graphics elements may render blank without intrinsic dimensions.
Extract Base64-encoded image data and wrap with explicit width and height:

```datagrok-exec
const b64 = await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', {mol: 'CCC', contours: 10});
return ui.image('data:image/png;base64,' + b64, 300, 200);
```
