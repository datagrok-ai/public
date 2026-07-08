---
name: datagrok-exec
description: Use whenever you need the Datagrok browser to actually execute JavaScript — adding viewers, filtering, modifying the view, or returning a result widget to the chat. Open this skill before calling the datagrok_exec tool.
---

# datagrok-exec

The ONLY way code runs in the Datagrok browser is via the **`datagrok_exec` tool**.
Regular markdown code blocks do NOT execute.

Call `datagrok_exec` to *perform* an action the user asked for (add a viewer, filter,
open a file, run a function). For informational questions — "how do I…", "what
is…", "explain…", "can you…" — answer in plain text; do not call the tool.

## Globals available inside the code

| Variable | Type            | Available                                         |
|----------|-----------------|---------------------------------------------------|
| `grok`   | module          | Always                                            |
| `ui`     | module          | Always                                            |
| `DG`     | module          | Always                                            |
| `view`   | `DG.ViewBase`   | Always — check `view.type` for specific view type |
| `t`      | `DG.DataFrame`  | Only when `view.type === 'TableView'`             |

The code runs in an async IIFE, so `await` works directly.

## Tool result

The tool returns `{success, returnValue?, error?}`.

**`success` only means the JS did not throw — NOT that the intended effect happened.** A typo'd
viewer name, a no-op, or an action that ran twice all return `success: true`. `returnValue` is
whatever object *your own code* chose to `return`; it populates the chat confirmation shown to the
user, but it is **not proof** — it is self-reported by the same code that took the action. Make your
code `return` a plain confirmation object so the user sees accurate details:

| Action              | Return                              |
|---------------------|-------------------------------------|
| open / load table   | `{name, rowCount}`                  |
| add viewer          | `{type}`                            |
| filter              | `{filteredRowCount, totalRowCount}` |
| upload DataFrame    | `{tableInfoId}`                     |
| save layout         | `{layoutId}`                        |
| open app / view     | `{viewName, viewType}`              |
| color / sort / pin  | `{column}`                          |

**Error**: `{success: false, error: "..."}` — report verbatim, do not silently retry.

**Never pre-announce** success before calling the tool.

## Verifying actions — required after every action

Proof does not come from your own action code — it comes from **re-reading live state in a separate
call**. After any action (`datagrok_exec` that changed something, or any MCP call that changed
platform state), call the **`datagrok_verify`** tool before reporting success. It runs your assertion in a
fresh scope (globals `grok`, `ui`, `DG`, `view`, `t`) that **cannot see your action's variables**, so
it can only pass by re-deriving from the real platform state.

The assertion is yours to design for whatever you changed: re-read that thing and `return` the
observed result (truthy = verified).
`{passed: true}` → report the observed value, not the value you intended. `{passed: false}` → the
action did NOT work; fix it and verify again — do not report success. The runtime blocks the turn
from ending until a verify passes (bounded retries — if they run out, the response is shown to the
user flagged as **unverified**, so never present unverified work as done).

## Returning a result to the chat

Only return an `HTMLElement` when the **explicit user goal** is to see output in the chat
("show me the molecule", "display the table", "draw a plot"). Never return HTMLElements for
intermediate results, confirmation data, or actions where the platform UI is the destination.
An HTMLElement return replaces the confirmation object — no `returnValue` is sent.

Raw scalars, strings, or `DG.DataFrame` will NOT render; convert via the table below:

| Result type                       | How to render                                  |
|-----------------------------------|------------------------------------------------|
| scalar                            | `ui.divText(String(value))`                    |
| key-value pairs                   | `ui.tableFromMap({key: value})`                |
| list of items                     | `ui.divV(items.map((x) => ui.divText(x)))`     |
| `DG.DataFrame`                    | `DG.Viewer.grid(df).root`                      |
| `DG.Viewer` / `DG.Widget`         | `obj.root`                                     |
| `DG.ViewBase` (incl. apps)        | open via `grok.shell.addView(v)` — do NOT return |
| molecule (SMILES / molblock)      | `grok.chem.drawMolecule(smiles, 300, 200)`     |
| macromolecule (HELM)              | see [HELM output](#helm-output) below          |
| graphics                          | see [Graphics output](#graphics-output) below  |

## Native top-menu commands

If an operation already exists as a top-menu command (aggregate, join, cluster,
add column, …), invoke it instead of building a custom `ui.dialog`. Walk the menu
with `find()` (chain through `|` groups), matching leaf text exactly:

```js
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

For task-specific API surface, open the matching skill before writing code:

| User intent                                    | Skill                          |
|------------------------------------------------|--------------------------------|
| Find / describe / add / remove / rename columns; set semType/units/format/friendly name; color coding | `datagrok-df-and-columns`   |
| Calculated (formula-driven) columns            | `datagrok-calc-column`         |
| Filter rows (range, categorical, contains, substructure, predicate) | `datagrok-filtering`     |
| Select rows; current row; selection ↔ filter   | `datagrok-selection`           |
| Add / configure / find / close viewers         | `datagrok-viewers`             |
| Grid sort, visibility, widths, pins, freeze; grid-only cell tint | `datagrok-grid-customization` |
| Cheminformatics: SMILES/MolBlock/InChI/canonicalize | `datagrok-chem-data` / `datagrok-chem-toolkit` |

## Multiple tool calls in one response

Each `datagrok_exec` call runs in its own scope (a fresh `new Function(...)` IIFE).
Calls execute sequentially, but JS variables declared in one call are NOT visible in the next.

Persists across calls:
- The `view` reference and the `t` DataFrame reference (column additions,
  filter changes, and other mutations are visible to later calls)
- Anything pushed into the platform: viewers added, dialogs opened, columns
  appended, server-side state

Does NOT persist: `const`/`let`/`var` bindings, helper functions, cached values.
A second call must re-derive from `t`.

Before consuming a query/function result, check the wrapper's return type in the apiRef:

```js
// CORRECT — wrapper returns Promise<string>
const helm = await grok.data.query('Biologics:GetBiologicsPeptideHelmByIdentifier',
  { peptideIdentifier: 'GROKPEP-000002' });

// WRONG — assuming a DataFrame because most queries return one
const df = await grok.data.query('Biologics:GetBiologicsPeptideHelmByIdentifier', {...});
const helm = df.columns.byIndex(0).get(0);
```

## Launching apps & views

Datagrok **apps** are functions with `meta.role: app`. Calling one via
`grok.functions.call('Pkg:appName')` may return a `DG.ViewBase` (or nothing if
the app opens its own view). A returned view is NOT yet attached to the
workspace — hand it to `grok.shell.addView()`. Do NOT `return` the view from
the code; open it.

```js
const result = await grok.functions.call('Pkg:appName');
if (result instanceof DG.ViewBase)
  grok.shell.addView(result);
```

## HELM output

```js
const helmInput = await ui.input.helmAsync('', { editable: false });
helmInput.stringValue = helmString;
return helmInput.root;
```

## Graphics output

Returned graphics elements may render blank without intrinsic dimensions.
Extract Base64-encoded image data and wrap with explicit width and height:

```js
const b64 = await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', {mol: 'CCC', contours: 10});
return ui.image('data:image/png;base64,' + b64, 300, 200);
```
