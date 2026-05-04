---
name: datagrok-exec
description: Use whenever you need the Datagrok browser to actually execute JavaScript — adding viewers, filtering, modifying the view, or returning a result widget to the chat. Plain javascript blocks do NOT execute. Open this skill before emitting a datagrok-exec block.
---

# datagrok-exec

The ONLY way code runs in the Datagrok browser is inside a fenced block tagged
`datagrok-exec`. Regular ` ```javascript ` blocks render as snippets only.

## Globals inside the block

| Variable | Type            | Available                                         |
|----------|-----------------|---------------------------------------------------|
| `grok`   | module          | Always                                            |
| `ui`     | module          | Always                                            |
| `DG`     | module          | Always                                            |
| `view`   | `DG.ViewBase`   | Always — check `view.type` for specific view type |
| `t`      | `DG.DataFrame`  | Only when `view.type === 'TableView'`             |

`grok`, `ui`, `DG` are already in scope — do not import them. The block runs in
an async IIFE, so `await` works directly.

## Multiple blocks in one response

Each `datagrok-exec` block runs in its own JavaScript scope (a fresh
`new Function(...)` IIFE). When you emit multiple blocks in a single response,
they execute sequentially — block N+1 awaits block N's promise — but JS
variables declared in earlier blocks are NOT visible in later ones.

State that **does** persist across blocks:
- The `view` object reference (same `DG.ViewBase` each time)
- The `t` DataFrame reference — column additions, filter changes, and any
  other mutations are visible to later blocks
- Anything you push into the platform itself: viewers added to the view,
  dialogs opened, columns appended, server-side state

State that does **not** persist:
- `const` / `let` / `var` bindings
- Local helper functions
- Cached values (e.g. `const ic50Col = t.col('IC50')`)

So if block 1 adds an `IC50` column and block 2 needs to read it, block 2
must re-derive `t.col('IC50')` from `t`. Never reference a variable from an
earlier block — that's a `ReferenceError`, not a "column doesn't exist yet"
problem.

When to split into multiple blocks vs. one big block:

- **Multi-block** is best when each step is independently meaningful to the
  user — they want to see step 1's result land before step 2 starts.
- **Single block** is best when steps share complex local state (intermediate
  arrays, helper functions, derived columns referenced multiple times) that
  would be tedious to re-derive every block.

Before consuming a query/function result inside `datagrok-exec`, look at the
wrapper's return type in the apiRef. Example — if the wrapper signature is
`Promise<string>`:

```js
// CORRECT — treat the awaited value as a string
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
| molecule (SMILES / molblock)      | `grok.chem.drawMolecule(smiles, 300, 200)`     |
| macromolecule (HELM)              | see [HELM output](#helm-output) below          |
| graphics                          | see [Graphics output](#graphics-output) below  |

## HELM output

For HELM macromolecule notation, render via the async HELM input.

```datagrok-exec
const helmInput = await ui.input.helmAsync('', { editable: false });
helmInput.stringValue = helmString;
return helmInput.root;
```

## Graphics output

Some Datagrok scripts return graphical elements. When you return such an element directly, it may render as blank because it lacks intrinsic dimensions. To display it correctly, extract the Base64-encoded image data and wrap it in an image component with explicit width and height:

```datagrok-exec
const b64 = await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', {mol: 'CCC', contours: 10});
return ui.image('data:image/png;base64,' + b64, 300, 200);
```