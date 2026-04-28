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
| graphics                          | see [Graphics output](#graphics-output) below  |

## Graphics output

Some Datagrok scripts return graphical elements. When you return such an element directly, it may render as blank because it lacks intrinsic dimensions. To display it correctly, extract the Base64-encoded image data and wrap it in an image component with explicit width and height:

```datagrok-exec
const b64 = await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', {mol: 'CCC', contours: 10});
return ui.image('data:image/png;base64,' + b64, 300, 200);
```