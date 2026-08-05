---
feature: network-diagram
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - network-diagram-spec.ts
related_bugs:
  - id: GROK-20617
    status: open
  - id: GROK-20618
    status: open
  - id: GROK-17125
    status: open
---

# Network diagram (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **Network diagram** from **Toolbox > Viewers**

Everything below is driven through the UI: the Viewers toolbox, the on-viewer node
selectors, the Context Panel property grid, and real mouse input on the canvas.
The JS API is used only to read back the selected- and filtered-row counts.

## Add the viewer

1. Click the **Network diagram** icon in **Toolbox > Viewers**
2. The on-viewer selectors read **SEX** and **CONTROL** — the node columns are
   picked automatically
3. The diagram is drawn (the canvas carries content)

## Suspend simulation

1. Go to **Context Panel > Misc** and check **Suspend Simulation**
2. The layout freezes — the canvas stops repainting on its own

## Node columns

1. Set **Node 1** to *RACE* with the on-viewer selector
2. The selector shows *RACE* and the diagram is redrawn

## Colour and size coding

1. On **Context Panel > Data**, set **Node1 Color** to *SEX* and **Node1 Size** to *AGE*
2. The property grid shows both columns and the diagram repaints
3. Set **Edge Color** to *AGE* — the edges are recoloured

## Clicking

1. Click a node — the rows behind it are selected and the diagram repaints
2. Uncheck **Select Rows On Click** *and* **Select Edges On Click** in **Misc**
3. Clicking the same places selects nothing
4. Re-check both

## Column selectors

1. Uncheck **Show Column Selectors** in **Misc** — the on-viewer selectors disappear
2. Check it again — they come back

## Arrows

1. Set **Show Arrows** to *to* in **Misc** — the arrow heads should be drawn
   straight away. They are not: the diagram only picks the setting up on the next
   rebuild (GROK-20617, guarded by `knownOpenBug`).
2. Change **Node 1** to force the graph to be rebuilt — the arrow heads are drawn
   and the setting survives the rebuild.

## Filtered out nodes

1. Filter the table to `SEX = F`
2. Check **Show Filtered Out Nodes** in **Misc** — the filtered-away nodes should
   come back. They do not, and nothing is repainted (GROK-20618, guarded by
   `knownOpenBug`).
3. Uncheck it and reset the filter

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`network-diagram-spec.ts`. Kept verbatim so no scenario is lost.

### Selection gestures on nodes and edges

> Manual

5. Shift+click another node, then Ctrl+click the first node.
   * Expected: Shift+click adds to selection; Ctrl+click toggles it.
6. Click an edge.
   * Expected: Only the rows backing that edge are selected.
7. Double-click empty canvas.
   * Expected: Selection clears on the first click; the view zooms to fit on the second.

---
{
  "order": 11,
  "datasets": ["System:DemoFiles/demog.csv"]
}
