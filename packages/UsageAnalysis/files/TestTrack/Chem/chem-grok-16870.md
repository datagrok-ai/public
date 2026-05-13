---
feature: chem
sub_features_covered: [chem.rendering, chem.rendering.molecule-cell, chem.rendering.rdkit-renderer]
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
produced_from: atlas-driven
produced_for: chem-grok-16870-spec.ts
authored_date: 2026-05-11
authored_by: claude-code-test-designer
related_bugs: [GROK-16870]
---

# Chem | RDKit cell renderer in non-Chem viewer tooltip (GROK-16870 regression-lock)

Bug-focused scenario locking the invariant: **the RDKit cell renderer MUST NOT
throw on null DataFrame when invoked from a non-Chem viewer tooltip rendering
path** (Box Plot, Scatter Plot, Histogram, Grid). Hovering over a data point of
a non-Chem viewer on a table containing a molecule column must not surface
`NullError: method not found: 'gS' on null` from `rdkit-cell-renderer.ts`. A
clean tooltip is the post-fix expected behavior; a console-error spew on hover
is the regression this scenario catches.

Per chain YAML (`scenario-chains/chem.yaml` rev 2
`bug_focused_candidates[chem-grok-16870-spec.ts]`): independent scenario
(`depends_on: []`), `pyramid_layer: bug-focused`, `target_layer: playwright`,
strategy `simple`. Parallel-coverage on `chem.rendering.*` with `info-panels.md`
(Context Panel rendering across notation formats, `coverage_type: regression`).
This spec owns the cross-viewer-tooltip rendering invariant — exercises a code
path NOT covered by `info-panels.md` (which is Context Panel only).

Bug-library reference: `references/bug-library/chem.yaml :: GROK-16870` —
title *Chem: Fix rdkit cell renderer error in scatter plot tooltip*, priority
p2, status `fixed`, `fixed_in: 1.22.0`, `test_coverage: needed`. Atlas critical
path overlap: `chem.cp.molecule-rendering-end-to-end` (p0 — RDKit renderer
must not crash in any rendering context).

## Setup

1. **No external provisioning.** The scenario consumes the bundled chem demo
   dataset `System:AppData/Chem/tests/smiles-50.csv` (validated in MCP recon 2026-05-11
   — 50 rows, `canonical_smiles` column auto-detects `semType: Molecule`,
   numeric columns also present for Box Plot value axis). No bespoke fixture,
   no project save, no share.
2. **Login** — standard test user. No special privileges required.
3. **No fixture consumed.** Independent scenario (`depends_on: []`).

## Scenarios

### Block A — Hover over Box Plot data point with a molecule column on the table — no rdkit-cell-renderer errors

Locks the invariant: hovering over data points in a non-Chem viewer (Box Plot
here as the canonical bug-report repro vector) on a table that contains a
`semType: Molecule` column MUST NOT surface RDKit cell renderer crashes
(`NullError`, `method not found`, `cellRenderer.render` thrown errors). The
viewer must remain responsive and the tooltip (if rendered) must show molecule
content cleanly OR no tooltip without crash.

1. **Open the chem demo dataset and add a Box Plot viewer.** Read
   `System:AppData/Chem/tests/smiles-50.csv`, add a table view, wait for semantic type
   detection (`canonical_smiles` becomes `semType: Molecule`), add a Box Plot
   viewer via `grok.shell.tv.addViewer('Box plot')`. The Box Plot renders against
   the default numeric column. Box Plot DOM container resolves at
   `[name="viewer-Box-plot"]`.
2. **Hook console.error to capture rdkit cell renderer errors.** Replace
   `console.error` with a function that pushes the stringified arguments to
   `window.__grok16870_errors` before delegating to the original handler.
   This captures the bug's stack signature
   (`NullError`, `method not found: 'gS' on null`, `rdkit-cell-renderer`,
   `cellRenderer.render`).
3. **Hover over the Box Plot at multiple positions.** Use Playwright's
   `page.mouse.move()` (real pointer events; programmatic
   `dispatchEvent('mouseover')` does not consistently trigger the tooltip
   pipeline in viewer canvases). Walk a 5-position sweep across the central
   region of the Box Plot bounding rect to maximize the probability of hitting
   a data point category and surfacing the tooltip render path.
4. **Wait the post-hover tooltip render window.** Allow up to 6 seconds of
   settle time after each hover position. The tooltip is async — the RDKit
   renderer fires on a microtask after the hover event.
5. **Assert no rdkit cell renderer errors.** Filter the captured
   `__grok16870_errors` array for entries matching
   `/rdkit-cell-renderer|method not found|gS|cellRenderer\.render|NullError/i`.
   The expected count is 0. Any matching error is the GROK-16870 regression.
6. **Assert Box Plot remains responsive.** Verify the Box Plot DOM container
   is still in the document, has non-zero size, and that its canvas elements
   are still mounted. This catches a secondary failure mode where a renderer
   crash poisons the viewer's render loop.

## Notes

- **`coverage_type: edge`** — non-Chem-viewer tooltip rendering edge case;
  this scenario asserts the failure-mode invariant for
  `chem.rendering.rdkit-renderer` outside its canonical grid-cell render
  context. Happy-path Context Panel rendering of molecule cells is owned by
  the parallel `info-panels.md` scenario (`coverage_type: regression`,
  `pyramid_layer: integration`).
- **Viewer choice — Box Plot as canonical repro.** The bug report names Box
  Plot as the first-observed crash vector ("Add a Box Plot viewer / Hover
  over a data point"). The bug edge_case_for_atlas notes the recurring class
  is "Box Plot, Scatter Plot, Grid, Histogram, etc.". This spec scopes to Box
  Plot only — single-viewer-vector single-invariant scoped per
  `strategy: simple`. A future cross-viewer cross-cutting spec could extend
  the same invariant to Scatter Plot / Histogram if a related regression
  resurfaces; under the current chain `bug_focused_candidates[GROK-16870]`
  span (`info-panels.md:Step 4`) Box Plot is the in-scope choice.
- **Real-pointer-event requirement.** The hover-triggered tooltip render in
  Datagrok viewer canvases is driven by the platform's internal pointer event
  handlers. Synthetic `dispatchEvent('mousemove' / 'mouseover')` events on
  DOM nodes do NOT reliably trigger the canvas hover-tooltip pipeline (MCP
  recon 2026-05-11: 5-position synthetic sweep produced 0 tooltips, 0
  errors). Playwright's `page.mouse.move()` produces real Chrome pointer
  events through the DevTools protocol and reliably engages the tooltip
  pipeline.
- **Lenient tooltip-presence assertion.** The spec asserts ABSENCE of crashes
  (the bug signature) rather than PRESENCE of a specific tooltip (which
  depends on data-density at the hovered position — and `smiles.csv` Box
  Plot may render an empty central region depending on the chosen
  category/value axes). A successful test = clean console + responsive
  viewer; tooltip-rendered-with-molecule-content is a stronger claim left
  to `info-panels.md` (cell render via Context Panel) and a future
  `chem.rendering.molecule-tooltip` cross-viewer spec if needed.
- **MCP-validated current behavior (2026-05-11 dev.datagrok.ai recon, qa-pw
  user, Chem package v1.17.6).** Box Plot rendering with smiles.csv produces
  no `rdkit-cell-renderer` console errors on initial render. GROK-16870 is
  `status: fixed`, `fixed_in: 1.22.0` per bug-library; spec exercises the
  post-fix invariant as a regression-lock. No `test.fixme()` — enabled.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` from
  `helpers-registry.yaml`. Inline patterns for dataset open + viewer add +
  hover sweep — no new helper authored under reuse threshold (single-use
  pattern; helper authoring threshold ≥3 not reached).
- **No JS API substitution for hover.** Steps 3-4 exercise the hover-tooltip
  pipeline that the bug lives in. JS API direct invocation of
  `rdKitCellRenderer.renderInternal` with a null DataFrame argument would
  reproduce the underlying API-level crash but NOT the tooltip-context-
  invocation surface that the bug-library entry describes
  (`Tooltip.getRowTooltip` → `GridCellRendererJsProxy.render` →
  `rdKitCellRenderer.renderInternal`). Real hover is the actionable path.
- **Order in chain.** Bug-focused candidate; placed parallel to
  `info-panels.md`. No `depends_on:` ordering — independent scenario.
