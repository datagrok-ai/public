---
feature: chem
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
priority: p1
realizes: [GROK-16870]
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

This spec runs parallel to `info-panels.md` (which covers Context Panel
rendering across notation formats) — it owns specifically the cross-viewer
tooltip-rendering invariant, a code path `info-panels.md` does not touch
(Context Panel only).

Bug reference: GROK-16870 — *Chem: Fix rdkit cell renderer error in scatter
plot tooltip* (priority p2, fixed in 1.22.0).

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
   `System:AppData/Chem/tests/smiles-50.csv`, add a table view, wait for
   semantic type detection (`canonical_smiles` becomes `semType: Molecule`),
   add a Box Plot viewer via `grok.shell.tv.addViewer('Box plot')`. The Box Plot renders against
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

- **Why Box Plot.** The bug report names Box Plot as the first-observed crash vector, though
  the same renderer class is shared by Scatter Plot, Histogram, and Grid tooltips. This
  scenario scopes to Box Plot only; a future cross-viewer spec could extend the same
  invariant to the other viewers if a related regression resurfaces.
- **Verification scope.** The scenario asserts absence of a console-error crash rather than
  presence of a specific tooltip — tooltip content depends on where the hover lands, which
  varies with the data. A clean console plus a still-responsive viewer counts as a pass;
  asserting that a tooltip actually renders molecule content is covered separately by
  `info-panels.md` (Context Panel).
