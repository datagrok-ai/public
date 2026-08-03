---
feature: grid
realizes_atlas:
  - grid.cp.appearance-summary-persist
realizes:
  - viewers.grid
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-19769
    status: fixed
  - id: GROK-19809
    status: fixed
realized_as:
  - grid-appearance-summary-persist-spec.ts
expected_results:
  - anchor: Step 4
    expectation: Applying Linear through the AGE column header menu leaves Linear as
      the column's colour-coding type. The gradient VALUES are owned by the flat
      colour-coding section; this scenario re-reads all four colours once, in
      the persistence tail.
  - anchor: Step 5
    expectation: Applying Conditional through the HEIGHT column header menu leaves
      Conditional as the column's colour-coding type with the configured named
      ranges stored on it.
  - anchor: Step 6
    expectation: Applying Categorical through the SEX column header menu leaves
      Categorical as the column's colour-coding type.
  - anchor: Step 7
    expectation: Applying Linked through the WEIGHT header menu's Edit dialog, with
      SEX chosen as the source column, leaves Linked as the column's
      colour-coding type and SEX recorded as its source.
  - anchor: Step 9
    expectation: After dragging the row height slider to a non-default value, the
      height of the cell bounds (grid.cell('AGE', 0).bounds.height) reflects the
      new row height and is measurably different from the default. The
      row-height prop itself is NOT the assertion target.
  - anchor: Step 10
    expectation: After setting a non-default missing value color, grid.cell().color
      on a cell whose underlying value is null returns a color that matches the
      configured missing-value color and differs from the default background.
      The missingValueColor prop itself is NOT the assertion target.
  - anchor: Step 12
    expectation: After adding all one-click summary column types (Sparklines, Bar
      Chart, Pie Chart, Radar, Smart Form, Tags, Confidence Interval),
      grid.columns.length has increased by the count of types successfully
      added, and each new column's cellType matches the type that was chosen.
  - anchor: Step 13
    expectation: Adding the min and max stats rows leaves every summary column from
      the previous step in place and raises no console error.
  - anchor: Step 15
    expectation: 'After re-applying the saved layout, the set of open viewers does
      not include the foreign viewer added in Step 14. The four colour codings
      all read back correctly: the AGE gradient gives different colours at the
      lowest and highest values and both differ from a plain cell; each HEIGHT
      cell resolves the colour configured for its range while an out-of-range
      cell stays plain; two different SEX values give two different colours; and
      each WEIGHT cell matches the SEX cell of its own row. The cell bounds
      height (Step 9) and missing-value cell colour (Step 10) are restored, and
      every summary column added in Step 12 is still present with its cellType
      intact. This is the GROK-19769 guard: a layout carrying summary columns
      re-applies without error.'
  - anchor: Step 17
    expectation: After Close All and reopening the saved project, the same assertion
      battery from Step 15 passes in full, and the console-error delta is 0 (no
      new errors since the reopen).
---

# Grid — Appearance, Summary Columns, and Persistence

Apply color coding on four columns with four different modes, adjust row style settings
through the gear panel, add all one-click summary column types, add stats rows, then confirm
the whole configuration survives a layout re-apply and a project close-and-reopen.

## Setup

- Close all open views.
- Open demog. The table opens with its grid; every step below acts on that grid.
- Note the AGE value in the first row of the table and identify one null cell in the dataset
  (any column with missing values — HEIGHT typically has a few) to use in Step 10.

## Scenarios

### 1. Linear color coding on AGE

1. Click any cell in the grid to give it focus.
2. Find the AGE column. Identify which visible row holds the minimum AGE value and which
   holds the maximum AGE value at the current sort order.
3. Right-click the AGE column header and open **Color Coding > Linear**. The menu path is
   on the column header context menu; accept the default linear settings.
4. Read the per-cell resolved color on the min-AGE row and the max-AGE row via the product
   API. The min-row and max-row colors must differ from each other, and both must differ from
   the background color of an uncolored column cell in the same rows.

### 2. Conditional color coding on HEIGHT

5. Right-click the HEIGHT column header and open **Color Coding > Conditional**. In the dialog
   that opens, configure at least two explicit ranges with distinct named colors (for example,
   range 1: values below 160 → blue; range 2: values above 180 → red). Accept the dialog.
   Read the resolved color on cells that fall within each configured range. Cells in different
   ranges must report different colors; each must match the color configured for that range.

### 3. Categorical color coding on SEX

6. Right-click the SEX column header and open **Color Coding > Categorical**. Accept the
   auto-generated palette. Read the resolved color on the first row whose SEX value is "M" and
   the first row whose SEX value is "F". The two rows must report distinct colors.

### 4. Linked color coding on WEIGHT with source SEX

7. Right-click the WEIGHT column header and open **Color Coding > Linked**. In the **Edit...**
   dialog that appears, locate the "Source column" input and select **SEX** as the source.
   Accept the dialog. (Categorical color coding on SEX is a precondition — it must already
   be enabled, as set in Step 6.) Read the resolved color on the WEIGHT cell and the SEX cell
   of the same row for at least two rows with different SEX values. The WEIGHT cell color
   must equal the SEX cell color of the same row.

### 5. Style settings via the gear panel

8. Open grid properties by clicking the **gear icon** in the grid's top-right area (or via
   the grid context menu). The Grid Properties panel opens on the right.
9. In the Grid Properties panel, locate the **Row Height** control (a slider combined with a
   text box). Drag the slider or type a value to set a row height that is visibly different
   from the default (for example, 40 if the default is 25). Close the panel. Measure the
   bounds height of a data cell. The height must reflect the new value and differ from the
   original. Do not assert the rowHeight property directly — assert the cell bounds height.
10. Re-open the Grid Properties panel. Locate the **Missing Value Color** hex input and enter
    a clearly distinct color (for example, #FFAAAA for a light red). Close the panel.
    Find the null cell noted during setup and read its resolved color. The color must match
    the configured missing-value color and differ from the default background. Do not assert
    the missingValueColor property directly — assert the cell's resolved color.

### 6. Summary columns — all one-click types

11. Click on any data cell in the grid so the column context menu is accessible.
12. Right-click a column header, navigate to **Add > Summary Columns**, and add each of the
    one-click types in turn: **Sparklines**, **Bar Chart**, **Pie Chart**, **Radar**,
    **Smart Form**, **Tags**, **Confidence Interval**. (The form-designer entries in that menu
    are manual-only — skip them.) After each addition, the grid's column count increments by one
    and the newly added column's cell type matches the type chosen. Confirm all seven types are
    present by the end of this step.

### 7. Stats rows — no-error floor

13. Right-click any column header, navigate to **Add**, and add a **min** stats row, then a
    **max** stats row (menu items appear in lowercase). After both additions:
    adding the min and max stats rows leaves every summary column from the previous step in
    place and raises no console error. (Stats rows themselves have no product-readable signal
    and are not asserted beyond the summary-column-survival and no-error checks — this is an
    honest no-error floor, not an invented visibility claim. This is the GROK-19809 guard.)

### 8. The arrangement survives a layout and a project round-trip

14. Leave all changes above in place and save the current view layout from the ribbon, noting
    its name. Then add a foreign viewer such as a Scatter Plot so the view holds more than
    just the grid.
15. Re-apply the saved layout. The Scatter Plot added in Step 14 is gone. The grid comes back
    with all four color-coding configurations intact (read each per-cell color as in Steps 4-7),
    the cell bounds height from Step 9 restored, the missing-value cell color from Step 10
    restored, and all seven summary columns from Step 12 still present with matching cell types.
    No new console errors appear during re-apply. This is the GROK-19769 guard.
16. Save the view as a project with the ribbon **Save** button under a unique name, then
    Close All.
17. Reopen that project. The full assertion battery from Step 15 holds again, and no new
    console errors appeared since the reopen. This is the final persistence guard.
18. Afterwards, delete the probe layout and the probe project so nothing is left behind, even
    if an earlier step failed.

## Automation notes
- The save step is the ONE red left here, and the sibling scenario grid-columns-style-persist
  already solved exactly this — copy its shape rather than inventing another: it imports
  knownOpenBug from helpers/known-open-bug.ts, declares a letter-agnostic benign pattern for the
  publish-chain pair (the minified symbol drifts, so never pin 'aZ'/'aY'), and treats that class as
  noise ONLY inside the save window. A bare expect(balloons).toEqual([]) across the whole run is
  what keeps this step red.

- PERSISTENCE-TAIL CANON (corpus-verified; follow it instead of inventing a path):
  save with saveProjectViaUI(page, name) from helpers/projects.ts — the real ribbon Save, the only
  path that carries the view layout; reopen with grok.dapi.projects.find(id) then open() inside
  page.evaluate, which is what every persistence tail in the corpus does and which works (do NOT
  use openProjectFromDashboards here — the gallery tile is not reliably reachable); tear down with
  deleteProjectWithCleanup(page, {projectId}) in finally.
- ERROR-CHANNEL CANON for the save window: the ribbon Save renders a publication preview by cloning
  the live view into an offscreen iframe, and the detached view emits a "cloned iframe" message plus
  a Dart NullError PAIR. Sibling sections gate a benign filter to the save window only — a
  LETTER-AGNOSTIC pattern (the minified symbol drifts: 'aZ' here, 'aY' in the sibling scenario, so
  never pin the literal) — and treat the same class OUTSIDE that window as a regression signal.
  Adopt that shape; do not assert an empty balloon list across the whole run.
- If a strict assertion must stay red because the product is genuinely broken, wrap it in
  knownOpenBug(bugId, assertion) from helpers/known-open-bug.ts: the reproduction reports green and
  the assertion self-flips loudly when the bug is fixed. A permanently red expect() is not an
  acceptable resting state.
- resetShell(page) from helpers/openers.ts is the only helper that strips leftover BALLOONS from the
  DOM; error balloons never auto-hide and closeAll does not remove them, so call it before any phase
  whose balloon or console channel you intend to assert.


- Step 16 is RED against live dev and it is a PRODUCT defect, not a spec defect: the ribbon Save
  raises an error balloon ("NullError: method not found: 'aY' on null") on every project save.
  Reproducible by hand: open demog, ribbon Save, name it, OK. Do not weaken the assertion — the
  save path is expected to be silent. Everything else, including the full appearance battery after
  the project round-trip, passes.
- The gear property panel can be rendered while the Context Panel is collapsed or scrolled out of
  the headless viewport, so its rows are waited for as ATTACHED, not visible, and the colour row's
  click is dispatched on the node. The panel is opened by trying four real gestures in order
  (trusted hover+click on the gear, the registered viewer-gear helper, a synthetic gear click,
  F4 on the focused grid) — the first that reveals the rows wins.
- The reopened cell height is compared against the height captured at peak, never a fixed pixel
  value: the summary columns' renderers also influence row height, so an absolute would be brittle.

- target_layer rationale: multi-step UI flow with colour-coding dialogs, the gear settings
  panel, summary column menus, and a full layout+project persistence round-trip — all require
  a real browser; target_layer: playwright.
- Colour-coding reads use grid.cell(col, row).color (product read through the renderer).
  Never assert the raw df column tag directly, and never count global non-white pixels
  (near-vacuous for a grid). The per-cell product read is the authoritative signal.
- Row height: assert grid.cell().bounds.height — the cell bounds height is the observable
  signal. Do NOT assert props.rowHeight: setting a prop and reading it back is a vacuous
  round-trip that fails the false-green guard.
- Missing value color: assert grid.cell().color on a null cell. Do NOT assert
  props.missingValueColor for the same reason.
- Stats rows (min / max) have no product-readable signal (recon-verified: not in
  props/look/df tags). Their expectation bullet is explicitly an honest no-error floor plus
  summary-column-survival (the GROK-19809 guard). Do NOT invent a visibility claim.
- Summary columns: assert grid.columns.length increments and each new column's cellType
  matches the chosen type. The seven one-click types are: Sparklines, Bar Chart, Pie Chart,
  Radar, Smart Form, Tags, Confidence Interval. Form-designer entries in the same menu are
  manual-only and are excluded.
- Linked colour coding on WEIGHT requires Categorical on SEX to be active first (it is the
  source column precondition). Set SEX categorical before opening the Linked dialog.
- The Linked dialog's source-column selector is labelled "Source column" (or similar) in the
  Edit dialog — it is NOT a submenu entry; it is an input field inside the opened dialog.
- Project save and reopen use the REAL ribbon Save button, not a JS-API project.save() call.
  Use saveProjectViaUI(page, name) from helpers/projects.ts, which clicks the ribbon Save,
  fills the name, and dismisses the Share dialog. Reopen via openProjectFromDashboards(page,
  name). Do NOT reopen via grok.dapi.projects.find(id).open() — find() returns undefined
  here and the call throws before anything is asserted.
- Teardown in finally: delete the probe layout and the probe project via
  deleteProjectWithCleanup(page, {projectId}) so neither leaks across test runs.
- GROK-19769 guard: the layout re-apply step MUST include asserting summary column survival
  (grid.columns.length and each cellType). A layout re-apply that drops summary columns is
  the regression.
- GROK-19809 guard: adding stats rows must NOT reduce grid.columns.length from the value
  after Step 12 — the regression was the special-rows insertion dropping virtual columns.

---
{
  "order": 11,
  "datasets": ["System:DemoFiles/demog.csv"]
}
