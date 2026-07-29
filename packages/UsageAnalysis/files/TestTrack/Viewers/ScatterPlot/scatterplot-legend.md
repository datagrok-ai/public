---
feature: scatterplot
realizes_atlas:
  - scatterplot.cp.legend
realizes:
  - viewers.scatter-plot
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-17227
    status: fixed
  - id: GROK-17225
    status: fixed
  - id: GROK-14940
    status: fixed
  - id: GROK-17230
    status: fixed
  - id: GROK-19083
    status: fixed
  - id: GROK-20228
    status: fixed
realized_as:
  - scatterplot-legend-spec.ts
expected_results:
  - anchor: "Color legend and the joint Color plus Marker legend"
    expectation: "Setting Color to a categorical column renders one legend entry per
      category of that column; adding a Marker column adds one glyph entry per
      marker category into the same legend container; switching Color to a
      numeric column keeps the marker glyph entries present; switching Color
      back to a categorical column leaves exactly one entry per color category
      plus one per marker category, with no duplicated legend and no divergence
      when Color and Marker are set to the same column (GROK-17225 guard)"
  - anchor: "Clearing the Marker column removes its glyph entries"
    expectation: "With Color and Marker both set to the same two-category column,
      clearing the Marker column drops the marker glyph entries from the legend
      to zero while the color entries stay — the legend is rebuilt, not left
      stale (GROK-19083 guard)"
  - anchor: "Filtered-out categories are absent from the marker legend"
    expectation: "With the demog RACE filter narrowed to Asian and Caucasian (Black
      and Other deselected) and Markers set to RACE, the marker legend carries
      exactly two glyph entries, reading Asian and Caucasian — both when the
      Marker column is set BEFORE the filter is applied and when it is set after
      (GROK-17230 guard)"
  - anchor: "Filtered-out categories are absent from the marker legend — filter
      cleared"
    expectation: "Clearing the RACE filter restores the marker glyph entries for all
      four RACE categories — the legend follows the filter in both directions"
  - anchor: "Clicking a legend category hides the other categories on the plot"
    expectation: "A real mouse click on the Asian legend entry drops the
      settle-gated ink on the plot canvas sharply below the unselected baseline
      (the other categories stop being drawn), the legend marks the clicked
      entry as current, dims the other color entries and shows it is filtering,
      and the table's filtered row count is UNCHANGED — the legend filters the
      viewer, not the table. Applying the SEX = F Filter Panel filter on top
      leaves the plot's ink strictly below the legend-only reference by a
      calibrated margin — the rows the panel removed never come back (GROK-14940
      guard) — the M entry is absent from the legend while the panel keeps only
      F, and applying panel and legend in either order reaches the same plot,
      ink for ink"
  - anchor: "Clicking a legend category hides the other categories on the plot —
      second click restores the plot exactly"
    expectation: "A second real click on the same legend entry restores the plot's
      settle-gated ink EXACTLY to its pre-click value (no epsilon — the canvas
      does not drift at rest), and the legend is left with no entry marked as
      current, no entry dimmed and no filtering state"
  - anchor: "Clicking a legend category hides the other categories on the plot — the
      entry for empty values"
    expectation: "With a Color column carrying empty values, the legend offers an
      empty-values entry and a real click on it raises no console error, marks
      that entry as current and puts the legend into its filtering state; the
      saturated colour of the valued categories collapses on the plot while the
      pale grey the empty-value markers are drawn in grows — the empty category
      filters the viewer like any other (GROK-20228 guard) — the table's
      filtered row count is unchanged, applying the SEX = F panel filter on top
      leaves the ink strictly below the legend-only reference by a calibrated
      margin, and releasing the selection with the panel reset restores the ink
      exactly"
  - anchor: "A column with conditional color coding still produces a legend"
    expectation: "On the spgi-100 dataset, with Conditional color coding defined on
      the CAST Idea ID column over the ranges 634783-634820 and 634820-634885,
      and that column set as the scatter plot's Color, a legend is present and
      carries one entry per defined range, labelled with those two ranges
      (GROK-17227 guard)"
  - anchor: "Legend visibility hides and restores the legend"
    expectation: "Setting Legend Visibility to Never removes the legend container
      from the viewer (its count drops to zero); setting it to Always brings the
      legend back with the same entry composition recorded before hiding — the
      same color categories and the same marker glyph entries; restoring the
      original value leaves the legend present with that composition, and the
      whole ladder raises no page or console errors"
  - anchor: "Layout and project persistence at peak configuration — layout round-trip"
    expectation: "After saving the view layout with Color RACE and Markers SEX
      configured, adding a Histogram viewer and clearing Color, and re-applying
      the saved layout, the view's viewer set equals the SAVED set (a Scatter
      plot is present AND the later-added Histogram is absent), the restored
      Scatter plot's color and marker columns equal the saved values, and the
      legend re-renders with the color entries plus the marker glyph entries;
      the probe layout is deleted afterwards even on failure"
  - anchor: "Layout and project persistence at peak configuration — project save /
      Close All / reopen"
    expectation: "After saving the view as a project through the real ribbon Save
      button, Close All, and reopening the saved project, a Scatter plot viewer
      is present in the reopened view, its color and marker columns equal the
      persisted values, and its legend is present with the same entry
      composition — a cross-session round-trip; the probe project is deleted
      afterwards even on failure"
---

# Scatter Plot — Legend Lifecycle, Filter Interplay, Persistence

## Purpose

Verifies the scatter plot's legend end to end on the demog and spgi-100
datasets: the legend that a categorical Color column produces, the joint
legend when a Marker column is added, its consistency as Color and Marker are
switched independently, its rebuild when the Marker column is cleared, its
agreement with the Filter Panel in both orderings, the behaviour of clicking a
legend category (including the category for empty values), the legend a column
with conditional color coding produces, and the survival of the configured
legend across a saved layout and a saved project. The legend is real markup,
so entry counts and entry kinds are read directly. A legend click filters what
the VIEWER draws, not the table, so its outcome is read from the plot's own ink;
the table's filtered row count is checked only to confirm the effect stays local
to the viewer.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table
   view to load.
3. Add a Scatter plot to the current table view via the Toolbox viewer icon.
4. On the viewer, set the X column selector to WEIGHT and the Y column
   selector to HEIGHT.
5. Record the full row count of the table with no filtering applied.

## Scenarios

### Scenario 1: Color legend and the joint Color plus Marker legend

Steps:
1. On the viewer, click the **Color** selector label text — not only its
   triangle — and set Color to RACE.
2. Verify a legend appears with one entry per RACE category (four on demog),
   each carrying the category name.
3. Open the viewer settings (gear icon on the Scatter plot title bar) and, in
   the **Marker** section, set the **Markers** column to SEX.
4. Verify the same legend now also carries one glyph entry per SEX category
   (two on demog) — the color entries and the marker glyph entries live in one
   legend.
5. Switch Color to AGE, a numeric column.
6. Verify the marker glyph entries are still present — the marker legend
   survives a numeric Color.
7. Switch Color back to RACE.
8. Verify the legend carries exactly one entry per RACE category plus one
   glyph entry per SEX category — no doubled legend and no leftover entries
   from the numeric Color.
9. Set the **Markers** column to RACE as well, so Color and Marker are the
   same column.
10. Verify the color entries and the marker glyph entries agree on the same
    category set — the two halves of the legend do not diverge.
11. Set the **Markers** column back to SEX (the state the next scenario
    starts from).

Expected:
- A categorical Color column renders one legend entry per category
- Adding a Marker column adds one glyph entry per marker category to the same legend
- Switching Color to a numeric column and back leaves the marker glyph entries intact and produces no doubled legend
- With Color and Marker on the same column the color and marker entries cover the same category set

### Scenario 2: Clearing the Marker column removes its glyph entries

Steps:
1. Set Color to SEX and the **Markers** column to SEX — both on the same
   two-category column.
2. Verify the legend carries one entry per category, and that both entries
   carry a marker glyph — with Color and Markers on the SAME column the legend
   is joint, so the glyph sits on the colored entries rather than adding a
   separate pair.
3. Clear the **Markers** column (pick the empty first row in the column
   picker).
4. Verify the marker glyph entries are gone from the legend while the color
   entries remain — the legend was rebuilt, not left stale.
5. Set Color back to RACE (round-trip revert; Markers stays cleared for the
   next scenario).

Expected:
- Clearing the Marker column drops the marker glyph entries to zero while the color entries stay

### Scenario 3: Filtered-out categories are absent from the marker legend

Steps:
1. Open the Filter Panel from the Toolbox **Filters** section and wait for its
   filter cards to finish building.
2. In the RACE filter, keep only **Asian** and **Caucasian** selected —
   deselect **Black** and **Other**.
3. Set the **Markers** column to RACE — that is, set Marker AFTER the filter.
4. Verify the marker glyph entries are exactly two and read Asian and
   Caucasian; the deselected Black and Other do not appear.
5. Clear the RACE filter, then clear the **Markers** column.
6. Reverse the ordering: set the **Markers** column to RACE FIRST, then again
   keep only Asian and Caucasian selected in the RACE filter.
7. Verify the marker glyph entries are again exactly Asian and Caucasian — the
   ordering does not matter.
8. Clear the RACE filter and verify the marker glyph entries come back for all
   four RACE categories (round-trip).

Expected:
- With the RACE filter narrowed to Asian and Caucasian, the marker legend lists exactly those two categories, both when the Marker column is set after the filter and when it is set before it
- Clearing the filter restores the marker glyph entries for all four RACE categories

### Scenario 4: Clicking a legend category hides the other categories on the plot

A legend click filters what the VIEWER draws, not what the table holds: the
points of the other categories stop being painted, while the table's own
filtered row count is untouched. The graded signal is therefore the amount of
ink on the plot canvas, measured with the pointer parked away from the plot
(a hovered point paints a mouse-over marker and spoils the frame) and read
only once the canvas has settled. That the table's filtered row count does not
move is asserted too — as the confirmation that the legend's effect is local
to the viewer, not as the main signal.

The guard is that a legend click must respect what the Filter Panel already
removed. It is graded by comparing against the right reference: the plot with
the legend category selected and NO panel filter. Comparing against the
panel-only plot would be wrong, because removing the other categories also
removes marker overlap, so the remaining category's own ink GROWS.

Steps:
1. With Color RACE and the **Markers** column SEX, no Filter Panel filter and
   no legend selection, park the pointer away from the plot and record the
   plot's ink once it settles; record the table's filtered row count.
2. Click the **Asian** entry in the legend with a real mouse click.
3. Verify the plot's ink drops sharply against step 1 — the points of the
   other categories are no longer drawn.
4. Verify the legend shows it is filtering: the clicked entry is marked as the
   current one and the other color entries are dimmed.
5. Verify the table's filtered row count is unchanged from step 1 — the legend
   filters the viewer, not the table.
6. Record this state's ink as the reference: legend selection, no panel
   filter.
7. Open the Filter Panel and, in the SEX filter, deselect **M**, leaving only
   **F**; verify the table's filtered row count drops below the full row
   count.
8. Verify the plot's ink is now strictly below the reference recorded in
   step 6 — the rows the Filter Panel removed did not come back when the
   legend narrowed the view further.
9. Verify the legend no longer offers the **M** entry while the panel keeps
   only F.
10. Reverse the order: release the legend selection, reset the Filter Panel,
    then apply the SEX filter FIRST and click **Asian** in the legend second.
11. Verify this state's ink equals the measurement from step 8 exactly — the
    two orders reach the same plot.
12. Release the legend selection and reset the Filter Panel (round-trip
    revert).

Expected:
- Clicking a legend category drops the plot's ink sharply — the other categories stop being drawn
- The legend marks the clicked entry as current, dims the other color entries, and shows it is filtering
- The table's filtered row count is unchanged by the legend click — the effect is local to the viewer
- With a Filter Panel filter applied on top, the plot's ink is strictly below the legend-only reference — the filtered-out rows never come back
- The M entry is absent from the legend while the panel keeps only F
- Applying the panel filter and the legend selection in either order reaches the same plot, ink for ink

#### Second click restores the plot exactly

Steps:
1. From the clean state (no panel filter, no legend selection), record the
   plot's ink.
2. Click the **Asian** legend entry, then click it again to release it.
3. Verify the plot's ink is exactly the value recorded in step 1 — the release
   restores the frame, not merely something close to it.
4. Verify the legend no longer marks any entry as current and no entry is
   dimmed.

Expected:
- Releasing the legend selection restores the plot's ink exactly to its pre-click value
- No legend entry is left marked as current and none is left dimmed

#### The entry for empty values

Steps:
1. Add a Color column that contains empty values (the demog categorical
   columns have none, so the fixture is a categorical column carrying missing
   values, added to the table) and set it as the Color column.
2. Verify the legend carries an entry for the empty values.
3. Park the pointer and record the plot's ink, separating the saturated colour
   of the valued categories from the pale grey the empty-value markers are
   drawn in.
4. Click the empty-values entry in the legend with a real mouse click.
5. Verify the click raised no console error, the legend marks that entry as
   current and shows it is filtering, and the saturated color on the plot
   collapses while the pale grey grows — only the empty-value points are left
   drawn.
6. Verify the table's filtered row count is unchanged — again the effect is
   local to the viewer.
7. Record this state's ink as the reference, then apply the SEX filter (only
   **F**) from the Filter Panel and verify the plot's ink is strictly below
   that reference.
8. Release the legend selection and reset the Filter Panel; verify the plot's
   ink returns exactly to the value recorded in step 3.
9. Remove the fixture column and set Color back to RACE and the **Markers**
   column to SEX (round-trip revert; the teardown runs even when a
   verification fails).

Expected:
- The legend carries an entry for empty values and clicking it raises no console error
- Clicking the empty-values entry collapses the saturated colour on the plot while the pale grey of the empty-value markers grows — the empty category behaves like any other
- The table's filtered row count is unchanged by the click
- With the panel filter applied on top, the plot's ink is strictly below the legend-only reference
- Releasing the selection and resetting the panel restores the plot's ink exactly

### Scenario 5: A column with conditional color coding still produces a legend

Steps:
1. Open the spgi-100 dataset: `System:AppData/Chem/tests/spgi-100.csv` — wait
   for the table view to load. Do not rely on the table's name; the loader
   does not name a table after its file.
2. On the grid, right-click the **CAST Idea ID** column header and set
   **Color Coding > Conditional**, defining two ranges that both occur in the
   data: **634783-634820** and **634820-634885**.
3. Add a Scatter plot to this table view via the Toolbox viewer icon and set
   the X column selector to **CAST Idea ID** and the Y column selector to
   **Idea ID** — the two columns the viewer auto-picks on this dataset.
4. Set the **Color** selector to **CAST Idea ID** — the column that now
   carries the conditional color coding.
5. Verify a legend is present and carries one entry per defined range,
   labelled 634783-634820 and 634820-634885 — a column coloured by rules
   produces a legend just like a plain categorical column does.
6. Close this table view and return to the demog view (round-trip revert).

Expected:
- With conditional color coding defined on CAST Idea ID and that column set as Color, a legend is present with one entry per defined range, labelled 634783-634820 and 634820-634885

### Scenario 6: Legend visibility hides and restores the legend

Legend Visibility is the one legend setting whose outcome is the presence of
the legend itself: the Never value removes the legend from the view entirely,
and Always or Auto bring it back. The ladder therefore runs on presence, and
the restored legend is compared against the composition recorded before it was
hidden.

Steps:
1. With Color RACE and **Markers** SEX still configured, read the current
   **Legend Visibility** value in the **Legend** section of the viewer settings
   (a freshly added viewer carries Auto) and record the legend's entry
   composition.
2. Set **Legend Visibility** to `Never`.
3. Verify the legend is gone from the viewer — no legend container remains.
4. Set **Legend Visibility** to `Always`.
5. Verify the legend is back and carries the same entries as recorded in
   step 1 — the same color categories and the same marker glyph entries.
6. Set **Legend Visibility** back to the value recorded in step 1 (round-trip
   revert) and verify the legend is still present with that composition.

Expected:
- Setting Legend Visibility to Never removes the legend from the viewer
- Setting it to Always brings the legend back with the same entry composition it had before
- Restoring the original value leaves the legend present and unchanged

### Scenario 7: Layout and project persistence at peak configuration

The persistence tail runs at the configuration the scenarios above reach —
Color RACE, Markers SEX, with the legend rendered — and nothing is reverted
before saving.

Steps:
1. Confirm the peak state on the demog view: Color RACE, **Markers** SEX, the
   legend showing the color entries plus the marker glyph entries. Record the
   color column, the marker column, and the legend entry composition.
2. Save the view layout.
3. Perturb the view: add a **Histogram** viewer to it through the Toolbox
   viewer icon and clear the **Color** column.
4. Re-apply the saved layout.
5. Verify the view's viewer set equals the SAVED set: a Scatter plot is
   present AND the later-added Histogram is absent.
6. Verify the restored Scatter plot carries the recorded color and marker
   columns, and that its legend re-rendered with the recorded entry
   composition.
7. Delete the probe layout.
8. Save the view as a project through the ribbon **Save** button, dismiss the
   Share dialog that follows, Close All, then reopen the saved project.
9. Verify a Scatter plot viewer is present in the reopened view, its color and
   marker columns equal the recorded values, and its legend is present with
   the same entry composition — a cross-session round-trip.
10. Delete the probe project.

Expected:
- Re-applying the saved layout restores the SAVED viewer set (Scatter plot present, the later-added Histogram absent) with the configured color and marker columns and the legend re-rendered
- Reopening the saved project restores a Scatter plot with the same color and marker columns and the same legend composition
- The probe layout and project are deleted even when a verification fails — they never leak

## Automation notes

- Narrowing a categorical filter (Scenario 3 step 2, Scenario 4 step 7) is
  driven through the Filter Panel's filter-group API rather than by clicking the
  card's canvas. The guard needs a DETERMINISTIC surviving category set, and the
  card's checkbox list is canvas-drawn: a coordinate click can toggle *a*
  category but cannot choose *which* one. Where a guard only needs "exactly one
  category left, whichever it is", the real coordinate click is used instead —
  see the labels-tooltip scenario. The narrowing is setup for the graded signal,
  never the signal itself.
- The viewer handle is
  `grok.shell.tv.viewers.find(v => v.type === 'Scatter plot')`; the viewer is
  added via the Toolbox icon `[name="icon-scatter-plot"]`. Resolve the viewer
  root as the `[name="viewer-Scatter-plot"]` element that is NOT inside a
  `.d4-dialog` — the Formula Lines dialog embeds a second viewer with the same
  name.
- Legend reads, all inside the viewer root: the container is
  `[name="legend"]`; `.d4-legend-item` counts all entries,
  `.d4-legend-item-coloring` the color categories (each with a
  `span.d4-legend-value` carrying the category text and a `span.d4-legend-cross`),
  `.d4-legend-item-extra` the marker glyph entries when Color and Markers are on
  DIFFERENT columns. CAVEAT: when they are on the SAME column the legend renders
  JOINTLY — `-extra` is 0 and the glyph sits on the coloring entries themselves,
  so a `-extra` count is vacuous there. Define marker entries uniformly as
  entries CARRYING a glyph (`i[name="legend-icon-color-picker"]`), which is
  correct in both configurations. Assert counts and the category-name set, never
  a screenshot. `legendVisibility: Never` removes the container from
  the markup, so presence is a valid assert for that property only.
- Scenario 6, Legend Visibility: the property is `legendVisibility`
  (`Auto | Always | Never`), settings-panel row
  `[name="prop-legend-visibility"]` with view cell
  `[name="prop-view-legend-visibility"]` in the collapsed `legend` category;
  clicking the view cell reveals a `select.property-grid-item-editor-spinner`.
  Enumerate the select's options and pick the value from them rather than
  assuming the wording; if the panel row does not yield a select, fall back to
  the context-menu mirror `div-Properties...---Legend---Legend-Visibility---
  <value>` and, again, enumerate the leaf before clicking it. Whichever route
  is used, read `legendVisibility` back afterwards — a route that silently does
  nothing must fail rather than let the presence assert pass vacuously. The
  presence signal is the count of `[name="legend"]` inside the viewer root
  (1 → 0 → 1), and the restored composition is compared against the counts and
  labels recorded before hiding.
- LEGEND CLICKS REQUIRE TRUSTED INPUT. Scripted clicks — `.click()` and a full
  `mousedown`+`mouseup`+`click` chain, on the value span and on the cross
  alike — do not actuate the legend at all. Use a real coordinate click
  (`page.mouse.click` over the entry's bounding box). A scripted click has no
  effect at all; a real click puts the legend into its filtering state and
  hides the other categories on the plot (see the Scenario 4 notes below).
- SCENARIO 4 — THE GRADED SIGNAL IS INK ON THE PLOT, NOT THE TABLE FILTER. A
  legend click filters what the viewer DRAWS; the table's filter bitset is not
  touched, and the viewer exposes no readable state for which categories are
  selected (nothing in its property list names them, `rowSource` stays
  `Filtered` and `filter` stays empty — the selection lives on the Dart side).
  The click also emits nothing but render events, which fire on any repaint and
  prove nothing. So the assert is a settle-gated per-colour ink measurement of
  `canvas[name="canvas"]` inside the viewer root that is not in a `.d4-dialog`,
  and `df.filter.trueCount` is asserted UNCHANGED as the locality confirmation.
- Ink measurement discipline: PARK THE POINTER away from the plot before every
  measurement — a hovered point paints its mouse-over indicator into the frame.
  Read the canvas twice and require the two reads to agree before using the
  value; the canvas does not drift at rest, so the settle gate can be strict
  and the restore-after-release assert can be exact equality with no epsilon.
- REFERENCE STATE — compare "panel filter + legend selection" against "legend
  selection with NO panel filter", never against "panel filter only". Removing
  the other categories also removes marker overlap, so the SELECTED category's
  own ink GROWS; a monotonic-ink assert on the remaining category would be
  wrong in both directions. The right statement is that the combined state
  carries strictly less ink than the legend-only reference — the rows the panel
  removed did not reappear. Applying the two in the reverse order reaches the
  same frame, which is asserted as equality.
- Colour buckets are NOT one-to-one with categories (one category lands in the
  neighbouring hue bucket, another is smeared across two because of overlap
  shading). Assert on total ink, on the saturated-versus-pale-grey split, or on
  measured state-to-state comparisons — never on an assumed
  "category → colour" map.
- Scenario 4, empty-value category: the fixture is a categorical column
  carrying missing values, added through the JS API and removed in a `finally`
  teardown. Empty-value markers are painted in a pale grey rather than a
  category colour, so they land in the pale bucket: classify pixels by
  saturation (a near-neutral pixel with a high channel value is an empty-value
  marker; the dark neutrals are axes and text) and expect the saturated ink to
  collapse while the pale ink grows. Legend markup is the same as for any other
  category, so the filtering-state and current-entry asserts carry over
  unchanged.
- Legend DOM is the secondary support for all of the above: `[name="legend"]`
  gains `d4-legend-filtering`, the clicked entry gains `d4-legend-item-current`,
  and the other colour entries drop to `opacity: 0.5`. Use it to confirm the
  click actuated; use the ink to grade what it did.
- Thresholds come from live calibration at spec-authoring time and are kept as
  margins, never as settled absolute values: a sharp-drop assert as a fraction
  of the measured baseline, and a fixed pixel margin for the
  strictly-less comparisons that is a fraction of the measured gap.
- Scenario 5, conditional color coding: the validated path is the column
  metadata API — `col.meta.colors.setConditional({'634783-634820': …,
  '634820-634885': …})` on `CAST Idea ID`, which sets the color-coding-type and
  conditional tags. Those two ranges are the live-verified ones on spgi-100 and
  they produce exactly two legend entries; ranges that match no data collapse
  into a single `no value` entry, so do not invent other bounds. The
  grid-column context-menu route to Conditional color coding was not exercised
  in recon and is assumed equivalent — drive it if it works, otherwise use the
  metadata API and say so in the spec. Mind the column-name casing: the live
  column is `CAST Idea ID` (the bug text spells it `Cast Idea Id`), and the Y
  column `Idea ID` is the live auto-pick partner.
- Scenario 7, the perturbing viewer is a **Histogram** — deliberately not
  another Scatter plot, so the "viewer set equals the SAVED set" check can tell
  the restored viewer from the perturbation by type
  (`[name="icon-histogram"]` in the Toolbox; the viewer set is read as
  `grok.shell.tv.viewers` types).
- The Filter Panel: the Toolbox `[name="div-section--Filters"]` entry can be
  present but invisible — guard on `offsetParent` and fall back to
  `tv.getFiltersGroup()`. The panel builds slowly (on spgi-100 it needs a long
  poll), so poll for its `.d4-filter` card count rather than waiting a fixed
  time. Reset is `[name="viewer-Filters"] [name="icon-arrow-rotate-left"]`,
  with no confirmation. Never query `[name="viewer-Grid"]` globally — the
  filter panel and the column popup both contain elements with that name.
- Column drives: the on-viewer selectors are
  `[name="div-column-combobox-x"]`, `-y`, `-color`, `-size` (lowercase role
  names, scoped to the viewer root — the property panel reuses the same
  names). They open on a synthetic `mousedown` on
  `.d4-column-selector-column` (a synthetic `.click()` does NOT open them),
  then real typing plus `Enter` commits; the popup grid is canvas-rendered, so
  the offered column names are not readable as text. The Markers column row is
  `[name="prop-markers"]` in the settings panel, carrying a nested
  `[name="div-column-combobox-markers"]`; the alternative route is the
  `Markers ▸ Marker Column` context-menu leaf. Clearing a column sets the
  property to an empty string, not to null — compare accordingly.
- Scenario 7: the layout is saved and re-applied via the JS API
  (`tv.saveLayout()` / `grok.dapi.layouts.save` / `tv.loadLayout`) — the
  View > Layout menu path has no headless handles and the round-trip end state
  is the same; the project is saved through the real ribbon Save button
  (`[name="button-Save"]` → `[name="dialog-Save-project"]` → `input#name` →
  OK), because only the UI Save captures the view layout, and the Share dialog
  that follows is dismissed via its CANCEL button. Probe names carry a
  timestamp suffix so concurrent runs never collide, and both the probe layout
  and the probe project are deleted in `finally` teardowns
  (`grok.dapi.layouts.delete` / `grok.dapi.projects.delete`).
- The project-publish preview clones the live view into an offscreen frame and
  emits a benign "Unable to find element in cloned iframe" plus a Dart null
  error pair against the detached view. Whitelist that class ONLY inside the
  ribbon project-save window; everywhere else — notably the legend-click
  console guards — the same error class is the regression signal and must
  reach the guard. Console guards must also filter the shared dev server's
  unrelated WebSocket reconnect noise and the WebGPU `powerPreference`
  warning.
- The scatter plot title bar has no close icon: reach it with
  `root.closest('.panel-base')` (not a fixed parent-hop count) and close the
  viewer from the hamburger menu or via `viewer.close()`.

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
