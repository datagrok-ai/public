---
feature: formsviewer
realizes_atlas:
  - formsviewer.cp.color-and-renderer-presentation
  - formsviewer.edge.molecule-column-render
  - formsviewer.edge.fit-semtype-promotes-renderer-size
realizes:
  - viewers.form
priority: p2
target_layer: playwright
realized_as:
  - forms-spec.ts
coverage_type: regression
related_bugs:
  - id: GROK-19339
    status: fixed
  - id: GROK-14963
    status: fixed
  - id: GROK-19473
    status: fixed
  - id: GROK-18814
    status: fixed
expected_results:
  - anchor: "Step 1a"
    expectation: >-
      With the AGE grid column colour-coded BY BACKGROUND, the computed
      background colour of the AGE field inside the current-row card equals the
      colour the AGE colour-coding scheme assigns to that same row. Both
      readings are taken at runtime — the field colour from the rendered
      element's computed style, the reference colour from the column's colour
      scheme — and are normalized to the same notation before comparison (the
      scheme yields "#rrggbb", the computed style yields "rgb(r, g, b)"). No
      expectation of any kind is placed on the field's TEXT colour.
  - anchor: "Step 1b"
    expectation: >-
      The AGE field's horizontal alignment and font equal the alignment and font
      the grid uses for the AGE column, both read at runtime from the grid
      column's own cell style rather than from constants. The alignment is set
      to a non-default value (Center) before the check, so a field that ignored
      the grid style could not pass by inheriting the same default.
  - anchor: "Step 1c"
    expectation: >-
      After turning Color Code OFF, the AGE field's computed background equals
      the background of a field whose column carries no colour coding, read from
      the same card in the same pass (the uncoloured default computed on the
      fly, never a hardcoded colour). After turning Color Code back ON, the
      background returns to the colour-scheme colour of Step 1a.
  - anchor: "Step 2a"
    expectation: >-
      With Renderer Size at its default small, the width and height ATTRIBUTES
      of the Structure field's canvas equal the renderer's declared default
      width and height scaled by 0.66 and rounded down, then multiplied by the
      display's pixel density. Where the renderer declares no default size of
      its own, the viewer's own small-step fallback size is used as the base.
      Both the base size and the pixel density are read at runtime; no pixel
      number is hardcoded and the drawn picture is never measured.
  - anchor: "Step 2b"
    expectation: >-
      After switching Renderer Size to normal, the same canvas's width and
      height attributes equal the same runtime-computed base size at factor 1,
      multiplied by the display's pixel density — and are strictly larger than
      the values recorded at Step 2a.
  - anchor: "Step 2c"
    expectation: >-
      After switching Renderer Size to large, the same canvas's width and height
      attributes equal the runtime-computed base size scaled by 1.5 and rounded
      down, multiplied by the display's pixel density — and are strictly larger
      than the values recorded at Step 2b.
  - anchor: "Step 3a"
    expectation: >-
      On spgi-100, the Structure field inside a card is a CANVAS element
      carrying the Structure column's name in its column attribute — not an
      INPUT holding a SMILES or molblock string. No error balloon appears and no
      console error is raised while the molecules paint.
  - anchor: "Step 3b"
    expectation: >-
      With the field set picked as Structure, Core and Primary Series Name, each
      card carries TWO separate molecular canvases — one per molecule column,
      each identified by its own column name — while Primary Series Name remains
      an ordinary text field. The count of molecular canvases per card equals
      the count of picked molecule columns, computed from the picked field set.
  - anchor: "Step 4a"
    expectation: >-
      On curves.csv, a Forms viewer added to a table that already carries a
      curve (semantic type "fit") column renders that field's canvas at the
      NORMAL step of the size ladder without the user touching Renderer Size:
      the canvas width and height attributes equal the runtime-computed base
      size at factor 1 times the pixel density, and differ from the small-step
      values computed the same way. The promotion is read through the canvas
      size, never through a read-back of the Renderer Size property.
  - anchor: "Step 4b"
    expectation: >-
      On curves.csv both the molecule column and the curve column render as
      CANVAS elements carrying their own column names — the curve column is not
      an INPUT holding raw JSON text. No error balloon and no console error is
      raised.
  - anchor: "Step 4c"
    expectation: >-
      With three rows selected, every drawn card carries its own molecule
      canvas and its own curve canvas; the number of cards holding a full set of
      canvases equals the drawn-card count (the current-row card plus one per
      selected row), computed at runtime from the selection.
  - anchor: "Step 5a"
    expectation: >-
      With a substructure filter driving the filter on spgi-100, the selected-row
      cards (the positions after the leading offset) still equal the
      runtime-recomputed intersection of the selection with the filter: a selected
      row that the substructure filter excludes produces no selected-row card, and
      the surviving selected rows still do. The intersection is recomputed inside
      the poll on every iteration, together with the card read, and the check
      requires that the filter actually excluded at least one selected row.
  - anchor: "Step 5b"
    expectation: >-
      The molecular field of a card that SURVIVES the substructure filter is
      rebuilt and repainted after the filter is applied. Proven structurally,
      not by a picture diff: the pre-filter canvas node (tagged before the
      filter) is gone from the DOM after the settled re-render, the same
      surviving row (identified by its id-field value) carries a FRESH Structure
      canvas, and that canvas has real drawn content — opaque, non-white pixels
      counted at runtime, never compared to a reference image — with a
      zero-console-error floor around the re-render. A same-row canvas-content
      diff is deliberately NOT used: a plain filter carries no substructure
      highlight, so the redrawn molecule is byte-identical for the same row (the
      length change the refdoc recorded was a fixed-position read, i.e. a
      different molecule once the card set shrank). Makes NO claim that the form
      picture matches the grid picture (out of scope — see forms-ui.md).
---

# Forms Viewer — Colour Coding and Renderer Presentation (p2)

## Setup

1. Close all open views and dialogs.
2. Open **System:DemoFiles/demog.csv** so the table view is active.
3. Add the **Forms** viewer from **Toolbox > Viewers** (click the Forms icon in
   the viewer toolbox on the right side of the screen).

Scenarios 2, 3 and 5 run on **System:AppData/Chem/tests/spgi-100.csv** (the
table is named `spgi-100`) and Scenario 4 on **System:DemoFiles/curves.csv**;
each of them opens its own dataset and adds its own Forms viewer, as stated in
its first step. The setup above is the entry state only — mounting the viewer,
the default field set and the persistence of the layout are covered by
`formsviewer-forms-core.md` and are not re-checked here. The card set =
selection ∩ filter invariant is likewise owned by that core scenario; Scenario 5
re-exercises the same invariant here only because the driver differs — a
MOLECULAR substructure filter rather than a plain value filter — and pairs it
with the repaint-by-node-identity check that is this scenario's own subject.

## Scenarios

### Scenario 1: Colour coding — the field picks up the grid column's background, alignment and font (GROK-19339)

Steps:

1. In the grid, right-click the **AGE** column header and apply **Color Coding >
   Linear**. Confirm the colouring lands on the cell **background**: the AGE
   cells show coloured backgrounds while their text keeps the default colour.
   This is a hard precondition for the whole scenario — a column coloured by
   TEXT instead produces no field background at all and there would be nothing
   to read.
2. Open the AGE grid column's settings and set two non-default values: its
   horizontal alignment to **Center** and its font to a distinct font (a bold
   italic serif, unlike the grid default). Both are deliberately non-default, so
   that a field which ignored the grid's style could not pass by inheriting the
   same default — the font half is what makes Step 1b's font comparison
   meaningful rather than a default-vs-default match.
3. In the grid, click row 0 to make it the current row, so the current-row card
   is populated. A card built for no row renders no field elements at all, so
   this must happen before any field is read.
4. Read the background colour rendered on the **AGE** field inside the card that
   carries the green current-row indicator. Compare it with the colour the AGE
   colour-coding scheme assigns to that same row, bringing both readings to the
   same colour notation first. Label this **Step 1a**.
5. Read the horizontal alignment and the font of that same AGE field and compare
   them with the alignment and font the grid uses for the AGE column. Label this
   **Step 1b**.
6. Do **not** read and do **not** assert the field's TEXT colour anywhere in
   this scenario — see the hard constraint in Automation notes.
7. In the Forms viewer's property panel (**Context Panel > Misc**), turn
   **Color Code** OFF. Wait for the re-render to settle. Read the AGE field's
   background again and compare it with the background of a field whose column
   carries no colour coding at all (read from the same card in the same pass) —
   the coloured background is gone and the field is back to the uncoloured
   default.
8. Turn **Color Code** back ON, wait for the re-render, and confirm the AGE
   field's background returns to the colour-scheme colour from Step 1a — the
   round-trip closes. Label steps 7–8 together as **Step 1c**.

Expected:

- Step 1a: the AGE field's computed background equals the colour the AGE
  colour-coding scheme assigns to that row, after both readings are normalized
  to the same notation. Nothing is asserted about the text colour.
- Step 1b: the field's horizontal alignment and font equal the grid column's,
  read at runtime from the grid column's cell style.
- Step 1c: with Color Code off the AGE field's background equals the uncoloured
  default read from another field in the same card; turning it back on restores
  the colour-scheme colour.

### Scenario 2: Renderer size ladder — small, normal, large (GROK-14963)

Steps:

1. Close all open views. Open **System:AppData/Chem/tests/spgi-100.csv** and add
   the **Forms** viewer from **Toolbox > Viewers**.
2. In the grid, click row 0 to make it the current row so the current-row card is
   populated.
3. In the property panel, read **Renderer Size** and confirm it is still at its
   default **small** — read it before writing anything, so the default is
   established rather than assumed.
4. Read the **width** and **height** attributes of the drawing surface rendered
   for the **Structure** field inside the current-row card. Compare them against
   the renderer's declared default size for that column, scaled by 0.66 rounded
   down and multiplied by the display's pixel density. Where the renderer
   declares no default size of its own, the viewer's own small-step fallback
   size is the base instead. Record the values. Label this **Step 2a**.
5. Set **Renderer Size** to **normal** in the property panel and wait for the
   re-render. Read the same two attributes and compare them against the same
   runtime-computed base size at factor 1, times the pixel density; confirm they
   are strictly larger than the values recorded at Step 2a. Label this
   **Step 2b**.
6. Set **Renderer Size** to **large** and wait for the re-render. Read the same
   two attributes and compare them against the base size scaled by 1.5 rounded
   down, times the pixel density; confirm they are strictly larger than the
   values recorded at Step 2b. Label this **Step 2c**.
7. Never measure the drawn picture, and never write a pixel number into the
   expectation — every number in steps 4–6 is computed at runtime from the
   renderer's declared defaults and the display's pixel density.

Expected:

- Step 2a: at small, the canvas width/height attributes equal the runtime-derived
  base size × 0.66 rounded down × pixel density.
- Step 2b: at normal, the same attributes equal the base size × 1 × pixel
  density, and exceed the small values.
- Step 2c: at large, the same attributes equal the base size × 1.5 rounded down
  × pixel density, and exceed the normal values.

### Scenario 3: Molecular fields render as drawing surfaces, several at once

Steps:

1. Continue in the spgi-100 view from Scenario 2 (or reopen it and add Forms;
   set Renderer Size back to **normal** so the surfaces are comfortably sized).
2. Confirm that the **Structure** field inside a card is a drawing surface
   carrying the Structure column's name, and **not** a text field holding a
   SMILES or molblock string. Confirm no error balloon is shown and no console
   error was raised while the molecules painted. Label this **Step 3a**.
3. Open the **Fields** picker from the property panel and set the field set to
   **Structure**, **Core**, **Primary Series Name** — in that order. `Core` sits
   beyond the viewer's 20-field cap on this 88-column table, so it only enters
   the form when it is picked explicitly. Then Ctrl+click two or three rows in the
   grid so several cards are shown — "each card carries two molecular canvases" is
   only a plural claim when more than the single current-row card is drawn.
4. Confirm each card now carries **two** separate molecular drawing surfaces,
   one per molecule column, each identified by its own column name, while
   **Primary Series Name** remains an ordinary text field. Compare the number of
   molecular surfaces per card against the number of molecule columns in the
   picked field set. Label this **Step 3b**.

Expected:

- Step 3a: the Structure field is a canvas carrying the column name, not a text
  input; no error balloon and no console error.
- Step 3b: two molecular canvases per card, one per picked molecule column, each
  carrying its own column name; the non-molecule field stays a text field.

### Scenario 4: Curve columns promote the renderer size on attach (GROK-19473)

Steps:

1. Close all open views. Open **System:DemoFiles/curves.csv**. Before adding the
   viewer, confirm the table carries a curve column — a column whose semantic
   type is `fit`, for example **multiple prefit**.
2. Add the **Forms** viewer from **Toolbox > Viewers**. Do **not** touch
   **Renderer Size**: the promotion this step checks happens while the viewer
   attaches, and only from the default small.
3. In the grid, click row 0 to make it the current row.
4. Open the **Fields** picker and set the field set to **smiles** and
   **multiple prefit**.
5. Read the **width** and **height** attributes of the curve field's drawing
   surface and compare them against the **normal** step of the size ladder,
   computed at runtime exactly as in Scenario 2 (base size × 1 × pixel density);
   confirm they differ from the small step computed the same way. The promotion
   is established through the drawn surface's size, never by reading the
   Renderer Size property back. Label this **Step 4a**.
6. Confirm that **smiles** renders as a drawing surface and that **multiple
   prefit** also renders as a drawing surface rather than as a text field
   holding raw JSON. Confirm no error balloon and no console error. Label this
   **Step 4b**.
7. In the grid, Ctrl+click three rows to select them. Wait for the re-render to
   settle. Confirm every drawn card carries its own molecule surface and its
   own curve surface, comparing the number of fully-drawn cards against the
   drawn-card count (the current-row card plus one per selected row) computed
   from the selection. Label this **Step 4c**.

Expected:

- Step 4a: the curve field's canvas is sized at the normal step of the ladder on
  a freshly attached viewer whose Renderer Size was never touched, and differs
  from the small step.
- Step 4b: both the molecule column and the curve column render as canvases
  carrying their column names; no raw JSON text field, no error balloon, no
  console error.
- Step 4c: each drawn card carries its own molecule and curve surfaces.

### Scenario 5: Substructure filter drives the card set and repaints the molecular field (GROK-18814)

Steps:

1. Close all open views. Open **System:AppData/Chem/tests/spgi-100.csv** and add
   the **Forms** viewer. Confirm **Show Selected Rows** is on in the property
   panel by reading it before writing anything.
2. In the grid, Ctrl+click four rows to select them, spread over the table so
   that a substructure query can plausibly exclude some of them. Wait for the
   re-render to settle and read a text field from each selected-row card so the
   starting card set is known.
3. Pick one selected row that will survive the filter and TAG its **Structure**
   surface's canvas node with a temporary marker attribute **before** the filter
   is applied (so Step 5b can prove the node was rebuilt).
4. Apply a substructure filter on the **Structure** column so that part of the
   dataset is filtered out. The Chem sketcher cannot be driven headless, so the
   filter is applied programmatically — a substructure query runs against the
   Structure column and its match bitset is pushed to the table filter, the same
   channel the filter panel's sketcher would ultimately drive.
5. Wait for the filtered re-render to settle, polling rather than sleeping.
   Confirm the selected-row cards (the positions after the leading offset) equal
   the intersection of the selection with the filter, recomputed on every poll
   iteration together with the card read: a selected row that the substructure
   filter excludes produces no selected-row card, and the surviving selected rows
   still do. Confirm the filter actually removed at least one selected row,
   otherwise the check proves nothing. Label this **Step 5a**.
6. Wait for the settled re-render and confirm the surviving row's card (found by
   its id-field value, since the card set shifts) now carries a FRESH **Structure**
   canvas — a different DOM node, the tagged one gone — that is actually drawn
   (non-empty, computed at runtime), with no console error raised. Label this
   **Step 5b**.
7. Do **not** compare the form's picture with the grid's picture, and do **not**
   use a non-white pixel count as a similarity or repaint measure. (A runtime
   draw-floor on the single field canvas — "> 20 drawn" opaque non-white pixels,
   well below the ~660 a real molecule paints, so it rejects a blank canvas
   without false-redding a drawn one — is a floor, not a similarity measure, and
   is allowed.) Whether the alignment and substructure highlighting visually agree
   is a manual check in `forms-ui.md`.

Expected:

- Step 5a: the selected-row cards equal the runtime-recomputed intersection of
  the selection with the filter, with the molecular filter as the driver.
- Step 5b: the molecular field of a surviving card is rebuilt (its canvas node is
  replaced) and drawn (fresh canvas non-empty, zero console errors); nothing is
  claimed about the two pictures agreeing.

## Automation notes

- **The TEXT colour is never assertable — hard constraint.** The source assigns
  the field's text colour the string `'<colour>!important;'`, which is not a
  valid CSS declaration value, so the browser discards it and the text keeps its
  inherited colour (forms-viewer.ts#L422-L426). Only the **background** is
  observable, and only when the grid column is colour-coded by BACKGROUND — the
  background is written in the `else` branch of `if (gc?.isTextColorCoded)`, so
  a text-colour-coded column produces no background at all. This is a documented
  source limitation, not a test defect, and it must not be "fixed" by asserting
  a text colour later.
- **Colour comparison needs format normalization**: the colour scheme yields
  `#rrggbb` (`DG.Color.toHtml`) while `getComputedStyle` yields `rgb(r, g, b)`.
  Parse both to the same representation before comparing. The colour
  `4294967295` (0xFFFFFFFF) means "no colour" and is excluded by the source.
- **The uncoloured default is read at runtime**, from a field whose column has no
  colour coding in the same card and the same pass — never a hardcoded colour
  string, which would fail the day the theme changes.
- **Renderer size is asserted on the canvas ATTRIBUTES, never on the picture.**
  `getRendererSize()` scales the renderer's declared `defaultWidth`/
  `defaultHeight` by `Math.floor(x * 0.66)` / `x` / `Math.floor(x * 1.5)`
  (forms-viewer.ts#L74-L87), and the canvas attributes are that size multiplied
  by `window.devicePixelRatio` (#L382-L386). The floor is in the source and an
  exact comparison without it breaks on odd default widths.
- **The spgi-100 molecular renderer declares no default size**, so the fallback
  `getSize()` applies: 120×60 / 200×100 / 300×150 CSS px for small / normal /
  large (#L89-L96). Measured live at `devicePixelRatio = 1.25` the attributes
  were 150×75, 250×125, 375×187 — note `150 × 1.25 = 187.5` lands as **187**
  because canvas width/height attributes are integers. Floor the product, or
  compare with a ±1 tolerance. Compute the base size from the renderer at
  runtime; never hardcode these numbers.
- **The 20-field cap bites on spgi-100** (88 columns): `Structure` is inside the
  first 20 and appears by default, `Core` and the `R*` columns are not and only
  enter the form when picked in the Fields dialog. The picker's rows are
  canvas-drawn and have no DOM handles — set an ordered subset through
  `setOptions({fieldsColumnNames: [...]})` and use the dialog only for its
  open/close and All/None affordances.
- **Curve promotion happens on ATTACH**: `onTableAttached` promotes
  `rendererSize` from `small` to `normal` when the dataframe already carries a
  column with semType `fit` (#L192-L193). The viewer must therefore be added to
  a table that already has the curve column, with `rendererSize` still at its
  default. Assert the promotion through the canvas size attributes, not by
  reading the property back — a property read-back is set-then-read.
- **Canvas fields have no text channel.** Their honest signals are: the element is
  a `CANVAS` carrying `[column="<Name>"]` (not an `INPUT` holding text), its
  width/height attributes, a runtime non-empty pixel count (opaque non-white px on
  that one field canvas — a draw-floor, never a similarity measure), and a
  zero-console-error floor. A `canvas.toDataURL()` content diff CANNOT prove a
  same-row repaint here: a plain `df.filter` bitset carries no substructure
  highlight, so the same row's molecule redraws byte-identical (measured live
  2026-08-11: 4162→4162, 664 px unchanged). The refdoc's 4454→4706 length change
  was a FIXED-POSITION read — after the filter that position shows a DIFFERENT
  molecule because the card set shrank and shifted — not the same row. Prove a
  same-row repaint by NODE IDENTITY (tag the canvas, confirm it is replaced) plus
  the non-empty floor, never by a content diff.
- **Never pixel-compare the form field with the grid cell**, and never take a
  global non-white pixel count. The form canvas is sized by `rendererSize`
  (× devicePixelRatio) while the grid cell is sized by the grid column width, so
  the two are fundamentally different rasters. Whether the highlight and the
  alignment look right is a visual judgement and lives in `forms-ui.md`.
- **Settle model**: every dataframe subscription is debounced by 50 ms
  (`DG.debounce(stream, 50)`, forms-viewer.ts#L196); filter changes reach the
  cards through `filter.onChanged` / `onMetadataChanged`, both debounced. Live
  reads after a filter write needed ≈1–2.5 s to be stable. Poll for the expected
  state; never sleep a fixed interval as the assert channel.
- **Read a card only after its row is current.** A card built for row `-1`
  renders one empty `div` per field with **no** `column` attribute, so the
  current-row card has zero `[column="…"]` elements until a row is made current.
  Every scenario here sets the current row before reading a field.
- **Scoping**: the Forms viewer's own root carries no `name` attribute; use the
  host `[name="viewer-Forms"]` as the scope and exclude `.temp` — `renderHeader`
  appends a throwaway row-0 form to `document.body` with class `temp`, and
  during that window a global `.d4-multi-form-form` / `[column=…]` query matches
  nodes outside the viewer. Ordinary cards are
  `[name="viewer-Forms"] #vlist .d4-multi-form-form:not(.temp)`.
- **Property panel actuation**: boolean rows (Color Code) toggle with a click on
  `[name="prop-view-<slug>"]`. Choice rows (Renderer Size) need the VALUE element
  clicked first — that click *creates* the `select.property-grid-item-editor-spinner`
  next to the label; a click on the `TR` row does nothing.
- **No automation tokens in step prose**: every reference to properties, the
  dataframe, selectors or page APIs lives in this section only. Step prose uses
  UI-action language and names the control.

### Spec must keep (hard constraints for the spec author)

- **No text-colour assertion, in any form, anywhere in this file.** The channel
  is dead in the product (invalid CSS value) — an assertion on it can only be
  false-green or permanently red.
- **Colour coding must be applied BY BACKGROUND** before Step 1a; a
  text-colour-coded column silently empties the signal and the step degrades to
  "background equals the default", which passes on a broken build.
- **Renderer-size expectations are computed at runtime** from the renderer's
  declared defaults (or the viewer's fallback size) and `window.devicePixelRatio`,
  with `Math.floor` applied where the source applies it. A hardcoded `150` /
  `250` / `375` is a false-green on any other display density.
- **The curve promotion is read through the canvas size, on a freshly attached
  viewer** whose Renderer Size was never written. Reading `rendererSize` back
  after attach is a property echo, not a check of the promotion's effect.
- **Step 5a recomputes selection ∩ filter INSIDE the poll**, on every iteration,
  in the same pass as the card read, and requires that the filter actually
  excluded at least one selected row (an invariant, not a hardcoded count).
  Capturing the expected set once right after the sketcher write races the
  debounced re-render and degrades the check to "cards == selection".
- **Step 5b proves the repaint by NODE IDENTITY, not a content diff.** Tag the
  survivor's Structure canvas before the filter; after the settled re-render assert
  the tagged node is gone, the same row (by id value) holds a fresh canvas without
  the tag, that canvas is drawn above a floor (runtime "> 20 drawn" opaque
  non-white pixels, well under the ~660 a real molecule paints), and the
  console-error floor around the re-render is zero. Do NOT reinstate a same-row toDataURL diff
  (byte-identical without highlight → false red), do NOT compare form vs grid, do
  NOT use a global non-white pixel count. If the node ever stops being replaced,
  no repaint happens at all — fix the atlas expectation, don't weaken the step.
- **Datasets**: demog is `System:DemoFiles/demog.csv`; the molecular dataset is
  `System:AppData/Chem/tests/spgi-100.csv` (table `spgi-100`); the curve dataset
  is `System:DemoFiles/curves.csv`. Those three paths are the only ones this file
  may use — an earlier edition of this section referenced an "SPGI" demo file
  under DemoFiles that does not exist on any server; do not resurrect it.
- **Every `page.evaluate` returns primitives or JSON** — never a DG object (the
  "object reference chain is too long" serialization failure).

### Spec must keep — realization anchors (forms-spec.ts)

- Colour coding applied via `df.col('AGE').meta.colors.setLinear()` (background scheme);
  step guards `gc.isTextColorCoded === false`. Alignment/font set on
  `grid.col('AGE').contentCellStyle`; the colour change fires onMetadataChanged which
  re-renders and re-reads the cell style. Colour compared as [r,g,b]; font normalized
  through a throwaway element's style.font.
- Renderer size actuated through the property panel (prop-view-renderer-size → the created
  select), never by writing rendererSize.
- Renderer-size expectations replicate getRendererSize/getSize at runtime (renderer
  defaultWidth/Height + devicePixelRatio); 120/200/300 only in the no-declared-size branch;
  no hardcoded pixel expectation, picture never measured.
- Curve promotion (4a) read through the curve canvas size only: == normal-step and
  != small-step; rendererSize never read back.
- Row selection (Scenario 4 step 7, Scenario 5 step 2) actuated through `df.selection.set(...)`
  on the dataframe, not the grid's Ctrl+click: the prose prescribes Ctrl+click to pick the rows,
  the spec sets the selection bitset directly. A control-means substitution — the resulting
  selected row set is identical, only the actuation channel differs.
- Substructure filter (Scenario 5) actuated via grok.chem.searchSubstructure → df.filter.init,
  not the Chem sketcher (not drivable headless). Pattern chosen at runtime as first proper
  subset; 4 rows = 2 matched + 2 excluded. 5a recomputes selection∩filter inside the poll and
  requires inter.length < selCount.
- 5b proves the repaint by NODE IDENTITY plus a draw floor, never by a content diff:
  the tagged canvas node is gone, the same row (found by id value, not position) holds
  a fresh untagged canvas, that canvas is non-empty at runtime, and the console-error
  floor around the re-render is zero. A vanished card cannot fake it and a blank canvas
  cannot pass it.

---
{
  "order": 9,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv", "System:DemoFiles/curves.csv"]
}
