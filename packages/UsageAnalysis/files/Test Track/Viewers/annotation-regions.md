### Annotation regions
Visual overlays on Scatter Plot and Line Chart viewers that highlight, select, and analyze specific areas of the plot.

#### 1. Region creation

1.1 Draw rectangle region
- Open `System:DemoFiles/demog.csv`.
- Add Scatter Plot viewer.
- Open viewer properties (gear icon). In the Tools section, ensure `Lasso Tool` is unchecked.
- Right-click the scatter plot canvas → **Tools** → **Draw Annotation Region**.
- Drag on the canvas to draw a rectangle.
- Expected result:
  - A rectangle region appears on the viewer.
  - Formula Lines dialog opens automatically with the newly created region selected.
  - The dialog's preview Scatter Plot uses the same X/Y axes as the source viewer.
- In the dialog's right-side editor, change: Title, Description, Fill color, Outline color, Outline width, Opacity.
- Expected result: all edits are reflected in the preview and on the source viewer immediately.

1.2 Draw Lasso (polygon) region
- On the same Scatter Plot, open properties → Tools → set `Lasso Tool` = true.
- Right-click → **Tools** → **Draw Annotation Region**.
- Drag a free-form polygon on the canvas.
- Expected result:
  - Polygon region is created.
  - Formula Lines dialog opens with the new region selected; preview axes match the source viewer.

1.3 Create formula region
- Right-click the Scatter Plot → **Tools** → **Formula Lines...**.
- In the dialog toolbar click **ADD NEW** → **Region - Formula Lines**.
- Expected result:
  - A formula region is created with two default boundary formulas referencing the viewer's current X and Y columns.
  - The region appears in the dialog preview and in the grid on the left side of the dialog.

#### 2. Visibility controls

2.1 Show/Hide viewer & dataframe regions separately
- Scatter Plot on demog with at least one **viewer** region (attached to the viewer in the Formula Lines dialog) and one **dataframe** region (same dialog, switch region source selector to the dataframe before adding).
- Open viewer properties → Annotation Regions section.
- Toggle `Show Viewer Annotation Regions`: only viewer region hides/shows.
- Toggle `Show Dataframe Annotation Regions`: only dataframe region hides/shows.
- Expected result: the two toggles act independently; no cross-effect.

2.2 Show/Hide all regions at once
- In viewer properties, uncheck both `Show Viewer Annotation Regions` and `Show Dataframe Annotation Regions`.
- Right-click the canvas → **Tools** → check `Show Annotation Regions`.
- Expected result:
  - Both viewer and dataframe regions become visible.
  - The global toggle overrides the two individual toggles.

#### 3. Region interaction

Prerequisite: Scatter Plot on demog with two overlapping rectangle regions (create two in the Formula Lines dialog, drag so they overlap in the preview, or draw two overlapping rectangles).

3.1 Hover interaction
- Hover over a single region: markers inside are highlighted as the mouse-over group (across all viewers of the table).
- Hover over the overlap area: markers from both regions are highlighted.
- If the region has a Title/Description, the header text is rendered near the region.

3.2 Click interaction
- Click a single region → markers inside become selected.
- Click the overlap area → selection is the union of markers from both regions.
- Ctrl+click on a selected region → markers of that region are removed from the selection.

#### 4. Region editing

4.1 Context menu editing
- Open Line Chart on demog and add one formula region (via right-click → Tools → Formula Lines... → ADD NEW → Region - Formula Lines).
- Right-click on that region → **Edit**.
- Expected result: Formula Lines dialog opens with the region pre-selected for editing.

4.2 Modify region properties
- With the dialog open on an existing region, change each of: Title text, Title color, Description text, Fill color, Outline color, Outline width, Opacity.
- Expected result: each edit is reflected in the dialog preview immediately; after closing the dialog, changes persist on the source viewer.

4.3 Update annotation font (viewer-level)
- Scatter Plot on demog with at least two regions having non-empty Title.
- Open viewer properties → find `Annotation Font` property (Misc / Style section).
- Change font family and size.
- Expected result: region headers on all regions of this viewer update immediately. Formula line labels (for plain formula lines, not region titles) are unaffected.

#### 5. Formula Lines dialog integration

5.1 Preview viewer
- Scatter Plot on demog with at least one region.
- Open the Formula Lines dialog. In the dialog, add one of each: rectangle region, lasso region, formula region.
- Expected result: preview accurately renders each region with the correct boundary shape, fill color, opacity, and (for formula regions) the correct curves.

5.2 Grid representation
- With the same dialog open, inspect the grid on the left side.
- Expected result: the grid lists all regions, with columns for type, title/header, fill color, outline color, formulas (for formula regions) and description.

#### 6. Line Chart specific

6.1 Multi-axis limitation
- Open demog → add Line Chart. In viewer properties set `Y Columns` to more than one numeric column (multi-axis mode).
- Right-click → Tools → try **Draw Annotation Region** and try **Formula Lines...** → ADD NEW → any region type.
- Expected result: annotation region creation is disabled (menu items disabled / drawing has no effect / dialog refuses to add). No partially-drawn region remains on the viewer.

6.2 Single-axis Line Chart region creation
- In the same Line Chart, set `Y Columns` to a single numeric column.
- Create one each: rectangle region, lasso region, formula region (same steps as 1.1–1.3).
- Expected result: all three regions are created successfully and behave identically to regions on a Scatter Plot (hover highlighting, click selection, editing).

---
{
  "order": 29,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
