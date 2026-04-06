### Annotation regions 
Visual overlays on scatter plot and line chart viewers that allow you to highlight, select, and analyze specific areas of data visualization.

#### 1. Region creation
  
1.1 Draw rectangle region
- Open demog dataset. 
- Add Scatterplot viewer. Check that lassoTool = false
- Right-click > Tools > Draw Annotation Region.
- Draw a rectangle in the viewer. 
- Observe Formula Lines Editor opens automatically with the region selected.
- Change all format and description properties.
  - Expected result:
  - Rectangle region is created.
  - Formula Lines Editor opens with the newly created region selected, viewer`s preview has proper axis set. 
  - All format and description property changes should be applied correctly.
  
1.2 Draw Lasso (polygon) region
- Set lassoTool = true
- Right-click > Tools > Draw Annotation Region.
- Draw free-form polygon in viewer.
- Verify Formula Lines Editor opens automatically.
  - Expected result:
  - Lasso region is created correctly. 
  - Editor opens with region selected. Viewer`s preview has proper axis set.

1.3 Create formula region
- Right-click > Tools > Formula lines. 
- Formula Lines dialog opens.
- Click ADD NEW > Region - Formula lines.
- Verify default formulas reference active columns.
  - Expected result:
  - Formula region created with two boundary formulas. Region appears in preview and grid.

#### 2. Visibility controls

2.1 Show/Hide viewer & dataframe regions separately
- Create Scatterplot with a couple of viewer and dataframe annotation regions.
- Open viewer`s properties. Find Annotation regions section.
- Toggle **Show Viewer Annotation Regions** checkbox.
- Toggle **Show Dataframe Annotation Regions** checkbox.
  - Expected result:
  - Viewer regions toggle independently from dataframe regions.
  - No cross-effect between the two toggles.

1.2 Show/Hide all regions at once
- In viewer`s properties uncheck both Annotation regions checkboxes.
- Right-click on the viewer > Tools > check Show Annotation Regions checkbox.
  - Expected result:
  - Both viewer and dataframe regions are toggled together. 
  - Additionally check that individual toggles are overridden by the global toggle.

#### 3. Region interaction

3.1 Hover interaction
- Demog. Scatterplot. At least two overlapping regions exists.
- Hover over a single region.
- Hover over region intersection.
  - Expected result:
  - Single region - markers inside highlighted as mouse-over group.
  - Regions intersection - markers from all intersected regions are highlighted.
  - Headers and descriptions displayed correctly (if present).

3.2 Click interaction
- Click on a single region. Click selects markers inside region.
- Click on an intersection of regions. Intersection selects union of markers.
- Ctrl+click on a selected region. Ctrl+click deselects markers from that region.

#### 4. Region editing

4.1 Context menu editing
- Demog. Line chart: One region is added.
- Right-click an annotation region. Select **Edit**.
  - Expected result:
  - Formula lines dialog opens with region loaded.

4.2 Modify region properties
- Open Formula lines dialog.
- Update:
  - Title text and color
  - Description text
  - Region and outline colors
  - Outline width
  - Opacity
- Expected result:
- Region updates in preview immediately. Changes persist after closing editor.

4.3 Update annotation font
- Demog. Scatterplot. One region is added.
- Right-click an annotation region. Change annotationFont property.
- Verify headers of all regions.
  - Expected result:
  - All headers update globally. No effect on formula line text.

#### 5. Formula Lines dialog integration

5.1 Preview viewer
- Demog. Scatterplot. One region is added.
- Open Formula Lines dialog. Add various regions.
  - Expected result:
  - Preview viewer accurately displays regions with correct boundaries, fill, opacity, and formulas.

5.2 Grid representation
- Demog. Scatterplot. At least one region is added.
- Open Formula Lines dialog. 
- Inspect grid entries.
  - Expected result:
  - Grid lists all regions with: type, header, fill/outline, formulas (if applicable), description.

#### 6. Line Chart specific
6.1 Multi-axis limitation
- Demog. Create multi-axis line chart.
- Attempt to create or edit regions.
  - Expected result:
  - Annotation regions disabled. No partial regions remain.

6.2 Single-axis Line Chart region creation
- Use line chart with a single Y-axis.
- Create rectangle, lasso, and formula regions.
  - Expected result:
  - All regions created correctly and behave like in scatter plot.



---
{
  "order": 29,
  "datasets": ["System:DemoFiles/demog.csv","System:DemoFiles/SPGI.csv"]
}
