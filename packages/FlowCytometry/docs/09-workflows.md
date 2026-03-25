# 09 — Workflows

## Primary Workflow: Import → Compensate → Gate → Analyze → Export

```
┌─────────────────────────────────────────────────────────────────┐
│  1. IMPORT        │  FCS files → Datagrok DataFrames            │
│  2. COMPENSATE    │  Spillover → Compensation matrix → Apply    │
│  3. GATE          │  Draw gates → Build hierarchy → Evaluate    │
│  4. ANALYZE       │  Statistics + Advanced algorithms           │
│  5. EXPORT        │  PDF report / CSV statistics / GatingML     │
└─────────────────────────────────────────────────────────────────┘
```

---

## Workflow 1: Single FCS File Analysis

### Step 1 — Import

**Trigger:** User drops a `.fcs` file onto Datagrok, or opens from Files browser.

**What happens:**
1. `importFcs()` is called with the file
2. Binary parser reads HEADER → TEXT → DATA
3. Event matrix becomes a `DG.DataFrame` (rows = events, columns = channels)
4. FCS metadata stored as DataFrame tags
5. File opens in a **Cytometry Workspace layout** (custom multi-panel layout):
   - Left panel: **Gating Tree** (empty, shows only "All Events")
   - Center: **CytometryScatterViewer** (FSC-A vs SSC-A by default)
   - Right: **Channel Statistics panel** (summary of all parameters)
   - Bottom: **Histogram strip** (mini histograms for each channel)

### Step 2 — Review Metadata

The channel panel shows:
- `$PnN` (detector name) + `$PnS` (marker name) for each channel
- Instrument name (`$CYT`), acquisition date (`$DATE`), event count (`$TOT`)
- Spillover matrix badge: green if `$SPILLOVER` found, yellow if not

### Step 3 — Compensation

**Sub-workflow A: Inline spillover (most common)**
1. Plugin detects `$SPILLOVER` or `$SPILL` in FCS TEXT
2. User clicks **"Apply Compensation"** button in toolbar
3. Compensation dialog opens showing the spillover matrix as an editable grid
4. User reviews/adjusts coefficients if needed
5. User clicks **"Apply"**
6. Compensated values added as new columns prefixed `comp_` (e.g., `comp_FITC-A`)
7. Viewer automatically switches to showing compensated channels

**Sub-workflow B: External compensation matrix**
1. No `$SPILLOVER` in FCS file (e.g., GvHD data)
2. User clicks **"Load Compensation Matrix"**
3. Uploads a CSV file or selects an already-open DataFrame
4. Same flow from step 5 above

**Sub-workflow C: Calculate from controls**
1. User has a set of single-stain control FCS files open
2. Clicks **"Calculate Spillover from Controls"**
3. Plugin runs server-side Python (FlowKit) to compute the spillover matrix
4. Result fed into dialog from step 5 above

### Step 4 — Transform Selection

- Default: **Logicle** applied to all fluorescence channels
- Scatter (FSC/SSC): always **linear**
- Time: always **linear**
- User can change per-channel via right-click on axis or channel panel

Transform is applied for **display only** — raw and compensated values are preserved in the DataFrame. Gate coordinates are stored in transformed space with a reference to the active transform.

### Step 5 — Gating

**Drawing a gate:**
1. User selects a gate tool from toolbar (polygon, rectangle, quadrant, range)
2. Clicks/draws on the scatter plot
3. Plugin evaluates which events fall inside the gate region
4. A new node appears in the Gating Tree panel with count and %Total
5. The DataFrame's filter/selection is updated to reflect the gated population

**Hierarchical gating:**
1. User clicks on a population in the Gating Tree to select it as the "active population"
2. The scatter plot dims events outside the active population
3. User draws another gate — it becomes a child of the active population
4. Statistics show %Parent (fraction of parent population passing this gate)

**Gate editing:**
- Drag vertices to reshape polygon gates
- Drag the whole gate to move it
- Right-click → Delete Gate

### Step 6 — Statistics

- Statistics panel updates live as gates are drawn
- Shows: Count, %Parent, %Total, MFI, Geometric Mean, CV for each gate × channel
- Click any gate row in the Gating Tree to see its statistics breakdown

### Step 7 — Export

- **Export Statistics CSV:** Flat table, one row per gate, columns = all statistics
- **Export PDF Report:** Plots + gating tree + statistics table
- **Export GatingML:** Gate hierarchy as GatingML 2.0 XML for use in other tools
- **Save as Datagrok Project:** Preserves DataFrame + gating state + viewer layout

---

## Workflow 2: Batch Analysis (Plate / Multi-Sample)

### Step 1 — Import Multiple Files

**Option A: Directory import**
- User connects to a folder (network share, S3) containing FCS files
- Selects all files → imports as a batch
- Each FCS file becomes one DataFrame in a Datagrok Project

**Option B: Plate-aware import**
- FCS files have `$WELLID` / `$PLATEID` keywords
- Plugin detects plate structure and offers to auto-organize by plate layout

### Step 2 — Create/Load Gating Template

**Option A: Create new template from one file**
- Analyze the first/representative FCS file fully (Workflow 1)
- Right-click Gating Tree → **"Save as Template"**
- Template stores gate definitions (coordinates + transforms + channel refs)

**Option B: Import FlowJo .wsp (Phase 2)**
- Drag in a `.wsp` file
- Plugin parses XML, extracts compensation matrices, transformations, and gate hierarchies
- Creates a GatingTemplate from the workspace

**Option C: Import GatingML (Phase 2)**
- Import a `.xml` GatingML 2.0 file as a template

### Step 3 — Apply Template to All Files

1. User clicks **"Apply Template to Batch"**
2. Plugin iterates over all FCS DataFrames in the project:
   - Apply compensation (from each file's `$SPILLOVER` or shared matrix)
   - Evaluate all gates in the template hierarchy
   - Compute statistics for each gate
3. Progress shown in Datagrok status bar (files processed / total)

### Step 4 — Cross-Sample Statistics Table

Results aggregated into a **summary DataFrame**:
```
filename | wellId | treatment | dose | %Lymphocytes | %CD4+ | %CD8+ | MFI_CD4 | ...
```

This DataFrame is a standard Datagrok DataFrame — can be:
- Visualized as a plate heatmap (dose-response view)
- Plotted in Datagrok's standard scatter plot / bar chart
- Fed into dose-response curve fitting
- Exported to CSV

### Step 5 — Plate Heatmap View

If FCS files have `$WELLID` / plate metadata:
1. **"Show as Plate"** button in batch results toolbar
2. Datagrok's existing HTS plate viewer renders selected statistics as a heatmap
3. Color scale = % of selected population; each well = one FCS file

---

## Workflow 3: FlowJo Workspace Migration (Phase 2)

1. User has existing FlowJo `.wsp` workspace file
2. Drops it onto Datagrok alongside the FCS files it references
3. Plugin parses `.wsp`: compensation matrices, biexponential transform parameters, gate hierarchies
4. Reconstructs all gates as Datagrok GatingTree objects
5. Applies to the FCS DataFrames
6. User sees their familiar gating hierarchy reproduced in Datagrok
7. Can then extend, modify, or export as GatingML

---

## Workflow 4: High-Dimensional Analysis (Phase 3)

1. User has a 20-50 parameter FCS file
2. After standard gating (singlets, live cells), user selects a population (e.g., "Live CD45+ cells")
3. Clicks **"Run UMAP"** or **"Run FlowSOM"**
4. Plugin sends the gated population to Datagrok Python scripting:
   - Server computes UMAP coordinates (umap-learn)
   - FlowSOM cluster labels added as a new column
5. Results returned as new columns in the DataFrame: `umap_x`, `umap_y`, `cluster`
6. CytometryScatterViewer offers to show `umap_x` vs `umap_y`, colored by `cluster`
7. User can draw manual gates on the UMAP plot to capture phenotypic subsets

---

## UI Layout Design

### Cytometry Workspace Layout (Datagrok split-pane)

```
┌────────────┬──────────────────────────────┬──────────────┐
│            │                              │              │
│  Gating    │   Main Scatter Plot          │  Channel     │
│  Tree      │   (CytometryScatterViewer)   │  Statistics  │
│  Panel     │                              │  Panel       │
│            │                              │              │
│  - All     ├──────────────────────────────┤              │
│    Events  │                              │              │
│  - Singlets│   Secondary Plot             │              │
│    - Live  │   (Histogram or 2nd scatter) │              │
│      - CD4 │                              │              │
│        CD8 │                              │              │
│            │                              │              │
└────────────┴──────────────────────────────┴──────────────┘
│              Histogram Strip (one per fluorochrome channel) │
└─────────────────────────────────────────────────────────────┘
```

### Viewer Toolbar Items

**CytometryScatterViewer toolbar:**
- X/Y channel selectors (dropdown, shows $PnN + $PnS)
- Scale selector: Logicle / Arcsinh / Log / Linear
- Gate tools: Polygon / Rectangle / Quadrant / Range
- Show/hide: density overlay, gate outlines, event count
- Color by: population / marker / metadata column

**Gating Tree panel:**
- Expandable tree nodes (gate name, count, %parent)
- Right-click menu: Rename, Delete, Export, Set as Active Population
- "Add Boolean Gate" button
- "Save as Template" button

---

## Integration with Other Datagrok Modules

| Integration | How it works |
|---|---|
| **Compound Registration** | Link FCS files to compound IDs via metadata CSV — enables "show %cell killing for compound X" |
| **Plate Reader / HTS** | Plate heatmap of population statistics; same well coordinate system as assay plates |
| **Scripting** | Python/R scripts appear in standard Datagrok script runner; results auto-added to DataFrame |
| **Files Browser** | `.fcs` files show a preview in the file info panel (mini scatter + metadata) |
| **Projects** | All analysis state (gating, compensation, transforms) saved in Datagrok Project format |
| **Data Connections** | Connect to instrument output folders (NFS, S3, SFTP) for live/automated import |
| **Audit Log** | All gating operations, compensation applications, and exports logged with user + timestamp |
