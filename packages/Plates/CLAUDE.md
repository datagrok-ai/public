# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Plates** (`@datagrok/plates`) is a Datagrok plugin for **assay plate** management — ingestion, storage, analysis, and visualization of 96-, 384-, and 1536-well experimental plates. It supports importing plates from CSV, Excel, and proprietary instrument formats, persisting plate data to a PostgreSQL database (`plts` schema), running pluggable analyses (DRC, Dose-Ratio, qPCR), and rendering plates as interactive heatmap widgets.

Category: **Visualizations**. Contains a full CRUD app accessible via `Apps | Plates`.

## Build Commands

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run test               # grok test
npm run build-all          # Builds js-api → utils → statistics → this package
npm run link-all           # Links datagrok-api, @datagrok-libraries/utils, @datagrok-libraries/statistics
```

## Key Dependencies

- `@datagrok-libraries/statistics` — `FitSeries`, `FitCurve`, `IFitSeries`, `IFitChartData`, `FIT_FUNCTION_4PL_REGRESSION`, `FitConstants` — used by DRC analysis for curve fitting
- `@datagrok-libraries/utils` — `SchemaEditor` (used in template property editing)
- `@datagrok-libraries/test` — `category`, `test`, `expect` test framework
- `@datagrok/curves` — **devDependency** — `AddStatisticsColumn` function called at runtime via `grok.functions.call('Curves:AddStatisticsColumn', ...)` to extract IC50/Hill/R²/AUC from fitted curves
- `exceljs` — Excel file parsing (loaded lazily via `DG.Utils.loadJsCss`)
- `jstat` — statistical functions (mean, median, std, geomean, etc.)
- `rxjs` — reactive subscriptions (Subject, Observable, operators)
- `wu` — lazy iteration over wells

## Architecture

### Entry Point — `src/package.ts`

`PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:

**Initialization:**
- `_initPlates` — `@init`: registers `PlateCellHandler` and `PlateTemplateHandler` as `ObjectHandler`s; auto-populates demo data if DB is empty
- `autostart` — `@autostart`: calls `createDummyPlateData()`

**Registered Functions:**

| Function | Decorator / Purpose |
|---|---|
| `assayPlatesDemo` | `@demo`, `demoPath: Plates \| Assay Plates` — loads demo xlsx plates |
| `platesFolderPreview` | `@folderViewer` — previews folders containing plate files |
| `previewPlate` | `@fileViewer(txt)` — renders text-format plates |
| `importPlate` | `@fileHandler(txt)` — imports text-format plates |
| `importPlateXlsx` | `@fileHandler(xlsx)` — imports Excel plates with DRC analysis view |
| `previewPlateXlsx` | `@fileViewer(xlsx)` — previews Excel plates |
| `checkExcelIsPlate` | Checks if Excel content contains plate data (size < 1MB, valid plate patterns) |
| `checkCsvIsPlate` | Checks if CSV first line contains well/position headers |
| `checkFileIsPlate` | Checks if text content matches any registered plate reader |
| `platesApp` | `@app` — main Plates application view |
| `platesAppTreeBrowser` | `@appTreeBrowser` — populates the app tree (Create, Search plates/wells/analyses, Templates) |
| `getPlateByBarcode` | Fetches a plate from DB by barcode |
| `createDummyPlateData` | Seeds demo data (templates + plates + analysis results) |

**Exports**: `_package`, auto-generated wrappers from `package.g.ts`

### Core Data Model — `src/plate/plate.ts`

`Plate` class — the central data structure:
- **Fields**: `rows`, `cols`, `data` (DG.DataFrame where each row = one well, columns = layers/properties), `details` (plate-level properties dict), `barcode`, `id`, `plateTypeId`, `plateTemplateId`
- **Well addressing**: `_idx(row, col)` → DataFrame row index; `rowIndexToExcel(dfRow)` → `[row, col]`
- **Outlier management**: `markOutlier(row, col, flag)`, `isOutlier(row, col)`, `markOutliersWhere(field, predicate)`. Emits via `onOutlierChanged` Subject
- **Data access**: `get(field, row, col)`, `values(fields, filter?)` with match/exclude filtering, `fieldValues(field, filter?)`, `getStatistics(field, stats, filter?)`
- **Well iteration**: `wells` getter yields `PlateWell` objects (row, col, all field values)
- **Normalization**: `normalize(field, f, inplace?)` — applies transformation function to create normalized column
- **Validation**: `validateWells(validators)` — runs `IPlateWellValidator[]` across all wells

**Static factory methods** (how plates are created):

| Method | Input |
|---|---|
| `fromGridDataFrame(df, field?)` | 8×12 / 16×24 / 32×48 grid DataFrame (columns = plate columns) |
| `fromTableByRow(df, options?)` | Tidy DataFrame with position column (e.g. "A1") or row/col columns |
| `fromDbDataFrame(df)` | Database query result (row, col, property_id, value columns) |
| `fromCsv(csv, options?)` | CSV string in grid format |
| `fromCsvTable(csv, field)` | CSV string parsed as grid DataFrame |
| `fromExcel(bytes, name?)` | Excel file bytes (uses `findPlatePositions` + `getPlateFromSheet`) |
| `fromPlates(plates[])` | Merges multiple plates (same dimensions) into one with multiple layers |
| `demo()` | Creates a demo plate from `grok.data.demo.wells(96)` |

### Plate Widget — `src/plate/plate-widget/plate-widget.ts`

`PlateWidget` (extends `DG.Widget`) — interactive plate visualization:
- Renders plate as a heatmap grid with circular wells
- **Tabs**: adds `DG.TabControl` for switching layers/analyses
- **Well click**: emits `onWellClick` Subject with `{row, col, dataIndex}`
- **Hover**: highlights hovered well, shows tooltip
- **Validation**: displays cross markers on wells with validation errors
- **Analysis coordinators**: `IAnalysisWidgetCoordinator` interface for bidirectional communication between plate widget and analysis views
- **Factory methods**: `fromPlate(plate)` (detailed view with side panel), `detailedView(plate)` (with well details panel)

### Plate Selection Controller — `src/plate/plate-widget/plate-selection-controller.ts`

`PlateSelectionController` — handles drag-to-select on the plate grid:
- Two modes: `'select'` (toggle selection) and `'passthrough'` (for outlier marking)
- Drag selection with majority-rule toggling (if most selected wells are already selected, unselect all; otherwise select all)
- Emits `onSelectionComplete` with list of selected wells

### Cell Renderer — `src/plate/plate-cell-renderer.ts`

`PlateGridCellRenderer` — renders `Plate` objects in grid cells:
- Cell type: `'Plate'`
- Default size: 240×160
- Renders plate heatmap inside the cell via `PlateWidget`
- Click opens detailed view in the property panel

`PlateCellHandler` — `ObjectHandler` that returns the cell renderer

### Plate Readers — `src/plate/plate-reader.ts`

`PlateReader` — dispatches to format-specific readers:

| Reader | File |  Detection |
|---|---|---|
| `DelfiaEnvision96MaxRluPlateReader` | `readers/delfia-envision-96-max-rlu.ts` | `"Plate information"` + `"DELFIA"` + `"96"` |
| `DelfiaEnvisionByRowPlateReader` | `readers/delfia-envision-by-row.ts` | Delfia Envision by-row format |
| `BmgPherastar96PlateReader` | `readers/bmg-pherastar-96.ts` | `"No. of Channels"` + `"No. of Cycles"` headers |
| `SpectramaxPlateReader` | `readers/spectramax.ts` | `"Plate:"` + `"##BLOCKS="` headers |

### Excel Plates — `src/plate/excel-plates.ts`

- `findPlatePositions(workbook)` — scans Excel sheets for plate patterns (sequential numeric column headers 1–12+, alphabetic row headers A–H+)
- `getPlateFromSheet(position)` — extracts plate data from sheet at detected position

### CSV Plates — `src/plate/csv-plates.ts`

- `parsePlates(df, identifierColNames)` — splits a tidy DataFrame into multiple plates based on identifier columns (e.g., barcode); returns `PlateFile[]` with common properties extracted
- `parsePlateFromCsv(csvContent)` — convenience wrapper using `"Destination Plate Barcode"` as identifier

### Utilities — `src/plate/utils.ts`

| Function | Purpose |
|---|---|
| `numToExcel(n)` / `excelToNum(s)` | Convert between 0-based index and Excel column letters (A, B, ..., Z, AA) |
| `parseExcelPosition(s)` / `toExcelPosition(r, c)` | Parse/format "A3" notation to/from `[row, col]` |
| `toStandardSize([r, c])` | Snap to nearest standard plate size (96/384/1536) |
| `standardPlateSizes` | Map: `{96: [8,12], 384: [16,24], 1536: [32,48]}` |
| `mapFromRow(row)` | Converts a DataFrame row to a display-friendly object (formats numbers, creates outlier checkbox) |
| `jstatStatistics` | Map of stat name → jStat function (mean, median, std, min, max, meansqerr, geomean) |
| `safeLog(n)` | Log10 that returns 0 for non-positive values |

### Well Validators — `src/plate/plate-well-validators.ts`

`IPlateWellValidator` interface with `validate(plate, row, col) → string | null`. Built-in validators:
- Concentration and volume must be both present or both absent
- Concentration must not be negative

---

## Analysis System

### Analysis Manager — `src/plate/analyses/analysis-manager.ts`

`AnalysisManager` — singleton registry of all analysis types:
- Registered analyses: `DrcAnalysis`, `DoseRatioAnalysis`, `QpcrAnalysis`
- `getAnalysis(name)` / `byFriendlyName(name)` — look up by internal or display name
- `init()` — registers analysis properties in the database

### Base Analysis — `src/plate/analyses/base-analysis.ts`

`IPlateAnalysis` interface and `AnalysisBase` abstract class:
- **Key methods**: `getRequiredFields()`, `createView(plate, plateWidget, mappings, ...)`, `queryResults(query)`, `formatResultsForGrid(df)`, `saveResults(plate, resultsDf, params, mappings)`
- **Save flow**: creates an `analysis_run` in DB, saves parameters and results per group
- Each analysis defines `parameters` (inputs) and `outputs` (computed results) as `IAnalysisProperty[]`

### Base Analysis View — `src/plate/analyses/base-analysis-view.ts`

`BaseAnalysisView` — standard UI pattern for analyses:
- Shows mapping panel if required fields are unmapped
- Shows results view when all required fields are mapped
- Uses `AnalysisMappingPanel` for column mapping

### DRC Analysis — `src/plate/analyses/drc/`

`DrcAnalysis` — Dose-Response Curve fitting:
- **Required fields**: Activity, Concentration, SampleID
- **Outputs**: Curve (semType=fit JSON), IC50, Hill Slope, R², Min, Max, AUC
- **Normalization**: auto-normalizes if High/Low Control wells are present
- **Curve fitting**: uses `FIT_FUNCTION_4PL_REGRESSION` from `@datagrok-libraries/statistics`
- **Statistics**: calls `Curves:AddStatisticsColumn` at runtime to extract IC50/Hill/R²/AUC from fitted curves
- `getDoseResponseSeries(plate, options)` — groups wells by sample, builds `FitSeries` sorted by concentration, respects outlier flags

`DrcAnalysisCoordinator` — bidirectional sync between plate widget and curves grid:
- Clicking a well on plate highlights the corresponding point on the curve
- Toggling outlier on the curve updates the plate and regenerates curves
- Listens to `fit-cell-outlier-toggle` custom events from the Curves package

### Dose-Ratio Analysis — `src/plate/analyses/dose-ratio/`

`DoseRatioAnalysis` — multi-curve agonist/antagonist analysis:
- **Required fields**: Agonist_Concentration_M, Antagonist_Concentration_M, Percent_Inhibition, (optional) Agonist_ID, Antagonist_ID
- Groups data by unique antagonist concentrations, fits 4PL curves per group
- Renders all curves overlaid on a single chart

### qPCR Analysis — `src/plate/analyses/qpcr/`

`QpcrAnalysis` — quantitative PCR delta-delta-Ct analysis:
- **Required fields**: Ct (cycle threshold)
- **Required roles** (assigned via well selection UI): Target Gene, Reference Gene, Control, Treated
- **Outputs**: Delta Delta Ct, Fold Change
- Two-phase UI: role assignment view (if roles not assigned) → results grid

---

## Plates App & CRUD Layer

### Plates App — `src/plates/plates-app.ts`

`platesAppView()` — main application view (DG.TableView listing all plates with a `Plate` cell type column)

`initPlatesAppTree(treeNode)` — populates app tree sidebar:
- **Create** — plate creation view
- **Search plates** — search by plate/well/analysis properties
- **Search wells** — search at well level
- **Search analyses** — search saved analysis results
- **Templates** — template management (create, edit, clone, delete)

### CRUD Layer — `src/plates/plates-crud.ts`

Central data access layer (~880 lines). Key types and functions:

**Types:**
- `PlateProperty` — property definition (id, name, type, scope, min, max, choices)
- `PlateTemplate` — template with plate/well properties
- `PlateType` — plate format (rows, cols, max volume)
- `PlateQuery` — search query with plate/well/analysis matchers
- `PropertyCondition` — property + matcher pair for search
- `AnalysisQuery` — analysis-specific search query

**Global state** (cached after `initPlates()`):
- `allProperties` — all registered properties from `plts.properties`
- `plateTemplates` — all templates with their plate/well properties
- `plateTypes` — standard plate formats (96, 384, 1536)

**Key functions:**

| Function | Purpose |
|---|---|
| `initPlates(force?)` | Loads properties, templates, plate types from DB. Caches globally |
| `savePlate(plate, options?)` | Persists plate + well values to DB. Auto-creates missing properties |
| `getPlateById(id)` | Fetches plate from DB by ID |
| `queryPlates(query)` | Searches plates by property conditions. Generates dynamic SQL with EXISTS subqueries |
| `queryWells(query)` | Searches at well level with property filters |
| `queryAnalysesGeneric(query)` | Searches analysis results with group/property filters |
| `createProperty(prop)` | Creates a new property in `plts.properties` |
| `createPlateTemplate(template)` | Creates template + links properties |
| `createAnalysisRun(plateId, type, groups)` | Creates analysis run record |
| `saveAnalysisRunParameter(params)` | Saves analysis input parameter |
| `saveAnalysisResult(params)` | Saves analysis output result per group |
| `getOrCreateProperty(name, type, scope)` | Get or auto-create a property |

**Events**: CRUD events emitted via `events: Subject<CrudEvent>` (before/after create/read/update/delete)

### Matchers — `src/plates/matchers.ts`

Search query matchers that generate SQL WHERE clauses:
- `StringMatcher` — string equality/pattern matching
- `StringInListMatcher` — `IN (...)` list matching
- `NumericMatcher` — supports `=`, `!=`, `>`, `>=`, `<`, `<=`, `in`, `not in`, range (`1-10`), `is null`, `is not null`. Parses user input strings like `">10"`, `"1.5-3.0"`

### Demo Data — `src/plates/plates-demo.ts`

`__createDummyPlateData()` — seeds the database with:
- **Cell counting** template (Imaging device, Status, Sample, Well cell count)
- **Dose-response** template (Project, Stage, Chemist, Biologist, QC, Z-Score, Sample, Role, Concentration, Volume, Activity)
- 12 demo plates per template with randomized data
- DRC analysis runs with fitted curves for dose-response plates

---

## Views & Components

### Create View — `src/plates/views/plates-create-view.ts`

Main plate creation workflow:
- Template/plate type selection
- CSV file import → auto-splits into multiple plates by identifier columns
- `PlateWidget` with interactive well selection and outlier marking
- Analysis tabs (DRC, Dose-Ratio, qPCR) with column mapping
- Uses `PlateStateManager`, `TemplatePanel`, `PlateGridManager`

### Search Views — `src/plates/views/plates-search-view.ts`

`searchPlatesView()` / `searchWellsView()` — property-based search:
- Dynamic form with inputs for template + other properties
- Analysis type selector with per-analysis sub-filters
- Results displayed in DG.TableView grid

### Analysis Search View — `src/plates/views/analyses-search-view.ts`

`searchAnalysesView()` — search saved analysis results:
- Filter by analysis type, group, and output properties
- Results formatted by each analysis type's `formatResultsForGrid()`

### Template Views

- `plates-templates-view.ts` — lists all templates with Edit/Clone buttons
- `plates-schema-view.ts` — template property editor using `SchemaEditor` from `@datagrok-libraries/utils`
- `template-editor-view.ts` — simplified template editor with prototype PlateWidget

### Shared State & Utilities

| File | Purpose |
|---|---|
| `shared/plate-state-manager.ts` | `PlateStateManager` — manages multi-plate state, template selection, identifier columns, analysis mappings. Emits `PlateStateChangeEvent` |
| `shared/types.ts` | `PlateFile`, `TemplateState`, `PlateStateChangeEvent`, `ValidationResult`, `MappingDialogOptions` |
| `shared/scopes.ts` | `MAPPING_SCOPES` constants: `TEMPLATE`, `DRC`, `DOSE_RATIO`, `QPCR` |
| `shared/mapping-utils.ts` | `createDynamicMappingRow()`, `createFormRow()` — UI helpers for column mapping |

### UI Components

| Component | File | Purpose |
|---|---|---|
| `TemplatePanel` | `components/template-panel/template-panel.ts` | Left sidebar: file import, template/plate type selection, plate/well property inputs, identifier column management |
| `PlateGridManager` | `components/plate-grid-manager/plate-grid-manager.ts` | Multi-plate grid showing barcodes, template match status, QC status, common properties |
| `MappingEditor` | `components/mapping-editor/mapping-editor.ts` | Generic column mapping UI: target properties ↔ source columns with dropdown selectors |
| `AnalysisMappingPanel` | `components/analysis-mapping/analysis-mapping-panel.ts` | Analysis-specific mapping panel wrapping `MappingEditor` with validation |
| `ValidationPanel` | `plates-validation-panel.ts` | Template property validation with mapping reconciliation |

### Plate Template Handler — `src/plates/objects/plate-template-handler.ts`

`PlateTemplateHandler` (extends `DG.ObjectHandler`) — context actions for templates:
- Edit, Clone, Delete operations
- Renders template schema view on selection

---

## Database Schema — `databases/plts/0000_init.sql`

PostgreSQL schema `plts` with forward-only migrations:

| Table | Purpose |
|---|---|
| `properties` | Property definitions (name, type, scope=plate/well, constraints, choices, format) |
| `templates` | Plate templates (name, description) |
| `template_plate_properties` | Template ↔ plate property junction (is_required, default_value) |
| `template_well_properties` | Template ↔ well property junction |
| `plate_types` | Plate formats (rows, cols, max_volume). Seeded: 96, 384, 1536 |
| `plates` | Plate instances (plate_type_id, barcode, description) |
| `plate_wells` | Individual wells (plate_id, row, col, details JSONB) |
| `plate_details` | Plate-level property values (value_num, value_string, value_bool, value_jsonb) |
| `plate_well_values` | Well-level property values |
| `analysis_runs` | Analysis execution records (plate_id, analysis_type, groups[]) |
| `analysis_run_parameters` | Analysis input parameters |
| `analysis_results` | Analysis output results per group_combination |
| `semantic_types` | Semantic type definitions (Molecule, Solvent, URL, Image) |
| `property_allowed_values` | Enumerated values for properties |

Seeded properties: Volume (uL), Concentration (uM), Sample, Well Role (Empty/Sample/DMSO/Low Control/High Control).

### SQL Queries — `queries/plates-crud.sql`

~30 parameterized queries for all CRUD operations. Key queries:
- `getPlates`, `getPlateByBarcode`, `getWellValuesByBarcode/ById`
- `getProperties`, `getPlateTypes`, `getPlateTemplates`
- `createProperty`, `createTemplate`, `addTemplatePlateProperty/WellProperty`
- `createAnalysisRun`, `saveAnalysisRunParameter`, `saveAnalysisResult`
- `queryAnalyses`, `queryAnalysesTemplate`, `getAnalysisRunGroups`
- `getPlatesCount` — used to detect empty DB for auto-seeding

---

## Tests — `src/tests/plate-tests.ts`

Tests using `@datagrok-libraries/test`:
- `fromCsvPlate` — CSV grid import, verifies dimensions
- `mergePlates` — merging concentration/layout/readout layers
- `fromFile` — loading from DemoFiles
- `fromExcel` — Excel plate import
- `normalization` — high/low control normalization
- `use case` — full pipeline: merge → normalize → DRC series
- `render` — PlateWidget rendering

---

## CSS Files

| File | Scope |
|---|---|
| `plate/plate-widget/plate-widget.css` | PlateWidget layout, tabs, details panel |
| `plate/analyses/plate-analyses.css` | Analysis view containers, warning/error messages |
| `plate/analyses/drc/drc-analysis.css` | DRC-specific styles (save button) |
| `plates/views/plates-create-view.css` | Create view layout |
| `plates/views/analyses-search-view.css` | Analysis search view layout |
| `plates/views/components/mapping-editor/mapping-editor.css` | Mapping editor rows, icons |
| `plates/views/components/template-panel/template-panel.css` | Template panel sidebar |
| `plates/views/components/template-panel/template-panel-and-mapping.css` | Combined template+mapping styles |
| `plates/views/components/plate-grid-manager/plate-grid-manager.css` | Multi-plate grid styles |

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Function/file handler registration | `src/package.ts` |
| Plate data model, factory methods | `src/plate/plate.ts` |
| Plate interactive widget | `src/plate/plate-widget/plate-widget.ts` |
| Grid cell renderer for plates | `src/plate/plate-cell-renderer.ts` |
| File format readers | `src/plate/plate-reader.ts` + `src/plate/readers/` |
| Excel import | `src/plate/excel-plates.ts` |
| CSV multi-plate parsing | `src/plate/csv-plates.ts` |
| Well position utilities | `src/plate/utils.ts` |
| DRC curve fitting analysis | `src/plate/analyses/drc/drc-analysis.ts` |
| Dose-Ratio analysis | `src/plate/analyses/dose-ratio/dose-ratio-analysis.ts` |
| qPCR analysis | `src/plate/analyses/qpcr/qpcr-analysis.ts` |
| Analysis registration | `src/plate/analyses/analysis-manager.ts` |
| Database CRUD, search queries | `src/plates/plates-crud.ts` |
| SQL query definitions | `queries/plates-crud.sql` |
| DB schema | `databases/plts/0000_init.sql` |
| Search matchers (numeric, string) | `src/plates/matchers.ts` |
| App tree, main view | `src/plates/plates-app.ts` |
| Plate creation workflow | `src/plates/views/plates-create-view.ts` |
| Multi-plate state management | `src/plates/views/shared/plate-state-manager.ts` |
| Template management | `src/plates/views/plates-templates-view.ts` + `plates-schema-view.ts` |
| Demo data seeding | `src/plates/plates-demo.ts` |
| Column mapping UI | `src/plates/views/components/mapping-editor/mapping-editor.ts` |
