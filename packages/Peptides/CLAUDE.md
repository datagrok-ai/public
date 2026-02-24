# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Peptides** (`@datagrok/peptides`) is a Datagrok plugin for **Structure-Activity Relationship (SAR) analysis** of peptide collections. It detects macromolecule columns automatically, renders amino acids with color-coded monomers, and provides interactive viewers to identify point mutations and residues causing major activity changes.

Category: **Bioinformatics**. Top menu: `Bio | Analyze | SAR...`.

## Build Commands

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run test               # grok test
```

## Architecture

### Entry Point — `src/package.ts`

`PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:
- `initPeptides` — `@init`: loads MonomerWorks, TreeHelper, and monomer library components
- `Peptides` — `@func`: app landing view with Simple/Complex/HELM demo buttons
- `peptidesDialog` (`Bio Peptides`) — `@func`: top-menu SAR dialog (`Bio | Analyze | SAR...`). Validates active table has Macromolecule + numerical columns, shows config dialog, launches analysis
- `peptidesPanel` — `@panel`: property panel widget for Macromolecule columns
- `manualAlignment` — `@panel`: property panel for Monomer semtype, manual sequence editing
- `macromoleculeSarFastaDemo` — `@func`: demo dashboard loading FASTA peptides with −lg scaling + MCL clustering
- `lstPiechartCellRenderer` — `@func`: custom grid cell renderer for Logo Summary Table pie charts

**Registered Viewers** (all `@func` with `role: 'viewer'`):
- `Sequence Variability Map` → `MonomerPosition`
- `Most Potent Residues` → `MostPotentResidues`
- `Sequence Mutation Cliffs` → `MutationCliffsViewer`
- `Logo Summary Table` → `LogoSummaryTable`
- `Sequence Position Statistics` → `SequencePositionStatsViewer`
- `Active peptide selection` → `ClusterMaxActivityViewer`

**Exports**: `_package`, `getMonomerWorksInstance()`, `getTreeHelperInstance()`

### Model — `src/model.ts`

`PeptidesModel` is the central controller for SAR analysis. Stored in `DataFrame.temp['peptidesModel']` as singleton per dataframe.

**Key responsibilities:**
- Manages analysis settings (sequence column, activity column, scaling, viewer toggles, sequence space params, MCL settings)
- Splits aligned sequences into per-position columns (`joinDataFrames`)
- Creates scaled activity column (none/lg/−lg)
- Manages collaborative selection across all viewers (WebLogo, invariant map, mutation cliffs, clusters)
- Builds the Properties panel accordion (Distribution, Mutation Cliffs pairs, Selection)
- Manages viewer lifecycle: add/close dendrogram, sequence space, cluster max activity, logo summary table, MonomerPosition, MostPotentResidues
- Computes and caches `MonomerPositionStats` — per-monomer-per-position statistics
- Handles grid rendering setup: WebLogo column headers, monomer cell renderers, tooltips

**Settings** are stored as a JSON tag on the DataFrame (`TAGS.SETTINGS`). Changing settings triggers selective updates (only affected viewers/computations refresh).

**Important fields:**
- `positionColumns` — columns representing individual positions in split sequences (tagged with `TAGS.POSITION_COL`)
- `webLogoSelection` — current monomer-position selection from WebLogo headers
- `_dm` — cached distance matrix
- `_sequenceSpaceViewer` — scatter plot for sequence space
- `_mclViewer` — MCL clustering viewer

### Types — `src/utils/types.ts`

Key types:
- `PeptidesSettings` — full analysis config: sequenceColumnName, activityColumnName, activityScaling, viewer toggles, sequenceSpaceParams, mclSettings
- `SequenceSpaceParams` — distance function, gap penalties, DBSCAN params, fingerprint type
- `MCLSettings` — MCL clustering params: inflation, threshold, iterations, WebGPU, min cluster size
- `MutationCliffs` — `Map<Monomer, Map<Position, Map<Index, Indexes>>>` — pairwise cliff data
- `Selection` — `{ [positionOrClusterType]: string[] }` — unified selection across viewers

### Constants — `src/utils/constants.ts`

- `COLUMNS_NAMES` — Activity, Monomer, Position, P-Value, Mean difference, Count, Ratio, etc.
- `TAGS` — DataFrame/column tags: SETTINGS, POSITION_COL, ANALYSIS_COL, CUSTOM_CLUSTER, MONOMER_POSITION_MODE, etc.
- `SEM_TYPES` — Monomer, MacromoleculeDifference
- `SCALING_METHODS` — none, lg, −lg
- `SUFFIXES` — Viewer prefixes for namespacing (LST, MP, MPR, WL)

## Viewers (`src/viewers/`)

| File | Class | Purpose |
|---|---|---|
| `sar-viewer.ts` | `SARViewer` (abstract) | Base for MonomerPosition and MostPotentResidues. Manages mutation cliffs computation, invariant map stats, selection, mode switching (Invariant Map / Mutation Cliffs) |
| `sar-viewer.ts` | `MonomerPosition` | Horizontal heatmap: monomers × positions. Two modes: Invariant Map (colored by aggregated stats) and Mutation Cliffs (circle size=count, color=mean diff). Has monomer search/filter |
| `sar-viewer.ts` | `MostPotentResidues` | Vertical viewer: one row per position showing the most potent monomer with mean difference, p-value, count, ratio |
| `logo-summary.ts` | `LogoSummaryTable` | Per-cluster summary grid: WebLogo, activity distribution histogram, members count, mean difference, p-value. Supports original + custom clusters, filtering small clusters |
| `mutation-cliffs-viewer.ts` | `MutationCliffsViewer` | Line chart of mutation cliffs at a selected position. Splits by series column, syncs selection with main dataframe |
| `cluster-max-activity-viewer.ts` | `ClusterMaxActivityViewer` | Scatter plot: cluster size vs max activity per cluster. Draws threshold lines, supports auto-selection of top quadrants |
| `position-statistics-viewer.ts` | `SequencePositionStatsViewer` | Box/violin plot of numerical values grouped by monomers at a selected position. Configurable motif overhang |

## Widgets (`src/widgets/`)

| File | Purpose |
|---|---|
| `peptides.ts` | Main SAR launch UI: `analyzePeptidesUI` (config dialog) + `startAnalysis` (creates model, adds viewers, starts analysis). This is the primary entry point for all SAR analysis |
| `distribution.ts` | Activity distribution panels with histograms. Supports breakdown by monomers, positions, or clusters |
| `mutation-cliffs.ts` | Mutation Cliffs panel: pairs grid + unique sequences grid. Filtering by monomer, shift-click multi-select |
| `selection.ts` | Selection summary grid mirroring selected rows with WebLogo headers and monomer-position stats |
| `manual-alignment.ts` | Text area for editing aligned peptide sequences with apply/reset |
| `settings.ts` | Settings dialog accordion: General (activity, scaling), Viewers (toggle each), Columns (aggregation), Sequence Space (distance, DBSCAN), MCL params |

## Utilities (`src/utils/`)

| File | Purpose |
|---|---|
| `algorithms.ts` | Core computations: `findMutations` (mutation cliff detection via web workers), `calculateMonomerPositionStatistics` (per-position stats), `calculateClusterStatistics`, `calculateMutationCliffStats` |
| `cell-renderer.ts` | Custom grid cell renderers: `renderMutationCliffs` (colored circles), `renderInvariantMap` (colored rectangles), `renderWebLogo` (stacked letter logos), `setWebLogoRenderer` (column header WebLogos with mouse events), `LSTPieChartRenderer` |
| `misc.ts` | Helpers: `scaleActivity`, `modifySelection` (Shift/Ctrl logic), `highlightMonomerPosition`, `initSelection`, `getSelectionBitset`, `isSelectionEmpty`, `extractColInfo`, `expandGrid` |
| `statistics.ts` | Stats engine: `calculateStats` (t-test, mean difference, p-value, ratio for monomer vs rest), `getAggregatedValue`, `getAggregatedValues`, `getFrequencies` |
| `tooltips.ts` | Rich tooltips for monomer-position cells: histogram + stats table + aggregated columns |
| `parallel-mutation-cliffs.ts` | `MutationCliffsCalculator` class: distributes pairwise comparisons across Web Workers for parallel computation |
| `types.ts` | Type definitions (see Types section above) |
| `constants.ts` | Constants (see Constants section above) |

### Web Worker — `src/workers/mutation-cliffs-worker.ts`

Receives a chunk of the upper-triangular pairwise comparison matrix. For each pair, checks activity delta threshold and max mutation count. Returns qualifying pairs (position, index1, index2).

### Helpers — `src/peptideUtils.ts`

`PeptideUtils` static class: lazily loads and caches monomer library (`IMonomerLib`) and sequence helper (`ISeqHelper`). Call `loadComponents()` before using.

### Demo — `src/demo/fasta.ts`

`macromoleculeSarFastaDemoUI`: loads sample CSV, configures FASTA/PT notation, scales IC50 with −lg, launches SAR with MCL clustering.

## Analysis Flow

1. User opens table with Macromolecule column → `initPeptides` loads MonomerWorks + TreeHelper
2. User triggers SAR via top menu (`Bio | Analyze | SAR...`) or column panel → `analyzePeptidesUI` shows config dialog
3. `startAnalysis` creates `PeptidesModel`, splits sequences into position columns, scales activity
4. Model adds viewers (MonomerPosition, MostPotentResidues, LogoSummaryTable, etc.) to the TableView
5. Viewers compute stats independently but share selection through the model's collaborative filtering
6. Property panel accordion shows Distribution, Mutation Cliffs pairs, and Selection based on current viewer/selection state

## Tests (`src/tests/`)

| File | Categories |
|---|---|
| `core.ts` | Start analysis (simple/complex), save/load project |
| `viewers.ts` | Viewer rendering, selection, interaction |
| `widgets.ts` | Settings dialog, distribution, mutation cliffs widgets |
| `model.ts` | Model initialization, settings changes, viewer management |
| `table-view.ts` | Grid interactions, selection sync, WebLogo behavior |
| `misc.ts` | Algorithm unit tests (mutation cliff detection) |
| `benchmarks.ts` | Performance tests: mutation cliffs, cluster stats, monomer-position stats at 5k–200k sequences |
| `utils.ts` | Test constants (`TEST_COLUMN_NAMES`) |

## Key Dependencies

- `@datagrok-libraries/bio` — macromolecule utilities, sequence splitting, monomer palettes, WebLogo, MonomerWorks
- `@datagrok-libraries/ml` — distance matrices, dimensionality reduction (UMAP/tSNE), MCL clustering, macromolecule distance functions
- `@datagrok-libraries/math` — DBSCAN worker, WebGPU utilities
- `@datagrok-libraries/statistics` — statistical functions
- `@datagrok-libraries/utils` — BitArray, UI helpers (`u2.appHeader`)
- `@datagrok-libraries/tutorials` — tutorial framework
- `wu` — lazy iteration
- `cash-dom` — jQuery-like DOM manipulation
- `uuid` — unique ID generation
- `rxjs` — reactive subscriptions

## Quick Lookups

| Looking for... | Check first |
|---|---|
| App/viewer/panel registration | `src/package.ts` |
| Analysis launch flow | `src/widgets/peptides.ts` → `startAnalysis` |
| Central model, settings, viewer management | `src/model.ts` |
| All type definitions | `src/utils/types.ts` |
| Column names, tags, scaling enums | `src/utils/constants.ts` |
| SAR heatmap viewers (MonomerPosition, MostPotentResidues) | `src/viewers/sar-viewer.ts` |
| Cluster summary viewer | `src/viewers/logo-summary.ts` |
| Mutation cliffs line chart | `src/viewers/mutation-cliffs-viewer.ts` |
| Custom grid cell rendering (WebLogo, heatmap cells) | `src/utils/cell-renderer.ts` |
| Statistics (t-test, p-value, mean diff) | `src/utils/statistics.ts` |
| Mutation cliff detection (parallel) | `src/utils/algorithms.ts` + `src/utils/parallel-mutation-cliffs.ts` |
| Selection logic (Shift/Ctrl) | `src/utils/misc.ts` → `modifySelection` |
| Tooltips for monomer-position cells | `src/utils/tooltips.ts` |
| Settings dialog UI | `src/widgets/settings.ts` |
| Demo data files | `files/aligned.csv`, `aligned_2.csv`, `aligned_3.csv` |
| Benchmark test data | `files/tests/` (5k–200k .d42 files) |
