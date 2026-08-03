# Scenario 2 — Tabular Transform

**Theme:** Working with tabular data
**Length target:** 2–3 min
**Audience lean:** Scientist (bench + comp chem)
**Recording risk:** Low

Mirrors the **Calculated Columns** tutorial (`packages/Tutorials/src/tracks/transform/tutorials/calculated-columns.ts`, 13 steps) and the **Filters**/**Scatter Plot** tutorials, compressed into one agentic flow.

## Prompt

> Add a column **LipE = pIC50 − cLogP**. Filter to rows where LipE > 5. Build a scatter plot of MW vs cLogP, color by LipE, with a regression line. Highlight the top 5 compounds by LipE.

## Setup

- A small medicinal-chem table loaded in the current view. Required columns: `pIC50` (numeric), `cLogP` (numeric), `MW` (numeric), and a molecule column. ~200–1000 rows.
- Suggested source: a curated subset of ChEMBL or a sample CSV bundled with the demo. Pre-load before the take.
- Chem package installed (for the molecule column type detection — purely cosmetic here).
- Empty chat history so the demo starts clean.

## Expected sequence

1. **Calculated column.** Grokky emits `datagrok-exec` that calls `df.columns.addNewCalculated('LipE', '${pIC50} - ${cLogP}')` — signature in `js-api/src/dataframe/column-list.ts:180`: `addNewCalculated(name, expression, type?='auto', treatAsString?, subscribeOnChanges?=true)`, returns `Promise<Column>`. The new column appears in the grid; values populate. Formula DSL with `${ColumnName}` references is documented in `help/transform/add-new-column.md` and demonstrated in `packages/Tutorials/src/tracks/transform/tutorials/calculated-columns.ts`.
2. **Filter.** Agent applies a numeric range filter on `LipE` (`min: 5`) via `tv.getFiltersGroup().updateOrAdd({...})` — the same filter-group API pattern that `Chem:Substructure Search` uses internally (`packages/Chem/src/package.ts:727`), applied here for a numeric column. Row count in status bar drops.
3. **Scatter plot.** `grok.shell.addViewer('Scatter plot', { x: 'MW', y: 'cLogP', color: 'LipE', showRegressionLine: true })`. Viewer type constant: `DG.VIEWER.SCATTER_PLOT = 'Scatter plot'` (`js-api/src/const.ts:691`). The regression-line property is `showRegressionLine` (`js-api/src/interfaces/d4.ts:579`), **not** `regressionLine`.
4. **Top 5 highlight.** Agent sorts grid by LipE descending, then sets the top 5 in the `df.selection` BitSet (`js-api/src/dataframe/data-frame.ts:161`; standard pattern is `df.selection.init(i => sortedRanks[i] < 5)` or per-row `df.selection.set(idx, true)`). Selection is linked across grid and scatter.

## Wow moment

The instant the `LipE` column appears in the grid AND the scatter colors update in lockstep. The audience sees the column populate row-by-row faster than they can read it, and the chart re-renders.

## Talking points

- **Scientist:** "The agent added LipE as a calculated column, set the filter, configured the scatter. The formula stays attached — change it and the column recomputes. Click in the scatter and the grid follows."
- **IT lead:** "The agent calls the same column and viewer operations as the UI. Governance, persistence, and sharing work as they would on a hand-built dashboard."

## Works today

- Adding a calculated column via `Execute` provider works — the LLM writes the JS that calls the DataFrame API.
- Adding scatter plot via `grok.shell.addViewer` works.
- Selection / sort works.

## Needs polish (see POLISH.md)

- **Structured `add-calculated-column` skill** (T1.2) — wrap `df.columns.addNewCalculated(name, expression, type?, treatAsString?, subscribeOnChanges?)` with explicit `{name, formula, type}` signature. Failure mode without the wrapper is silently-wrong formulas (LLM misremembers operator vs function-call syntax).
- **`configure-viewer` skill** (T1.1) — set viewer properties by name with validation. Confirmed gotchas: `showRegressionLine` not `regressionLine`; viewer type string is `'Scatter plot'` not `'Scatterplot'`.
- **Filter API skill** (T1.3) — wrap `TableView.getFiltersGroup().updateOrAdd({type, column, min, max})` for numeric/categorical ranges. Same `updateOrAdd` pattern the platform code uses internally for chem substructure filtering.
- **Select-top-N skill** (T1.3) — sort by column + write top N to `df.selection` BitSet. The agent shouldn't be writing index loops or remembering BitSet semantics.

## Backup plan

- The expression DSL is `df.columns.addNewCalculated(name, formula)`. If the agent's first formula syntax attempt fails (e.g. tries an operator the parser rejects), retry with explicit function names — `Sub(${pIC50}, ${cLogP})` instead of `${pIC50} - ${cLogP}`. Both forms are valid per `help/transform/add-new-column.md`.
- If regression line toggling fails, drop it from the demo and lean on the color encoding for the wow.
- Pre-record a clean take to play back if live execution misbehaves.

## Variant — non-chem audience

Same prompt structure, BMI from height/weight, color by BMI, top 5 by BMI. Use the patient demographics table from the Tutorials package. Useful for non-pharma side-events and general analytics audiences.
