# Scenario 3 — Substructure SAR

**Theme:** Working with tabular data + chem analysis
**Length target:** 3–4 min
**Audience lean:** Chemists, comp chem, drug-discovery leads
**Recording risk:** Medium

Mirrors two tutorials in the Tutorials package, ground truth for every step here:
- `packages/Tutorials/src/tracks/chem/tutorials/substructure-search-filtering.ts` (10 steps)
- `packages/Tutorials/src/tracks/chem/tutorials/r-groups-analysis.ts` (14 steps)

## Prompt (multi-turn)

The demo is most natural as a short conversation. Drop them all at once or step through:

> 1. Load the kinase compound set with activity values.
> 2. Filter the molecule column by substructure — indole: `c1ccc2[nH]ccc2c1`.
> 3. Run R-groups analysis using MCS to detect the scaffold.
> 4. In the resulting trellis plot, switch the inner viewer from Pie chart to Histogram and set Value to the activity column.

## Setup

- A kinase compound set with a molecule column (SMILES, auto-detected as `Molecule` semType) and an activity column (e.g. `pIC50` or `In-vivo Activity`). For staging, `packages/Tutorials/files/sar_small-R-groups.csv` is what the R-Groups tutorial uses and works directly; or stage an equivalent kinase CSV.
- `Chem` package installed and registered. Sketcher available.
- Pre-warm Chem once on this server so RDKit WASM is cached.

## Expected sequence

1. **Load.** Agent opens the table. Chem detectors flag the molecule column as `Molecule` (`packages/Chem/src/detectors.js` / sem-type detection).
2. **Substructure filter.** Agent invokes `Chem:Substructure Search` (top-menu `Chem | Search | Substructure Search...`, impl `SubstructureSearchTopMenu` at `packages/Chem/src/package.ts:724`). The function opens the table view's filter group with a chem substructure filter on the molecule column; the sketcher is empty (`DG.WHITE_MOLBLOCK`). To populate it the agent goes through the filter group's API (the same `getFiltersGroup().updateOrAdd({...})` that the function itself uses at `:727`), setting `molBlock` from the indole SMARTS. Grid updates: only indole-containing rows remain, scaffold highlighted in each.
3. **R-Groups Analysis.** Agent invokes `Chem:R-Groups Analysis` (top-menu `Chem | Analyze | R-Groups Analysis...`, impl `rGroupsAnalysisMenu` at `:1095`). The function finds the molecule column by semType and calls `rGroupAnalysis(col)`, which opens a dialog with a sketcher and an **MCS** button. Agent (or presenter) clicks MCS, then OK. The analysis adds `R1`, `R2`, … columns to the grid and opens a **Trellis Plot** (`DG.VIEWER.TRELLIS_PLOT`) with pie charts faceted by R-group positions.
4. **Trellis viz switch.** Agent flips the trellis's inner viewer from Pie chart to **Histogram** via the combo-popup in the trellis's top-left (`v.root.querySelector('.d4-combo-popup')` in the tutorial; the event is `d4-trellis-plot-viewer-type-changed`). Then sets the histogram's **Value** to the activity column via the trellis's Context Panel. Each cell now shows an activity distribution per R-group combination.

## Wow moment

Step 3 lands two beats at once: R-group columns appear in the grid AND the Trellis Plot opens with one chart per substituent combination. Step 4 turns those charts into activity histograms — the SAR map.

## Talking points

- **Chemist:** "Substructure search, then R-groups analysis with MCS. The agent ran the same Chem menu items as the tutorials in the help center."
- **IT lead:** "Each step is a registered platform function — `Chem:Substructure Search`, `Chem:R-Groups Analysis`. The agent invokes them by name through the MCP tool layer; the implementations live in the Chem package."
- **Drug-discovery lead:** "R-groups analysis adds the R-columns and renders the trellis in one operation. The agent kicked it off; everything else is the Chem workflow your chemists already use."

## Works today

- All three referenced top-menu functions exist and are callable through MCP `call_function`:
  - `Chem:Substructure Search` (`packages/Chem/src/package.ts:719–745`) — adds a chem substructure filter to the current table view.
  - `Chem:R-Groups Analysis` (`:1091–1102`) — opens dialog; on OK, adds R-columns and a Trellis Plot.
  - (`Chem:Similarity Search` at `:544–550` is also available if you want a similarity beat as a variant — see below.)
- The Trellis Plot is the native Datagrok `DG.VIEWER.TRELLIS_PLOT`; switching the inner viewer type uses standard viewer events.
- The two tutorials referenced above drive the exact same UI paths the agent will trigger — the demo is the agent doing what those tutorials walk a user through.

## Needs polish (see POLISH.md)

- **Scaffold injection into Substructure Search** (T2.1) — `Chem:Substructure Search` opens the filter with an empty `molBlock`. The agent needs to inject the indole SMARTS. Two options: (a) a Grokky-side skill that opens the search and then calls `tv.getFiltersGroup().updateOrAdd({type: SUBSTRUCTURE, column, molBlock})` with the parsed SMARTS — exactly what `SubstructureSearchTopMenu` itself does internally; (b) bypass the menu and set the filter directly. Either way, wrap as one named skill so the LLM stops improvising.
- **Headless R-Groups Analysis with core** (T2.2) — today `rGroupsAnalysisMenu` opens a dialog and waits for the user to click MCS + OK. The agent needs either to drive the dialog or to call a wrapper that accepts the core SMARTS (or `mode: 'MCS'`) and produces the R-columns + Trellis Plot non-interactively. Without this, the demo has one unavoidable human click; with it, the demo is fully agent-driven.
- **Trellis inner-viewer config** (T1.1 / provisional T2.3) — switching Pie → Histogram and setting Value to an activity column is the most fragile step. A `configureViewer(viewer, {inner: {type: 'Histogram', value: 'pIC50'}})` skill closes this; if a single `viewer.setOptions({inner: {...}})` works it collapses into T1.1.
- **Step-chaining reliability** — three sequential agent actions; partial failures should be retryable per-step.

## Backup plan

- **One human click is acceptable.** If the R-Groups Analysis dialog isn't driven headlessly, the presenter clicks **MCS** and **OK** when the dialog appears. The agent did the open; the human did one click. Voiceover: "MCS finds the common scaffold." Still convincing.
- If the trellis inner-viewer switch (step 4) fails, end on step 3 — R-columns + default pie-chart trellis are already the SAR moment. The histogram is a nice-to-have.
- If substructure injection fails, the presenter draws the indole in the sketcher manually after the agent opens it. Loses the "all from prompt" beat but the rest of the flow is intact.
- Pre-record a fallback take with everything working.

## Variant — add a similarity beat

If you want a fourth Chem capability on screen, insert between steps 2 and 3:

> Open similarity search next to the filtered grid.

Agent invokes `Chem:Similarity Search` (`:548`). The function adds a `Chem Similarity Search` viewer to the table view; the viewer follows the current row by default. The Similarity-Diversity tutorial walks the same UI: `packages/Tutorials/src/tracks/chem/tutorials/similarity-diversity-search.ts`. Adds about 30 seconds and one Chem package surface area to the demo; not load-bearing.

## Note on substructure input

SMARTS in chat is the demo default — fastest, fully scriptable. Drawing in the sketcher is the tutorial path and looks better in person; if Grokky's sketcher-context handoff is wired (presenter draws, agent picks up the molfile), prefer that for the live booth demo and keep the SMARTS prompt for recordings.
