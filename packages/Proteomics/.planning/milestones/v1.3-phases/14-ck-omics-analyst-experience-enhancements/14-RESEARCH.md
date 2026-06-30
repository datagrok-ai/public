# Phase 14: CK-omics Analyst-Experience Enhancements — Research

**Researched:** 2026-05-28
**Domain:** Datagrok TS package extension — viewer polish, info-panel extension, parser post-processing, external REST integration (Ensembl), pathway-result filtering
**Confidence:** HIGH (all critical findings verified against in-tree Datagrok JS-API source, official Ensembl docs, and the locked CK-omics reference functions)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Scope packaging & sequencing
- **D-01:** **One phase, all in.** R1–R5 + G1–G4 land together in Phase 14. Volcano polish (G1+G2+G3) folds naturally into R2 since both touch the same viewer; G4 (Filters scoping) is a 1-task fix.
- **D-02:** **MUST-have for BP DMD/WT client deliverable:** G1 (volcano visual parity), R1 (gene-label resolution), R2 (live counters). Planner sequences these first. Nice-to-haves: R3, R4, R5, G2, G3, G4.

#### Volcano polish (G1 + G2 + G3 + R2)
- **D-03:** **Default point labels:** auto-label top ~15–20 proteins by the active significance metric (adj.p-value default; p-value when toggled). Recompute on metric/filter/selection change.
- **D-04:** **Direction-label semantics from `proteomics.groups` tag.** Render legend categories as `"Enriched in <group1>"` / `"Enriched in <group2>"` / `"Not significant"`. **Color map LOCKED:** magenta = group1 (numerator) `0xFFFF00FF`, cyan = group2 (denominator) `0xFF00FFFF`, gray = NS `0xFFAAAAAA`. Single code path for both Spectronaut Candidates and non-Candidates flows.
- **D-05:** **Protein search via native Datagrok Filters viewer**, NOT a custom textarea. Filter-viewer matches wire to `df.selection.set(...)` — NOT `df.filter` — so highlighted points stay visible in the cloud. **Unifies with G4:** add `Gene name` / `Protein ID` + `Comparison` to one Filters viewer, explicitly exclude `Flags`.
- **D-06:** **Live counter overlay** anchored bottom-right of the volcano (floating). Shows Total / Enriched-in-g1 / Enriched-in-g2 / NS, and by-location breakdown in location mode. Recomputes on `df.filter`, `df.selection`, and viewer-property changes.

#### Volcano polish — Claude's discretion (planner picks)
- Contrast-aware title synthesized from `proteomics.de_method` + group names.
- Axis label rewriting on metric toggle.
- G2 dialog state preload mechanism (sp.getOptions vs df.tags vs cache).
- G3 progress widget choice (TaskBarProgressIndicator vs ui.setUpdateIndicator).

#### Gene-label resolution (R1)
- **D-07:** **Eager resolution at parse time, all rows.** New `Display Name` column populated by parser post-parse step before any viewer renders.
- **D-08:** **Inline provenance markers in `Display Name` string** (`*` grouped, `†` predicted/reclassified). Raw ID kept in `Source ID` column. Single source-of-truth — no separate Provenance column.
- **D-09:** **Ensembl REST POST `/lookup/id` batched by species (~1000 IDs per request).** Group IDs by detected prefix → one POST per species → cross-session cache via `grok.userSettings` (mirror Phase-13 D-02 UniProt cache, including `__schema_v` invalidation). **Three-level fallback:** Ensembl name → Ensembl description → raw ID. All traffic via `grok.dapi.fetchProxy()`.
- **D-10:** **Duplicate descriptions: disambiguate by appending source ID.** Emit `grok.shell.warning('N duplicate gene names disambiguated')` once after import.

#### Per-protein panel + correlation + smart filter
- **D-11:** **Per-group quantities in UniProt panel: compact bar chart with mean ± SD whiskers.** Extends `src/panels/uniprot-panel.ts`.
- **D-12:** **Group-mean correlation scatter at `Proteomics | Visualize | Group-Mean Correlation…`** — native Datagrok scatter viewer with `Numerator Mean` / `Denominator Mean` columns. Color by significance category (D-04 palette). Pearson + Spearman as title annotation. `createGroupMeanCorrelation()` factory in `viewers/`. Distinct from QC dashboard.
- **D-13:** **Smart pathway filtering: port CK-omics `apply_smart_pathway_filtering` verbatim.** Match-for-match port. GO (BP/CC/MF), KEGG, Reactome, WP.
- **D-14:** **Smart pathway filtering default-on with dialog toggle to disable.**

### Claude's Discretion
- Exact wording/placement of axis label rewrites and the counter overlay.
- Storage mechanism for G2 dialog state preload.
- Progress indicator widget choice for G3.
- Ensembl batch size cap if 1000 turns out too aggressive.
- Whether `Display Name` column replaces or augments the existing `Gene name` column.
- How `*` vs `†` marker rules map to CK-omics cases.
- Per-bar styling in the per-group quantity panel.

### Deferred Ideas (OUT OF SCOPE)
- Custom volcano JsViewer.
- Per-protein waterfall chart in UniProt panel (alternative to D-11 bar chart).
- Configurable per-source pathway cap exposed in Enrichment dialog (alternative to D-14 binary toggle).
- Group-mean correlation fold into QC dashboard (alternative to D-12 standalone menu).
- Provenance-driven color or filter in viewers (alternative to D-08 inline-marker).
- Methodology-check / power-user mode for raw enrichment (alternative to D-14 default-on).
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| R1 | Resolve predicted/uncharacterized gene IDs (ENS*, MGP_, LOC, RGD, AABR) via Ensembl; provenance markers; duplicate disambiguation | §"Standard Stack > Ensembl REST"; §"Code Examples > R1 verbatim port"; §"CK-omics Reference Functions" |
| R2 | Live filter-aware counters on volcano (up/down/NS + by-location) | §"Architecture Patterns > Live Counter Overlay"; §"Code Examples > Counter overlay wiring" |
| R3 | Click-point → per-protein group-quantity bars in UniProt panel | §"Architecture Patterns > UniProt Panel Extension"; §"Code Examples > Bar chart SVG" |
| R4 | Group-mean correlation scatter (Pearson + Spearman) | §"Standard Stack > Statistics"; §"Code Examples > Correlation factory" |
| R5 | Smart hierarchical pathway filtering (drop generic GO parents) | §"Code Examples > R5 verbatim port"; §"CK-omics Reference Functions" |
| G1 | Volcano visual parity (title, axes, default labels, magenta/cyan, search) | §"Architecture Patterns > Volcano polish"; §"Common Pitfalls > Axis labels" |
| G2 | Dialog state preload from current viewer state | §"Architecture Patterns > Dialog State Preload" |
| G3 | Color→Location progress indicator + perf | §"Architecture Patterns > Progress UX" |
| G4 | Filters viewer `columnNames` strict scoping | §"Architecture Patterns > Filters Viewer Scoping"; §"Common Pitfalls > columnNames auto-include" |
</phase_requirements>

## Summary

Phase 14 is a *layered enhancement* phase on top of the Phase 13 parity foundation. It is **not** a greenfield phase. Every requirement extends an existing file:

- **Volcano polish (G1+G2+G3+R2+D-03/D-04/D-06):** extends `src/viewers/volcano.ts` (the `createVolcanoPlot` factory + `recomputeVolcano` recompute function) and the `volcanoOptions` handler in `src/package.ts:274-355`. No fork — Phase 13 explicitly chose the property-toggle-on-one-viewer model.
- **Gene-label resolution (R1+D-07..D-10):** new shared module `src/utils/gene-label-resolver.ts` invoked from all five parsers (`maxquant-parser.ts`, `spectronaut-parser.ts`, `spectronaut-candidates-parser.ts`, `fragpipe-parser.ts`, `generic-parser.ts`). Mirrors the Phase 13 `src/analysis/subcellular-location.ts` cache architecture (already in tree: `grok.dapi.userDataStorage` + `__schema_v` + bounded-concurrency worker pool + per-chunk POST via `grok.dapi.fetchProxy`).
- **Filters viewer unification (D-05 + G4):** extends `dockComparisonFilterIfMultiContrast` in `src/package.ts:117-131`. The root-cause fix for G4 is the typed `filters` property (`Array<{type, column}>`), not the higher-level `columnNames` allowlist that Phase 13 observed silently extending itself.
- **UniProt panel + per-group bars (R3+D-11):** extends `renderUniProtWidget` in `src/panels/uniprot-panel.ts:109-176` with a compact SVG bar chart inserted after the GO Terms section.
- **Group-Mean Correlation (R4+D-12):** new factory `createGroupMeanCorrelation` in `src/viewers/group-mean-correlation.ts`. Inline Pearson + Spearman (the @datagrok-libraries/statistics package only exposes Kendall's tau; verified in `libraries/statistics/src/correlation-coefficient.ts:68`).
- **Smart pathway filtering (R5+D-13+D-14):** verbatim TypeScript port of CK-omics `apply_smart_pathway_filtering` invoked after `gGOSt` in `src/analysis/enrichment.ts`, gated by a new checkbox in `showEnrichmentDialog`.

**Primary recommendation:** Sequence per D-02. Wave 1 ships the BP-DMD-vs-WT-blocking trio: G1 (volcano visual parity), R1 (gene-label resolution), R2 (live counters). Wave 2 ships R3+R4+R5+G2+G3+G4. The gene-label resolver should land before G1 — the volcano's default top-N labels (D-03) want `Display Name` over `Gene name`.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Ensembl `/lookup/id` POST (R1/D-09) | API (external REST via `grok.dapi.fetchProxy`) | Browser cache (`grok.userSettings`) | External REST owns the gene-name resolution; cache lives client-side to survive sessions. Mirrors Phase 13 UniProt pattern. |
| Parse-time gene-label resolution (D-07) | Browser (parser modules run in-browser) | API (Ensembl batch fetch) | Parsers already run in-browser today; the new shared resolver is a post-parse step that calls Ensembl once per species. |
| Display Name / Source ID columns (D-08/D-10) | DataFrame (in-memory columns added by parsers) | — | New semantic types; pure data-shape extension. |
| Volcano counter overlay (D-06/R2) | Browser (DOM/RxJS subscriptions on `df.onFilterChanged` / `df.selection.onChanged` / `sp.onPropertyChanged`) | — | Pure client-side computation; no I/O. |
| Default top-N labels (D-03/R2) | Browser (scatter prop binding via `sp.props.labelColumnNames` + `showLabelsFor: 'Selected'`) | DataFrame (drive top-N via selection bitset) | Uses native Datagrok scatter label engine; selection state is the lever. |
| Filters viewer unification (D-05/G4) | Browser (DG.Viewer.filters config with typed `filters` array) | DataFrame (`df.selection.set` for search-match highlight) | Filters viewer is a platform viewer; selection is the highlight-not-hide mechanism. |
| UniProt per-group bars (R3/D-11) | Browser (inline SVG in panel widget) | DataFrame (read per-group columns synchronously) | Panel is an existing widget; bars are pure DOM, no new viewer. |
| Group-Mean Correlation viewer (R4/D-12) | Browser (DG.ScatterPlotViewer factory) | DataFrame (derived `Numerator Mean` / `Denominator Mean` columns) | Same shape as volcano — native scatter + derived columns. |
| Smart pathway filter (R5/D-13) | Browser (post-g:GOSt result transform in `analysis/enrichment.ts`) | — | Pure transform of the response array before `buildEnrichmentDf`. |
| Enrichment dialog checkbox (D-14) | Browser (ui.input.bool in `showEnrichmentDialog`) | — | Existing dialog gains one input. |

## Standard Stack

### Core (all already in the package — no new deps)

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `datagrok-api` | ^1.25.0 [VERIFIED: package.json] | Platform primitives (`DG.Viewer.*`, `ui.*`, `grok.dapi.*`) | Platform-bundled; only TS API for Datagrok packages |
| `@datagrok-libraries/statistics` | ^1.2.12 [VERIFIED: package.json] | Existing stats (Kendall's tau, multiple-tests, box-plot) | Already used by other Proteomics modules |
| `rxjs` | ^6.5.5 [VERIFIED: package.json] | Subscription pattern for viewer event wiring | Platform-bundled webpack external |
| `wu` | ^2.1.0 [VERIFIED: package.json] | Iterator helpers (existing utility) | Already used elsewhere in the codebase |

### Supporting (in-tree, used by name)

| Module | Path | Purpose | When to Use |
|--------|------|---------|-------------|
| `grok.dapi.fetchProxy` | `datagrok-api/grok` [CITED: js-api/src/dapi.ts] | CORS-safe external REST | EVERY external call (Ensembl, UniProt, g:Profiler) |
| `grok.dapi.userDataStorage` | `datagrok-api/grok` | Cross-session per-user KV cache | Mirror Phase 13 UniProt cache for R1 Ensembl cache |
| `DG.Viewer.filters` | `js-api/src/viewer.ts:235` [VERIFIED: in-tree] | Native Filters viewer | D-05 / G4 — unified search + scoping |
| `IFiltersSettings.filters` | `js-api/src/interfaces/d4.ts:1708` [VERIFIED: in-tree] | Typed per-column filter spec array | D-05 / G4 root-cause fix — explicit per-column filter objects, NOT `columnNames` allowlist |
| `IFiltersSettings.columnNames` | `js-api/src/interfaces/d4.ts:1706` [VERIFIED: in-tree] | Higher-level columnNames allowlist | Phase 13 observed silently extending; planner verifies whether `filters` array supplants this |
| `DG.RowSet.Selected` | `js-api/src/interfaces/d4.ts:325` [VERIFIED: in-tree] | `showLabelsFor` enum value for selection-driven labels | D-03 default top-N labels |
| `DG.TaskBarProgressIndicator.create(label)` | already used at `package.ts:337` [VERIFIED: in-tree] | Top-bar progress UX | G3 Color→Location progress |
| `df.onFilterChanged` (Observable) | `js-api/src/dataframe/data-frame.js:389` [VERIFIED: in-tree] | Fires on filter mutation | D-06 counter recompute trigger |
| `df.onSelectionChanged` (Observable) | `js-api/src/dataframe/data-frame.js:387` [VERIFIED: in-tree] | Fires on selection mutation | D-06 counter recompute trigger |
| `viewer.onPropertyChanged` | DG.Viewer base [CITED: pattern used in Datagrok docs] | Fires on viewer-prop edit | D-06 counter trigger for metric/color toggle |
| `sp.getOptions()` | `js-api/src/viewer.ts:158` [VERIFIED: in-tree] | Returns `{look: {colorColumnName, yColumnName, ...}}` | G2 dialog state preload source of truth |
| `df.meta.formulaLines` | already used at `volcano.ts:178-188` [VERIFIED: in-tree] | Threshold/diagonal lines | D-12 correlation diagonal reference |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Inline Pearson/Spearman in `viewers/group-mean-correlation.ts` | `@datagrok-libraries/statistics` Pearson/Spearman | Statistics lib only exposes `kendallsTau` [VERIFIED: `correlation-coefficient.ts:68` is the sole export]; adding Pearson/Spearman would require an upstream PR. Inline is the established pattern. |
| Ensembl REST POST `/lookup/id` (D-09 locked) | MyGene.info, biomaRt, manual lookup CSV | LOCKED by D-09; CK-omics already validated this path. Alternatives would diverge from client contract. |
| `ui.div` + canvas for per-group bars (D-11) | New `DG.JsViewer` for bar chart | Overkill — panel widget is the right surface; SVG is theme-portable without extra plumbing. |
| `df.filter` for search-match (rejected) | `df.selection` (LOCKED D-05) | LOCKED — selection keeps NS cloud visible, matches CK-omics behavior. |

**No new npm dependencies for Phase 14.** All capabilities use platform primitives or existing in-tree libraries. The Package Legitimacy Audit section is consequently empty — no install step in this phase.

**Version verification:** Verified the four runtime deps via `package.json` inspection above. `datagrok-api` is a webpack external — version pin in `package.json` is informational only; runtime resolves to the host platform.

## Package Legitimacy Audit

| Package | Registry | Age | Downloads | Source Repo | slopcheck | Disposition |
|---------|----------|-----|-----------|-------------|-----------|-------------|
| — | — | — | — | — | — | No new packages installed in Phase 14 |

**Packages removed due to slopcheck [SLOP] verdict:** none — phase installs nothing.
**Packages flagged as suspicious [SUS]:** none — phase installs nothing.

*Phase 14 is a pure source-code modification phase: extends existing files, adds new files under `src/utils/` / `src/viewers/`, and adds two new SEMTYPE constants. No `npm install` step in any task. No legitimacy gate needed.*

## Architecture Patterns

### System Architecture Diagram

```
Parser entry (5 sources: maxquant / spectronaut / spectronaut-candidates / fragpipe / generic)
  |
  v
parseXText() — returns DG.DataFrame with Gene Name / Protein ID columns + semTypes
  |
  v
[NEW Phase 14] resolveGeneLabels(df)  ──── src/utils/gene-label-resolver.ts (R1/D-07..D-10)
  |    ├── scan Gene Name column for predicted-pattern prefixes (ENS*, MGP_, LOC, RGD, AABR)
  |    ├── group IDs by species (ENSG=human, ENSMUSG=mouse, ENSRNOG=rat, ENSDARG=zebrafish, MGP_=mouse)
  |    ├── cache lookup via grok.dapi.userDataStorage[STORE_GENE_LABELS] keyed by ID
  |    ├── misses → batched POST https://rest.ensembl.org/lookup/id (≤1000 IDs / batch)
  |    │           → grok.dapi.fetchProxy (CORS-safe)
  |    │           → bounded concurrency (reuse runWithConcurrency from subcellular-location.ts:219)
  |    ├── three-level fallback: external_name → description → raw ID
  |    ├── markers: '*' (grouped) appended, '†' (predicted/reclassified) suffix
  |    └── writes two new columns: Display Name (semType=DISPLAY_NAME), Source ID (semType=SOURCE_ID)
  |
  v
grok.shell.addTableView(df)
  |
  v
[EXTENDED Phase 14] dockComparisonFilterIfMultiContrast(tv, df)  ──── src/package.ts:117 (D-05/G4)
  |    ├── builds typed filters array: [{type:'categorical', column:'Comparison'},
  |    │                                {type:'free-text', column:'Display Name'},
  |    │                                {type:'free-text', column:'Source ID'}]
  |    ├── docks DG.Viewer.filters(df, {filters: [...]}) — REPLACES the columnNames-only allowlist
  |    └── subscribes onFilterChanged → convert filter-matches to selection (D-05 wiring)
  |
  v
Volcano flow:  createVolcanoPlot(df) ─── src/viewers/volcano.ts (G1/D-03/D-04/D-06)
  |    ├── ensureNegLog10Column / ensureDirectionColumn (existing — direction-string change)
  |    ├── DIRECTION_COLORS = {[`Enriched in ${g1}`]: 0xFFFF00FF,
  |    │                       [`Enriched in ${g2}`]: 0xFF00FFFF,
  |    │                       'Not significant':    0xFFAAAAAA}
  |    ├── sp.props.labelColumnNames = ['Display Name']  (was 'Gene name')
  |    ├── sp.props.showLabelsFor = 'Selected' + driver: top-N indices → df.selection
  |    │   (preserves hover via separate mouseOverRow render path — see Pitfall 4)
  |    ├── synthesized title via sp.setOptions({title: `Volcano Plot: ${g1} vs ${g2}`})
  |    ├── DOM overlay for axis labels (Y label: -Log10(Q-value) / -Log10(p-value),
  |    │   X label: Log2 Fold Change) — platform has no xAxisCustomTitle prop
  |    └── attaches Live Counter Overlay (D-06):
  |         floating div anchored bottom-right via sp.root.style.position=relative
  |         subscribed to df.onFilterChanged, df.onSelectionChanged, sp.onPropertyChanged
  |         recompute from df.filter bits + active direction/location column
  |
  v
Volcano Options dialog re-open (G2):
  |    ├── source-of-truth: sp.getOptions().look.colorColumnName + .yColumnName
  |    ├── derive 'metric' from yColumnName via mapping (negLog10P + active metric tag on df)
  |    └── derive 'colorDim' from colorColumnName (significance = direction col, location = LOCATION_COL)

UniProt panel click (R3/D-11):  uniprotPanel(proteinId)  ─── src/panels/uniprot-panel.ts:183
  |    ├── existing fetch via fetchUniProtData → renderUniProtWidget
  |    └── [NEW] append per-group bar chart section AFTER GO Terms (line 173):
  |         reads getGroups(df) + df.col(intensityCol).get(currentRowIdx) per group column
  |         renders SVG with 2 bars (mean ± SD whiskers), text below: "mean=x.xx  SD=x.xx"
  |         empty state: "No per-group quantities available for this protein"

Group-Mean Correlation (R4/D-12):  Proteomics | Visualize | Group-Mean Correlation…
  |    └── createGroupMeanCorrelation(df) ─── src/viewers/group-mean-correlation.ts (NEW)
  |         ├── ensureFreshFloat(df, 'Numerator Mean')  — mean of group1.columns per row
  |         ├── ensureFreshFloat(df, 'Denominator Mean')— mean of group2.columns per row
  |         ├── DG.Viewer.scatterPlot(df, {x, y, color: direction})
  |         ├── df.meta.formulaLines.addLine(y = x diagonal at #888888)
  |         ├── inline pearson(num, den) + spearman(rank(num), rank(den))
  |         └── sp.setOptions({title: `Group-Mean Correlation — r=${r.toFixed(2)} (Pearson), ρ=${rho.toFixed(2)} (Spearman)`})

Enrichment dialog (R5/D-13/D-14):  showEnrichmentDialog(df) ─── src/analysis/enrichment.ts:371
  |    ├── [NEW] ui.input.bool('Apply smart pathway filter', {value: true})
  |    └── post-gGOSt transform:  applySmartPathwayFilter(gostResults, maxPerSource=15)
  |         ├── partition by source: GO:BP vs other
  |         ├── sort GO:BP by p_value ascending
  |         ├── for each GO:BP row: drop if name matches generic-parent terms AND
  |         │   already-kept set has specific-child term
  |         ├── cap at maxPerSource per source
  |         └── if kept < total, set df.setTag('proteomics.enrichment_smart_filtered','true') +
  |             render banner above grid (see Copywriting in 14-UI-SPEC.md)
```

### Recommended Project Structure

```
src/
├── utils/
│   ├── proteomics-types.ts           # ADD: DISPLAY_NAME, SOURCE_ID, NUMERATOR_MEAN, DENOMINATOR_MEAN to SEMTYPE
│   ├── column-detection.ts           # No change
│   └── gene-label-resolver.ts        # NEW (R1) — shared Ensembl resolver + cache, called from all 5 parsers
├── parsers/
│   ├── maxquant-parser.ts            # +1 line: call resolveGeneLabels(df) post-process
│   ├── spectronaut-parser.ts         # +1 line: call resolveGeneLabels(df) post-process
│   ├── spectronaut-candidates-parser.ts  # +1 line: call resolveGeneLabels(df) post-process
│   ├── fragpipe-parser.ts            # +1 line: call resolveGeneLabels(df) post-process
│   ├── generic-parser.ts             # +1 line in dialog onOK: call resolveGeneLabels(df)
│   └── shared-utils.ts               # No change
├── analysis/
│   ├── enrichment.ts                 # MODIFY (R5/D-13/D-14): dialog checkbox + smart-filter transform
│   ├── experiment-setup.ts           # No change
│   └── subcellular-location.ts       # No change (reference pattern for R1 cache)
├── viewers/
│   ├── volcano.ts                    # MODIFY (G1/D-03/D-04/D-06): direction strings + counter overlay
│   ├── group-mean-correlation.ts     # NEW (R4/D-12)
│   ├── enrichment-viewers.ts         # MAYBE MODIFY (R5 banner — if banner lives in this layer)
│   └── ...                           # No change
├── panels/
│   └── uniprot-panel.ts              # MODIFY (R3/D-11): append per-group bars after GO Terms
├── package.ts                        # MODIFY: dockComparisonFilterIfMultiContrast (D-05/G4),
│                                     #         volcanoOptions handler (G2 preload + G3 progress already partly there)
└── tests/
    ├── gene-label-resolver.ts        # NEW: fixture-driven Ensembl parse + cache + species detection
    ├── volcano.ts                    # MODIFY: assert magenta/cyan/gray + new direction strings + counter math
    ├── enrichment.ts                 # MODIFY: assert smart filter dropped parents + cap per source
    └── group-mean-correlation.ts     # NEW: assert Pearson/Spearman correctness + diagonal line
detectors.js                          # MIRROR: add detectors for DISPLAY_NAME, SOURCE_ID, NUMERATOR/DENOMINATOR_MEAN
```

### Pattern 1: Filters viewer with typed per-column filters (D-05/G4 root-cause fix)

**What:** The Phase-13 observation — `DG.Viewer.filters(df, {columnNames: [...]})` silently included `Flags` — is consistent with the platform's documented filter-auto-detection behavior: boolean columns get auto-combined into a combined-bool filter independent of `columnNames`. The robust scoping mechanism is the typed `filters` property (an array of `{type, column}` objects), not `columnNames`.

**When to use:** D-05 (unified protein search + comparison filter) and G4 (drop Flags).

**Example:**
```typescript
// Source: js-api/src/interfaces/d4.ts:1708 — `IFiltersSettings.filters: Array<{[index:string]: any}>`
// Source: datagrok.ai/help/visualize/viewers/filters — filter type names: 'categorical', 'free-text', 'histogram'
const filtersViewer = DG.Viewer.filters(df, {
  filters: [
    {type: 'categorical', column: 'Comparison (group1/group2)'},
    {type: 'free-text', column: 'Display Name'},
    {type: 'free-text', column: 'Source ID'},
  ],
  showBoolCombinedFilter: false,    // explicit — keeps auto-combined-bool filter from re-introducing Flags
  showHeader: true,
  showSearchBox: true,
});
tv.dockManager.dock(filtersViewer, DG.DOCK_TYPE.RIGHT, null, 'Filters', 0.3);
```

**Verification gate the planner should add:** when both `filters` and `columnNames` are provided, observe whether the platform deduplicates or stacks them. Local test: create a DataFrame with `Comparison`, `Display Name`, `Source ID`, **AND** a `Flags` boolean column, dock the viewer with only the typed `filters` array, then read `viewer.getOptions().look` and assert `Flags` is absent. If it still appears, fall back to `viewer.setOptions({filters: [...], columnNames: ['Comparison (group1/group2)', 'Display Name', 'Source ID']})` post-create.

### Pattern 2: Live counter overlay anchored to a viewer (D-06)

**What:** A floating DOM overlay rendered into the viewer's `root` element with `position: absolute; bottom: 8px; right: 8px`. Subscribed to three observables.

**When to use:** D-06 volcano counter.

**Example:**
```typescript
// Source: in-tree volcano.ts pattern + js-api d4.ts onFilterChanged/onSelectionChanged observables
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';

function attachCounterOverlay(sp: DG.ScatterPlotViewer, df: DG.DataFrame): rxjs.Subscription[] {
  // Container — sp.root is the viewer root; relative-position the parent so absolute child anchors correctly.
  const host = sp.root;
  host.style.position = host.style.position || 'relative';
  const overlay = ui.divV([], {
    style: {
      position: 'absolute', bottom: '8px', right: '8px',
      padding: '8px 16px', background: 'rgba(255,255,255,0.85)',
      borderRadius: '4px', fontSize: '0.9em', pointerEvents: 'none',
      zIndex: '5',  // above canvas, below dialogs (platform dialogs are z-index 1000+)
    },
  });
  host.appendChild(overlay);

  const recompute = () => {
    const total = df.filter.trueCount;
    const colorCol = df.col(sp.props.colorColumnName);
    const counts = new Map<string, number>();
    if (colorCol) {
      const cats = colorCol.categories;
      for (const cat of cats) counts.set(cat, 0);
      // Iterate via getRawData for string columns is not viable — use get() guarded by isNone.
      for (let i = 0; i < df.rowCount; i++) {
        if (!df.filter.get(i)) continue;
        const cat = colorCol.get(i) as string;
        counts.set(cat, (counts.get(cat) ?? 0) + 1);
      }
    }
    // Re-render in place — clear children + rebuild
    overlay.replaceChildren();
    overlay.appendChild(ui.divText('Visible Proteins', {style: {fontWeight: 'bold'}}));
    overlay.appendChild(ui.divText(`Total: ${total.toLocaleString()}`));
    for (const [cat, n] of counts) overlay.appendChild(ui.divText(`${cat}: ${n.toLocaleString()}`));
  };

  recompute();
  const subs: rxjs.Subscription[] = [
    df.onFilterChanged.pipe(debounceTime(50)).subscribe(recompute),
    df.onSelectionChanged.pipe(debounceTime(50)).subscribe(recompute),
    sp.onPropertyChanged.pipe(debounceTime(50)).subscribe(recompute),
  ];
  return subs;
}
```

**Note:** Phase 13 enrichment-viewers.ts:9 uses a module-level `activeSubscriptions` for cleanup. Mirror that — store subscriptions and dispose on viewer detach OR keep them tied to the viewer object via a closure.

### Pattern 3: Default top-N point labels driven by selection (D-03)

**What:** Use `sp.props.showLabelsFor = 'Selected'` + drive `df.selection` from the top-N significance metric. Hover labels still work because the scatter shows ALL labelable rows under MouseOverRow regardless of `showLabelsFor`.

**When to use:** D-03 default top-N labels on volcano + D-12 same logic on correlation scatter.

**Example:**
```typescript
// Source: js-api/src/interfaces/d4.ts:322 (RowSet enum: 'Selected', 'MouseOverRow', etc.)
function setTopNLabels(df: DG.DataFrame, sp: DG.ScatterPlotViewer, n: number = 15): void {
  // Top-N by ASCENDING metric over filtered set; tiebreak by |log2FC| desc.
  const yCol = df.col('negLog10P')!;  // larger = more significant
  const fcCol = df.col('log2FC')!;
  const yRaw = yCol.getRawData() as Float32Array | Float64Array;
  const fcRaw = fcCol.getRawData() as Float32Array | Float64Array;
  const idx: number[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (df.filter.get(i) && yRaw[i] !== DG.FLOAT_NULL && fcRaw[i] !== DG.FLOAT_NULL)
      idx.push(i);
  }
  idx.sort((a, b) => {
    const dy = yRaw[b] - yRaw[a];  // higher -log10(p) first
    if (dy !== 0) return dy;
    return Math.abs(fcRaw[b]) - Math.abs(fcRaw[a]);
  });
  const topN = idx.slice(0, n);
  // CK-omics labels these in addition to MouseOverRow — drive via selection.
  df.selection.setAll(false, false);
  for (const i of topN) df.selection.set(i, true, false);
  df.selection.fireChanged();
  sp.props.labelColumnNames = ['Display Name'];
  sp.props.showLabelsFor = 'Selected';  // values from RowSet enum: 'All' | 'Filtered' | 'Selected' | 'SelectedOrCurrent' | ...
}
```

**Caveat:** D-05 ALSO uses `df.selection` for search-match highlights. The two are not exclusive — search adds to selection on top of the top-N seed. Planner picks the merge policy. Recommended: top-N labels seed selection on initial render and on metric/filter toggle; search-match union-adds. Clearing search restores top-N seed via a `recompute()` call.

### Pattern 4: Cross-session per-user cache (Phase 13 carry-forward, applied to R1)

**What:** Mirror the EXACT shape used by `src/analysis/subcellular-location.ts:STORE` and `SCHEMA_V`. The reason: the Volcano Options handler at `package.ts:325` already does `grok.dapi.userDataStorage.get(SUBCELL_STORE)` to drive cold-vs-warm-cache UX wording. R1's resolver should expose the same shape so analogous UX can drive a similar pre-OK toast for re-imports.

**When to use:** D-09 Ensembl cache.

**Example:**
```typescript
// Source: in-tree src/analysis/subcellular-location.ts:165-253
export const STORE_GENE_LABELS = 'proteomics-gene-labels';
const SCHEMA_V_GENE_LABELS = '14-r1-1';  // bump only if the resolver contract changes
const SCHEMA_KEY = '__schema_v';

async function loadGeneLabelCache(): Promise<Record<string, GeneLabel>> {
  try {
    const raw = (await grok.dapi.userDataStorage.get(STORE_GENE_LABELS)) ?? {};
    if (raw[SCHEMA_KEY] !== SCHEMA_V_GENE_LABELS)
      return {[SCHEMA_KEY]: SCHEMA_V_GENE_LABELS};
    return raw as Record<string, GeneLabel>;
  } catch {
    return {[SCHEMA_KEY]: SCHEMA_V_GENE_LABELS};
  }
}

async function flushGeneLabelCache(cache: Record<string, any>, fetched: Record<string, any>): Promise<void> {
  try {
    await grok.dapi.userDataStorage.put(STORE_GENE_LABELS,
      {...cache, ...fetched, [SCHEMA_KEY]: SCHEMA_V_GENE_LABELS});
  } catch (e: any) {
    console.warn(`Gene-label cache write failed: ${e?.message ?? e}`);
  }
}
```

**Provenance note:** the existing `grok.dapi.userDataStorage` is fresh-instance-per-access (see project memory `reference_dapi_fresh_instance_patch`). Do not patch the instance. Pattern in `subcellular-location.ts:289` (single timer-driven flush, not in worker bodies) is the correct race-free model.

### Pattern 5: Dialog state preload from viewer (G2)

**What:** Read `sp.getOptions()` at dialog open time and seed every input from the resolved values, with `df.tags` as the secondary source and the last-input closure cache as the absolute fallback.

**When to use:** G2 — Volcano Options re-open shows current state.

**Example:**
```typescript
// Source: js-api/src/viewer.ts:158 — Viewer.getOptions() returns {id, type, look: {...}}
function readVolcanoState(df: DG.DataFrame, sp: DG.ScatterPlotViewer): VolcanoState {
  const opts = sp.getOptions();
  const look = (opts as any).look ?? {};
  // colorColumnName: 'direction' → significance; 'Subcellular Location' → location
  const colorDim: 'significance' | 'location' =
    look.colorColumnName === 'Subcellular Location' ? 'location' : 'significance';
  // yColumnName is always 'negLog10P' (stable name); active metric is on a df.tag we add.
  const metric = (df.getTag('proteomics.volcano_metric') as MetricKind) ?? 'adj.p-value';
  return {colorDim, metric};
}
```

**Counterpart:** when the OK handler runs `recomputeVolcano(...)`, persist the metric to the df tag: `df.setTag('proteomics.volcano_metric', metric)`. This makes the dialog idempotent across re-opens AND across page reloads of the same DataFrame.

### Pattern 6: Inline SVG bar chart in a panel (D-11)

**What:** Pure SVG, no canvas, no new viewer. Renders into the existing panel container.

**When to use:** R3/D-11 per-group bars in UniProt panel.

**Example:**
```typescript
// Source: extends src/panels/uniprot-panel.ts:109-176 render path
function renderPerGroupBars(df: DG.DataFrame, accession: string): HTMLElement | null {
  const groups = getGroups(df);
  if (!groups) return null;
  // Find this protein's row by accession
  const idCol = findColumn(df, SEMTYPE.PROTEIN_ID, ['primary protein id', 'protein id']);
  if (!idCol) return null;
  let rowIdx = -1;
  for (let i = 0; i < df.rowCount; i++) {
    const raw = idCol.get(i) as string | null;
    if (raw && parseAccession(raw) === accession) { rowIdx = i; break; }
  }
  if (rowIdx < 0) return null;

  const stats = [groups.group1, groups.group2].map((g) => {
    const vals = g.columns.map((n) => df.col(n)?.get(rowIdx) as number)
      .filter((v) => typeof v === 'number' && !isNaN(v) && v !== DG.FLOAT_NULL);
    const n = vals.length;
    if (n === 0) return {name: g.name, mean: NaN, sd: NaN, n: 0};
    const mean = vals.reduce((s, v) => s + v, 0) / n;
    const variance = n > 1 ? vals.reduce((s, v) => s + (v - mean) ** 2, 0) / (n - 1) : 0;
    return {name: g.name, mean, sd: Math.sqrt(variance), n};
  });
  if (stats.every((s) => s.n === 0))
    return ui.divText('No per-group quantities available for this protein');

  // SVG render — 2 bars, fixed width, mean ± SD whiskers, magenta/cyan
  const W = 200, H = 120, BAR_W = 24, GAP = 8, PAD = 24;
  const maxVal = Math.max(...stats.map((s) => (isNaN(s.mean) ? 0 : s.mean + (isNaN(s.sd) ? 0 : s.sd))));
  const scale = (v: number) => H - PAD - (v / maxVal) * (H - 2 * PAD);
  const COLORS = ['#FF00FF', '#00FFFF'];  // D-04 palette
  const svgNs = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNs, 'svg');
  svg.setAttribute('width', String(W));
  svg.setAttribute('height', String(H));
  stats.forEach((s, i) => {
    const x = PAD + i * (BAR_W + GAP * 2);
    if (s.n === 0) return;
    const y = scale(s.mean);
    const rect = document.createElementNS(svgNs, 'rect');
    rect.setAttribute('x', String(x));
    rect.setAttribute('y', String(y));
    rect.setAttribute('width', String(BAR_W));
    rect.setAttribute('height', String(H - PAD - y));
    rect.setAttribute('fill', COLORS[i] ?? '#AAAAAA');
    svg.appendChild(rect);
    if (s.n > 1 && !isNaN(s.sd)) {
      // Whiskers from mean-sd to mean+sd
      const whiskerX = x + BAR_W / 2;
      const yTop = scale(s.mean + s.sd);
      const yBot = scale(s.mean - s.sd);
      const line = document.createElementNS(svgNs, 'line');
      line.setAttribute('x1', String(whiskerX));
      line.setAttribute('x2', String(whiskerX));
      line.setAttribute('y1', String(yTop));
      line.setAttribute('y2', String(yBot));
      line.setAttribute('stroke', '#333');
      line.setAttribute('stroke-width', '1');
      svg.appendChild(line);
    }
    const label = document.createElementNS(svgNs, 'text');
    label.setAttribute('x', String(x + BAR_W / 2));
    label.setAttribute('y', String(H - 4));
    label.setAttribute('text-anchor', 'middle');
    label.setAttribute('font-size', '10');
    label.textContent = `${s.name} (n=${s.n})`;
    svg.appendChild(label);
  });
  const container = ui.divV([svg]);
  stats.forEach((s) => {
    if (s.n > 0) {
      const txt = ui.divText(`${s.name}: mean=${s.mean.toFixed(2)}  SD=${s.sd.toFixed(2)}`);
      txt.style.fontSize = '0.85em';
      container.appendChild(txt);
    }
  });
  return container;
}
```

The DataFrame the panel needs to read from is the protein DataFrame, but the panel only receives a `proteinId` string. Solution: walk `grok.shell.tables` to find the table whose `getGroups()` is set AND whose primary-protein-id column resolves the accession; if multiple tables match, prefer the one whose `tv.dataFrame === grok.shell.tv?.dataFrame`. This logic should live in a tiny `findHostDataFrameForProtein(accession)` helper.

### Pattern 7: Smart pathway filter — verbatim CK-omics port (R5/D-13)

**What:** Match the CK-omics two-phase filter exactly: GO:BP sorted by p_value ascending → for each row check generic-parent flag + already-kept specific-child → keep if it survives → cap at maxPerSource; other sources just take `.head(maxPerSource)`.

**When to use:** R5 / D-13.

**Reference:** `~/Downloads/ck/CKomics_tool2.py:4685-4736`. Generic parent terms (literal CK-omics list): `'localization'`, `'cellular component organization'`, `'transport'`, `'cellular process'`, `'biological process'`, `'metabolic process'`. Specific child terms (literal CK-omics list): `'actin'`, `'vesicle'`, `'endocytosis'`, `'cytoskeleton'`. Default `max_per_source = 15`.

**Example port:**
```typescript
// Source: ~/Downloads/ck/CKomics_tool2.py:4685-4736 — apply_smart_pathway_filtering
// LOCKED CLIENT CONTRACT — match-for-match port. Do not re-derive heuristics.
const GENERIC_PARENT_TERMS = [
  'localization', 'cellular component organization', 'transport',
  'cellular process', 'biological process', 'metabolic process',
];
const SPECIFIC_CHILD_TERMS = ['actin', 'vesicle', 'endocytosis', 'cytoskeleton'];

export interface SmartFilterStats {
  total: number; kept: number; droppedParents: number; cappedAtN: number;
}

export function applySmartPathwayFilter(
  results: GostResult[], maxPerSource = 15,
): {kept: GostResult[]; stats: SmartFilterStats} {
  if (results.length === 0)
    return {kept: [], stats: {total: 0, kept: 0, droppedParents: 0, cappedAtN: maxPerSource}};

  const goBp = results.filter((r) => r.source === 'GO:BP');
  const other = results.filter((r) => r.source !== 'GO:BP');

  let droppedParents = 0;
  const filteredGoBp: GostResult[] = [];
  if (goBp.length > 0) {
    const sorted = [...goBp].sort((a, b) => a.p_value - b.p_value);
    for (const r of sorted) {
      const name = r.name.toLowerCase();
      const isGeneric = GENERIC_PARENT_TERMS.some((g) => name.includes(g));
      if (isGeneric) {
        const hasSpecificChild = filteredGoBp.some((existing) =>
          SPECIFIC_CHILD_TERMS.some((s) => existing.name.toLowerCase().includes(s)));
        if (hasSpecificChild) { droppedParents++; continue; }
      }
      filteredGoBp.push(r);
      if (filteredGoBp.length >= maxPerSource) break;
    }
  }

  // Other sources: top-N by FDR ascending, per source (CK-omics groups them in one .head(maxPerSource)
  // but the GSD ROADMAP wording — "cap per source" — implies per-source. CK-omics line 4731
  // does .head(maxPerSource) on the combined other_data; mirror VERBATIM to honor D-13 contract.)
  const otherSorted = [...other].sort((a, b) => a.p_value - b.p_value);
  const otherCapped = otherSorted.slice(0, maxPerSource);

  const kept = [...filteredGoBp, ...otherCapped].sort((a, b) => a.p_value - b.p_value);
  return {
    kept,
    stats: {
      total: results.length,
      kept: kept.length,
      droppedParents,
      cappedAtN: maxPerSource,
    },
  };
}
```

**Where to invoke:** after `gGOSt(...)` returns inside `runEnrichmentPipeline` at `src/analysis/enrichment.ts:344-348`, BEFORE `buildEnrichmentDf`. The dialog checkbox controls whether the transform runs at all; if disabled, the raw results flow through unchanged.

### Pattern 8: Gene-label resolver — verbatim CK-omics port (R1/D-07..D-10)

**What:** Mirror CK-omics `improve_gene_labels_with_ensrnog_marking` (~/Downloads/ck/CKomics_tool2.py:895-1059) end-to-end. Critical contract elements:

1. **Predicted prefix list (in order):** `['ENSRNOG','ENSMUSG','LOC','ENSG','ENSDARG','ENSRNO','MGP_','RGD','AABR']` — line 898-909.
2. **Ensembl IDs eligible for /lookup/id POST:** ONLY those starting with `'ENS'` OR `'MGP_'` — line 933: `[pid for pid in all_predicted_ids if pid.startswith(('ENS', 'MGP_'))]`. LOC/RGD/AABR/ENSRNO are NCBI/RGD/Affymetrix — Ensembl will not resolve them. CK-omics keeps them with `†` marker only.
3. **Resolution priority (lines 800-804 in get_ensembl_annotations):** `external_name → display_name → description.split('[')[0].strip() → gene_id`. The chosen name is the `best_name`.
4. **Acceptance gate (lines 977-988):** a candidate name is REJECTED if it equals `clean_label` (the raw ID) OR if it starts with any predicted-pattern prefix. If rejected, fall through to description.
5. **Description cleanup (lines 1062-1095, `extract_readable_description`):** strip `[organism]` suffix, regex-replace `^Predicted to (enable|be involved in|be located in|be part of|be)\s+` (case-insensitive), take first sentence by splitting on `.`, strip `_RAT`/`_MOUSE`/`_HUMAN` suffixes, capitalize first letter.
6. **ProteinDescriptions fallback (lines 997-1005):** if Ensembl yields nothing, look for a `ProteinDescriptions` column on the row and apply `extract_readable_description` to its value.
7. **Marker rules (line 1011):** `final_label = improved_name + ('*' if has_asterisk else '') + '†'`. The `†` is appended UNCONDITIONALLY to every predicted-improved label. The `*` is preserved from the input only if it was already there.
8. **Original ID column (lines 938-940, 1014):** stored in `Original_ID` for cases where the candidate name was rejected and we kept the raw ID. Empty string when a real gene name was found.
9. **Duplicate warning (lines 1047-1057):** group by final `GeneLabel` for renamed proteins; warn if any label maps to >1 protein. CK-omics text: `"WARNING: {n} descriptions map to multiple proteins"`. D-10 modifies this — instead of warn-only, **disambiguate** with `"{name} ({source_id})"` suffix.

**Critical species-detection mapping:**

| Prefix | Species | Ensembl species code | Resolves via /lookup/id? |
|--------|---------|----------------------|--------------------------|
| `ENSG` | Human | `homo_sapiens` | YES |
| `ENSMUSG` | Mouse | `mus_musculus` | YES |
| `ENSRNOG` | Rat | `rattus_norvegicus` | YES |
| `ENSDARG` | Zebrafish | `danio_rerio` | YES |
| `MGP_` | Mouse Genomes Project | `mus_musculus` | YES |
| `ENSRNO` | Rat (legacy transcript/protein form) | `rattus_norvegicus` | maybe — keep but treat as fallback |
| `LOC` | NCBI Entrez ID | (NCBI, not Ensembl) | NO — `†` only, raw ID kept |
| `RGD` | Rat Genome Database | (RGD, not Ensembl) | NO — `†` only, raw ID kept |
| `AABR` | Affymetrix/predicted | (not Ensembl) | NO — `†` only, raw ID kept |

The /lookup/id endpoint does NOT need a `species` parameter — IDs are globally unique within Ensembl. The species code helps if the planner later wants to add a `species` filter for safety, but CK-omics omits it (line 779-787). [CITED: rest.ensembl.org/documentation/info/lookup_post — `species` is optional]

**Why this matters:** D-09 said "POST per species" — that's CK-omics's *batching convenience*, not a REST requirement. The planner can collapse the per-species batching down to "one POST per ≤1000 Ensembl-eligible IDs" since IDs are unique. This simplifies the resolver. Document both options in the plan.

### Pattern 9: Group-Mean Correlation viewer (R4/D-12)

**What:** `createGroupMeanCorrelation(df)` factory parallel to `createVolcanoPlot`. Mean columns are derived per-row; Pearson/Spearman computed over filtered rows; diagonal reference at y=x; magenta/cyan/gray by significance.

**Key details:**

- **Derived columns:** `Numerator Mean`, `Denominator Mean`. Use `ensureFreshFloat(df, name)` pattern from `viewers/qc-computations.ts:28-32` to guard re-runs. Set semType so detectors.js can mirror.
- **Inline Pearson:** `r = Σ((xi-x̄)(yi-ȳ)) / sqrt(Σ(xi-x̄)² · Σ(yi-ȳ)²)`. Use `getRawData()` for both columns. Skip rows where either is `DG.FLOAT_NULL` or `NaN`.
- **Inline Spearman:** rank both arrays (tie via fractional ranking), Pearson on ranks.
- **Diagonal:** `df.meta.formulaLines.addLine({formula: '${Numerator Mean} = ${Denominator Mean}', color: '#888888', width: 1})`. Same idiom as `volcano.ts:184-188`.

### Anti-Patterns to Avoid

- **Re-binding the scatter color column on every counter recompute.** The counter overlay reads `colorColumnName` and counts categories; do NOT call `sp.setOptions({color: ...})` from inside the recompute — that triggers a re-render and a chain of property-change events. Read-only access.
- **Per-row `col.set()` for the new Display Name / Source ID columns.** Use `col.init((i) => arr[i])` bulk path. Memory: `feedback_dg_column_bulk_init` measured 137–255× speedups on the volcano/DE writeback paths.
- **Hardcoded `df.col('Gene names')`.** Use `findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol'])`. The resolver writes to a separate `Display Name` column — it does NOT mutate `Gene name`.
- **Forking `createVolcanoPlot`.** D-04/D-05/D-06 extend it in place. The Phase-13 architecture intentionally chose ONE volcano with property toggles. Forking is a regression.
- **Treating `columnNames` as the Filters scoping API.** Phase 13 documented `columnNames` silently extending. The typed `filters: [{type, column}]` array is the actual scoping mechanism. See Pattern 1.
- **Using `df.filter.setAll(false)` for D-05 search-match.** The contract is highlight-in-place. Use `df.selection`, leave `df.filter` alone.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Cross-session per-user key/value cache | localStorage / IndexedDB / cookies | `grok.dapi.userDataStorage.get/put` with `__schema_v` key | Phase 13 already proved the pattern; analyst sessions span devices; userDataStorage handles auth + per-user partitioning |
| CORS-safe external REST | `fetch()` / `axios` | `grok.dapi.fetchProxy(url, opts)` | Datagrok platform handles auth, retries, CORS; raw fetch fails behind enterprise proxies |
| Filters viewer column scoping | Custom textarea + manual highlight | `DG.Viewer.filters(df, {filters: [...]})` with typed filter specs | D-05 LOCKED; platform Filters viewer gives free-text search + categorical multi-select + selection-not-filter wiring |
| Top-N point labels | Custom canvas overlay or manual SVG | `sp.props.labelColumnNames` + `showLabelsFor: 'Selected'` + driving `df.selection` | Scatter has a native label layout engine with collision avoidance |
| Progress UX during long fetch | Custom progress bar | `DG.TaskBarProgressIndicator.create(label)` (already used at `package.ts:337`) | Platform top-bar progress is the established pattern; works with the cache hit/miss UX in place |
| Bounded-concurrency worker pool | `Promise.all` over chunks | `runWithConcurrency(items, limit, work)` already in `subcellular-location.ts:219` | Already in tree and tested; reuse not duplicate |
| Per-row Mean computation in correlation viewer | Loop with bridge hops | `ensureFreshFloat(df, name)` then `col.init((i) => arr[i])` after bulk-read source columns via `getRawData()` | 137–255× faster (memory feedback_dg_column_bulk_init) |
| RxJS subscription cleanup | Manual `.unsubscribe()` scattered | Pattern in `enrichment-viewers.ts:9` (`activeSubscriptions` array, dispose on re-open) | Established convention |

**Key insight:** Phase 14 is heavy on capability extension and light on novel mechanism. Every new capability has a Phase-9-or-13 in-tree precedent or a platform primitive. The planner's main job is wiring, not invention.

## Runtime State Inventory

Phase 14 is an enhancement phase — not a rename, refactor, or migration. **Skip this section.**

For the new `userDataStorage` STORE_GENE_LABELS entry (R1/D-09): this is forward-additive. No existing user-data-storage records need migrating. Existing users will see one Ensembl batch fetch on their first re-import after the upgrade; subsequent imports hit the warm cache.

## Common Pitfalls

### Pitfall 1: Filters viewer `columnNames` is not a strict allowlist (G4 root cause)

**What goes wrong:** Phase 13 used `DG.Viewer.filters(df, {columnNames: ['Comparison (group1/group2)']})` and observed the `Flags` boolean column appearing anyway.

**Why it happens:** Documented platform behavior — `IFiltersSettings.showBoolCombinedFilter` defaults to true; multiple boolean columns auto-combine into one filter regardless of `columnNames`. Spectronaut Candidates files carry a `Flags` (Valid/significant) column that gets folded into the combined-bool filter.

**How to avoid:** Use the typed `filters` property (`Array<{type, column}>`) AND set `showBoolCombinedFilter: false`. See Pattern 1. The planner should verify in a local test that the `filters` array reliably excludes columns NOT listed.

**Warning signs:** When opening the dialog, the Filters viewer shows MORE columns than the developer requested. Test: count `viewer.getOptions().look.filters?.length` vs the expected list length.

### Pitfall 2: `sp.getOptions()` returns a snapshot, not a live binding

**What goes wrong:** The G2 dialog reads `sp.getOptions().look.colorColumnName`, the user changes the color via the property panel WHILE the dialog is open, then OK uses the stale value.

**Why it happens:** `getOptions()` serializes at call time. It does not return an observable.

**How to avoid:** Read the options at OK time, not at dialog open time. Dialog open seeds the form; OK re-reads (or trusts the form inputs since the user just edited them). Specifically: snapshot at open, populate inputs, on OK use INPUT values (not re-read).

### Pitfall 3: Ensembl `/lookup/id` POST 1000-ID cap is the hard limit

**What goes wrong:** A single dataset with ~3000 ENSRNOG IDs sends one POST → 413 Payload Too Large OR truncated response.

**Why it happens:** Documented limit. [CITED: rest.ensembl.org/documentation/info/lookup_post: "Maximum POST Size: 1000 IDs per request"]

**How to avoid:** Chunk to 1000 IDs per batch. Reuse `chunk()` helper from `subcellular-location.ts:205`.

**Warning signs:** HTTP 4xx from `/lookup/id`. Datagrok's `fetchProxy` surfaces the body — log it.

### Pitfall 4: Ensembl rate limit (55,000 requests/hour, ~15 req/sec)

**What goes wrong:** A bounded-concurrency worker pool with limit=20 + a re-import after cache invalidation = brief burst over 15 req/sec = 429 Too Many Requests with `Retry-After: 40.0`.

**Why it happens:** [CITED: github.com/Ensembl/ensembl-rest/wiki/Rate-Limits] — 55k/hour, sub-second Retry-After header on 429.

**How to avoid:** Cap concurrency at the same `FETCH_CONCURRENCY = 6` value used by subcellular-location.ts. Respect `Retry-After` header on 429: parse the value, sleep, retry once. Reuse `runWithConcurrency` from `subcellular-location.ts:219`.

**Warning signs:** Sporadic null entries in the resolver result map. Check the response status — `console.warn` on 429s would reveal the pattern.

### Pitfall 5: Color column rename breaks the volcano/enrichment cross-link

**What goes wrong:** D-04 changes the direction column's category strings from `'up'/'down'/'not significant'` to `'Enriched in g1'/'Enriched in g2'/'Not significant'`. The `enrichment-viewers.ts:wireEnrichmentToVolcano` doesn't use those strings (it uses the Intersection column), so it's safe — but ANY test asserting on the literal `'up'`/`'down'` strings breaks.

**Why it happens:** `ensureDirectionColumn` in `volcano.ts:62` sets the value strings. Tests in `src/tests/` may assert against them.

**How to avoid:** Grep `'up'`, `'down'`, `'not significant'` across `src/tests/` BEFORE the migration; update assertions to use the new strings.

**Warning signs:** Phase-13 test suite goes red on direction-string assertions.

### Pitfall 6: Volcano counter overlay attached to a re-created viewer

**What goes wrong:** User runs `Visualize | Volcano Plot` twice; second invocation creates a NEW `DG.ScatterPlotViewer`, the old subscriptions still fire but point at a detached overlay.

**Why it happens:** No automatic teardown when the viewer is replaced.

**How to avoid:** Mirror the `enrichment-viewers.ts:9` `activeSubscriptions` cleanup pattern. Better: tie subscriptions to the `sp` viewer instance by stashing them in a `WeakMap<DG.ScatterPlotViewer, rxjs.Subscription[]>`, OR observe `sp.onDestroy` / `tv.viewerRemoved` events to dispose.

**Warning signs:** Multiple overlays stacking on the same volcano; memory leak warnings; counter values frozen on stale data.

### Pitfall 7: D-03 top-N labels collide with D-05 search-match selection

**What goes wrong:** Top-N labels SEED the selection bitset. User types `MYH7` in the Filters search box → planner subscribes to filter-change → wants to set selection on the match. If top-N selection logic re-runs first, it clobbers the search match.

**Why it happens:** Both mechanisms write to `df.selection`. Order of subscriptions matters.

**How to avoid:** Union the two sources — `setTopNLabels` accepts a `mode: 'replace' | 'union'` parameter; default `'replace'` on metric/color/filter change, `'union'` when re-applying after a search match. Alternative: track top-N indices separately and recompute selection as `topN ∪ searchMatches`.

**Warning signs:** Search match disappears when filter category is toggled.

### Pitfall 8: Smart pathway filter changes intersection-column gene names

**What goes wrong:** After smart filter drops rows, the cross-link in `enrichment-viewers.ts:wireEnrichmentToVolcano` highlights based on `Intersection` column values. If the filter is invoked BEFORE `buildEnrichmentDf`, the Intersection column already reflects only kept rows — fine. If invoked AFTER, the index `i` of `intersections[i]` no longer aligns with the kept row's `i`.

**Why it happens:** `buildEnrichmentDf` builds the per-row intersection string from `r.intersections` indexed by query position. Smart filter operates on `GostResult[]`, BEFORE the per-row strings are built. So intersection-by-index alignment is preserved.

**How to avoid:** Run `applySmartPathwayFilter` BEFORE `buildEnrichmentDf`. Order in `runEnrichmentPipeline`: gGOSt → applySmartPathwayFilter (gated on dialog checkbox) → buildEnrichmentDf.

### Pitfall 9: Display Name column changes break Volcano label binding silently

**What goes wrong:** `Display Name` semType (new) overlaps with `Gene name` semType for the label binding. `sp.props.labelColumnNames = ['Gene name']` (current) needs to become `['Display Name']`. If only some parsers populate `Display Name`, the volcano falls back to no labels for those files.

**Why it happens:** Five parsers, one resolver. Skipping a parser is a partial-coverage bug.

**How to avoid:** EVERY parser calls `resolveGeneLabels(df)`. The resolver MUST always create the `Display Name` column even if no predicted IDs are detected — populate it with the unmodified gene name as a no-op. Then the volcano can unconditionally bind to `'Display Name'`.

**Warning signs:** Volcano labels work for MaxQuant but disappear for FragPipe.

### Pitfall 10: G3 cache check hits a stale promise

**What goes wrong:** Per `feedback_dapi_fresh_instance_patch`, `grok.dapi.userDataStorage` returns a fresh instance per access. The `package.ts:325` code already accounts for this. Phase 14 must do the same.

**How to avoid:** Always invoke `grok.dapi.userDataStorage.get(...)` fresh. Do not cache the storage reference. Already correct in the codebase — just don't regress.

## Code Examples

### Ensembl /lookup/id POST request shape

```typescript
// Source: https://rest.ensembl.org/documentation/info/lookup_post
// Verified: rate limit 55000/hour ~ 15 req/sec; max 1000 IDs/request; 429 returns Retry-After (sub-second float)
async function lookupEnsemblBatch(ids: string[]): Promise<Map<string, EnsemblEntry>> {
  if (ids.length === 0) return new Map();
  if (ids.length > 1000) throw new Error('Ensembl /lookup/id supports ≤1000 IDs per request');
  const url = 'https://rest.ensembl.org/lookup/id';
  const resp = await grok.dapi.fetchProxy(url, {
    method: 'POST',
    headers: {'Content-Type': 'application/json', 'Accept': 'application/json'},
    body: JSON.stringify({ids}),  // expand:0 is the default; CK-omics passes expand:1 — we don't need expanded data
  });
  if (resp.status === 429) {
    const retryAfter = parseFloat(resp.headers.get('Retry-After') ?? '1');
    console.warn(`Ensembl rate-limited; retrying after ${retryAfter}s`);
    await new Promise((r) => setTimeout(r, retryAfter * 1000));
    return lookupEnsemblBatch(ids);  // single retry
  }
  if (!resp.ok) {
    console.warn(`Ensembl /lookup/id returned ${resp.status}`);
    return new Map();
  }
  const json = await resp.json();
  // Response shape: { "ENSG00000157764": {display_name, external_name, description, species, biotype, ...},
  //                   "ENSG00000999999": null  // unknown ID → null entry }
  const out = new Map<string, EnsemblEntry>();
  for (const [id, entry] of Object.entries(json)) {
    if (entry && typeof entry === 'object')
      out.set(id, entry as EnsemblEntry);
  }
  return out;
}

interface EnsemblEntry {
  display_name?: string;       // e.g., "BRAF"
  external_name?: string;      // e.g., "BRAF" — usually same as display_name for genes
  description?: string;        // e.g., "B-Raf proto-oncogene, serine/threonine kinase [Source:HGNC...]"
  species?: string;            // e.g., "homo_sapiens"
  biotype?: string;            // e.g., "protein_coding"
  object_type?: string;        // e.g., "Gene"
}
```

### Verifying the volcano color-string migration

```typescript
// Source: extend src/viewers/volcano.ts:62-86 — ensureDirectionColumn
// LOCKED: D-04 group names from getGroups(df), color hex from CONTEXT.md
import {getGroups} from '../analysis/experiment-setup';

export const DIRECTION_COLORS_BASE = {
  enrichedG1: 0xFFFF00FF,  // magenta
  enrichedG2: 0xFF00FFFF,  // cyan
  notSig:     0xFFAAAAAA,  // gray (matches existing volcano.ts:84 'not significant')
} as const;

export function ensureDirectionColumn(
  df: DG.DataFrame, fcThreshold: number, pThreshold: number, metric: MetricKind = 'adj.p-value',
): string {
  const fcCol = df.col('log2FC');
  if (!fcCol) throw new Error('log2FC column not found');
  const pCol = pickMetricColumn(df, metric);
  const groups = getGroups(df);
  const g1Label = groups ? `Enriched in ${groups.group1.name}` : 'Enriched in group1';
  const g2Label = groups ? `Enriched in ${groups.group2.name}` : 'Enriched in group2';
  const nsLabel = 'Not significant';

  const col = df.col(DIRECTION_COL) ?? df.columns.addNewString(DIRECTION_COL);
  col.init((i) => {
    if (fcCol.isNone(i) || pCol.isNone(i)) return nsLabel;
    const fc = fcCol.get(i) as number;
    const p = pCol.get(i) as number;
    if (p <= pThreshold && fc > fcThreshold) return g1Label;
    if (p <= pThreshold && fc < -fcThreshold) return g2Label;
    return nsLabel;
  });
  col.meta.colors.setCategorical({
    [g1Label]: DIRECTION_COLORS_BASE.enrichedG1,
    [g2Label]: DIRECTION_COLORS_BASE.enrichedG2,
    [nsLabel]: DIRECTION_COLORS_BASE.notSig,
  });
  return DIRECTION_COL;
}
```

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework | `@datagrok-libraries/test` ^1.1.0 (already a devDep) + `grok test` runner |
| Config file | none — discovered via webpack two-entry-point pattern |
| Quick run command | `grok test --test "<TestName>"` (single test by name prefix) |
| Full suite command | `grok test` (from package dir, requires running Datagrok instance) |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| R1 (parse-time resolution) | Display Name populated, * and † markers correct, predicted prefixes detected, LOC/RGD/AABR keep raw ID + † | unit | `grok test --test "geneLabelResolver"` | Wave 0 (NEW `src/tests/gene-label-resolver.ts`) |
| R1 (cache) | `__schema_v` bump invalidates cache; cache miss triggers Ensembl POST; cache hit short-circuits | unit (mocked fetchProxy) | `grok test --test "geneLabelResolverCache"` | Wave 0 |
| R1 (Ensembl POST shape) | One POST per ≤1000 IDs; 429 with Retry-After triggers single retry | unit (mocked fetchProxy) | `grok test --test "geneLabelEnsemblPost"` | Wave 0 |
| R1 (duplicate disambig) | Two identical Display Names with different Source IDs → disambiguated with `(Source ID)` suffix; warning toast count | unit | `grok test --test "geneLabelDuplicates"` | Wave 0 |
| R2 (counter math) | `df.filter` mask → category counts; `df.selection` change doesn't affect counts; metric toggle drives recompute | unit | `grok test --test "volcanoCounter"` | extend `src/tests/volcano.ts` |
| R3 (per-group bars) | Bars render for protein with valid group quants; empty state when groups unset; mean/SD match manual computation | unit | `grok test --test "uniprotPanelPerGroup"` | Wave 0 (new test file or extend) |
| R4 (Pearson/Spearman) | Known fixture (linear data → r=1, ρ=1; uncorrelated → r≈0); rank-tie handling for Spearman | unit | `grok test --test "groupMeanCorrelation"` | Wave 0 NEW `src/tests/group-mean-correlation.ts` |
| R4 (derived cols) | Numerator Mean / Denominator Mean populated; semTypes set; re-run replaces (ensureFreshFloat) | unit | `grok test --test "groupMeanColumns"` | Wave 0 |
| R5 (smart filter) | Generic GO:BP parent with specific child → parent dropped; cap at maxPerSource; other sources passthrough | unit | `grok test --test "smartPathwayFilter"` | Wave 0 |
| R5 (dialog default-on) | Checkbox default checked; uncheck → raw flow; stats banner triggers | unit | `grok test --test "enrichmentDialogSmartFilter"` | Wave 0 |
| G1 (color hex) | Direction column color map contains exact `0xFFFF00FF`/`0xFF00FFFF`/`0xFFAAAAAA` after migration | unit | `grok test --test "volcanoDirectionColors"` | extend `src/tests/volcano.ts` |
| G1 (legend strings) | Direction column categories include `"Enriched in ${g1}"`/`"Enriched in ${g2}"`/`"Not significant"` derived from `getGroups(df)` | unit | `grok test --test "volcanoDirectionStrings"` | extend `src/tests/volcano.ts` |
| G1 (top-N labels) | Top-15 rows by ascending metric AND tiebreaker |log2FC| desc match | unit | `grok test --test "volcanoTopN"` | extend `src/tests/volcano.ts` |
| G1 (visual parity) | Side-by-side comparison against `~/Downloads/ck/DMD_vs_WT/volcano_plots/` | manual-only | n/a — human UAT | n/a |
| G2 (preload) | After volcano metric toggle + dialog re-open, the metric input shows the active state | manual-only | n/a — human UAT (DOM read at dialog-show time difficult to assert automatically) | n/a |
| G3 (progress) | TaskBarProgressIndicator created before fetch begins; closed in `finally` even on error | unit | `grok test --test "locationProgressLifecycle"` | extend `src/tests/spectronaut-parser.ts` or `volcano.ts` |
| G3 (perf) | Warm-cache re-toggle Color→Location completes <1s on a 5000-protein fixture | perf-bench (manual) | n/a — wall-clock measure | n/a |
| G4 (Flags excluded) | After `dockComparisonFilterIfMultiContrast(tv, dfWithFlags)`, the docked Filters viewer's `getOptions().look.filters` does NOT reference Flags | unit | `grok test --test "filtersScopingNoFlags"` | Wave 0 (NEW test or extend `src/tests/spectronaut-parser.ts`) |

### Sampling Rate

- **Per task commit:** `grok test --test "<area>"` (e.g., `--test "volcano"`, `--test "geneLabel"`)
- **Per wave merge:** `grok test` (full suite)
- **Phase gate:** Full suite green BEFORE `/gsd:verify-work`; manual-only items (G1 visual parity, G2 preload, G3 perf wall-clock) handled in `14-HUMAN-UAT.md`

### Wave 0 Gaps

- [ ] `src/tests/gene-label-resolver.ts` — NEW; covers R1 contract (port + cache + Ensembl POST mock + duplicate disambig)
- [ ] `src/tests/group-mean-correlation.ts` — NEW; covers R4 (Pearson, Spearman, derived columns)
- [ ] Mock helper for `grok.dapi.fetchProxy` — Wave 0 may need a fixture pattern; check `src/tests/spectronaut-parser.ts` for existing fetch-mock style first
- [ ] Extend `src/tests/volcano.ts` for D-04 color migration, D-03 top-N selection logic, D-06 counter math
- [ ] Extend `src/tests/enrichment.ts` for R5 smart filter cap+drop assertions
- [ ] Fixture data: small (~20-row) synthetic `Display Name` test corpus covering each species prefix + LOC/RGD/AABR cases; fixture g:Profiler response with both generic GO:BP parents and specific children

### Note: existing volcano tests will fail until D-04 migration completes

Tests in `src/tests/volcano.ts` (and possibly `enrichment.ts`) that assert literal `'up'`/`'down'`/`'not significant'` direction strings will go red the moment `ensureDirectionColumn` flips to the new `"Enriched in g1"` form. Plan: update the assertions as part of the SAME task that migrates the color strings — do not let red tests linger across task boundaries.

## Security Domain

> Phase 14 `security_enforcement` status not explicit in config.json. Project-level config has no `security_enforcement: false` toggle. Treating as enabled per the default.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | No auth surface — package only calls public REST (Ensembl, UniProt, g:Profiler) and platform's authenticated `grok.dapi.userDataStorage` |
| V3 Session Management | no | Platform owns session; package contributes nothing |
| V4 Access Control | no | No new access boundaries |
| V5 Input Validation | yes | Ensembl response untrusted; planner must validate field types before write; description sanitization (CK-omics regex strips, take-first-sentence) prevents pathological input |
| V6 Cryptography | no | Plain-text REST; no PII; gene IDs and locations are public reference data |
| V11 Business Logic | yes | Smart filter changes downstream analysis results — must be flagged via banner to analyst (CONTEXT.md explicitly requires) |
| V13 API & Web Service | yes | External REST calls — must use `grok.dapi.fetchProxy()` not raw `fetch()` (CORS + auth handled) |
| V14 Configuration | no | No new config |

### Known Threat Patterns for {Datagrok TS package + external REST}

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Untrusted Ensembl response with malformed strings | Tampering | Type-guard every field read; ignore non-string `display_name`/`description`; HTML-escape if rendered as innerHTML (we don't — `ui.divText` and SVG `textContent` are safe) |
| Cache poisoning via shared userDataStorage | Tampering | `__schema_v` invalidation; cache entries are scoped per-user by platform; planner cannot bypass this |
| Excessive Ensembl request rate (DoS to upstream) | Denial of Service | Bounded concurrency at 6 (Phase 13 pattern); respect `Retry-After` on 429 |
| Silent transformation of enrichment results (R5) | User trust / Repudiation | Banner above grid when `kept < total`; dialog checkbox makes it opt-out; df tag `proteomics.enrichment_smart_filtered` persists state |
| Leaking analyst protein lists to third parties | Information Disclosure | Ensembl POST sends ONLY ID strings (no analyst quant data); g:Profiler already receives gene lists (Phase 13 carry-forward, accepted) |

No new external endpoints beyond Ensembl `/lookup/id`. UniProt and g:Profiler are Phase-13 inheritances. Threat surface change is minimal.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Hardcoded gene-name labels with no Ensembl resolution | Eager Ensembl resolution at parse time with `__schema_v`-keyed cache | Phase 14 (this phase) | Closes BP DMD/WT analyst readability gap; analyst never sees raw ENSRNOG IDs |
| Red/blue/gray volcano coloring | Magenta/cyan/gray with group-name legend categories | Phase 14 D-04 | CK-omics visual parity; LOCKED color contract |
| `columnNames` allowlist for Filters viewer | Typed `filters: Array<{type, column}>` array | Phase 14 G4 | Robust scoping that drops `Flags` reliably |
| Per-row `col.set()` writes | Bulk `col.init((i) => arr[i])` + `getRawData()` reads | Phase 11 carry-forward | 137–255× speed improvement on column writes (project memory) |
| Per-tooltip raw `fetch()` | `grok.dapi.fetchProxy()` with `__schema_v` cache | Phase 13 D-01/D-02 | CORS + cross-session caching; mirror this exact pattern for R1 |

**Deprecated / outdated within Phase 14 scope:**
- `'up'` / `'down'` / `'not significant'` direction strings — replaced by group-name-aware strings under D-04.
- `sp.props.showLabelsFor = 'MouseOverRow'` as the sole label mode — D-03 adds default top-N selection-driven labels (MouseOverRow stays for hover).
- `Gene name` column as the volcano label source — D-08 makes `Display Name` the canonical label everywhere.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Setting `IFiltersSettings.filters` with explicit per-column objects + `showBoolCombinedFilter: false` is sufficient to drop `Flags` from the docked Filters viewer | Pattern 1 / G4 | Phase 13 observation confirms `columnNames` alone is insufficient; the typed `filters` array is the documented [CITED: datagrok.ai/help/visualize/viewers/filters] but the planner must verify locally. Fallback: post-create `viewer.setOptions(...)` with both fields. |
| A2 | Ensembl `/lookup/id` returns the IDs that ARE found as nested objects keyed by ID, and IDs NOT found as `null` or omitted | Code Examples / R1 verbatim port | [CITED: rest.ensembl.org docs say "no data" for unknown IDs] but exact null-vs-omit shape needs a smoke test. Planner should hit the endpoint with one good + one bad ID and capture the response. |
| A3 | Reusing `runWithConcurrency` from `subcellular-location.ts:219` works without modification for Ensembl batched POSTs | Pattern 4 | LOW — the function is generic over `T` items and async work. |
| A4 | The CK-omics `apply_smart_pathway_filtering` line 4731 `.head(maxPerSource)` over combined `other_data` (not per-source) is intentional, NOT a bug | Pattern 7 / R5 | Reading the function: `other_data` is everything NOT GO:BP. CK-omics applies a single combined cap. If the planner instead caps per-source (one per KEGG, REAC, WP, GO:MF, GO:CC), the result diverges from CK-omics. D-13 says "verbatim port" — mirror lines 4730-4733 literally. If client later wants per-source, that's a follow-up phase. |
| A5 | `df.selection.set(...)` does NOT trigger `df.filter.onChanged` (D-05 highlight-not-hide invariant) | Pattern 3 / D-05 | The observables are independent (`onSelectionChanged` vs `onFilterChanged`) per js-api source. The Filters viewer's own search-box behavior is the open question — does typing into a free-text filter mutate `df.filter` or just visually narrow? Local test required. |
| A6 | Two-source-of-truth for G2 (sp.getOptions for color/y, df.tag for metric) is sufficient — no need to track all dialog inputs in a closure cache | Pattern 5 / G2 | If a future Volcano Options input has no viewer-side or tag representation, this falls through. For Phase 14 there are only two inputs (metric + colorDim), both representable. |
| A7 | The /lookup/id endpoint does not require `species` parameter for cross-species lookup; IDs are globally unique within Ensembl | Pattern 8 footnote | [CITED: rest.ensembl.org/documentation/info/lookup_post — species is OPTIONAL]. CK-omics omits it (line 779-787). Safe. |
| A8 | LOC/RGD/AABR predicted IDs are NOT Ensembl IDs and CK-omics intentionally never sends them to Ensembl — they keep the raw ID + `†` marker only | Pattern 8 / R1 contract | [VERIFIED: CKomics_tool2.py:933 explicitly filters `[pid for pid in all_predicted_ids if pid.startswith(('ENS', 'MGP_'))]`]. Critical for matching CK-omics behavior. |
| A9 | Ensembl rate limit of 55,000 req/hour ≈ 15 req/sec accommodates any single Proteomics import — a 50k-protein dataset chunks to 50 batched POSTs, well under the per-hour budget | Pitfall 4 / R1 | Documented limit. The bigger risk is concurrent users on the same Datagrok instance sharing the per-IP quota. Mitigation: cap at 6 concurrent (matches Phase-13 UniProt cap). |
| A10 | `df.col('Subcellular Location').categories` returns the live category list from the Datagrok DataFrame, suitable for driving the counter overlay's by-location breakdown | Pattern 2 / D-06 location-mode counter | Datagrok categorical columns expose `.categories`. Planner should verify the API name (`.categories` vs `.uniqueValues` etc.) — js-api docs needed for exact form. |

**This is a long assumptions list because Phase 14 spans 9 different surfaces.** The MUST-have trio (G1, R1, R2) hinges on A1, A2, A4, A5. Resolve those four via smoke tests before locking the plan.

## Open Questions (RESOLVED)

1. **Does the Filters viewer's free-text search-box mutate `df.filter` or surface a separate "matched rows" event?**
   - What we know: D-05 requires highlight-not-hide. `df.selection` is the highlight mechanism.
   - What's unclear: whether typing in the free-text box hides non-matching rows via `df.filter` (which would break D-05's invariant) or just emits a "match" signal we can intercept.
   - RESOLVED: 30-minute spike at plan-start. Dock a free-text filter on a 100-row DataFrame, type a substring, inspect `df.filter.trueCount`. If it changes, planner must override the default behavior (e.g., subscribe to the change, restore `df.filter`, drive selection from the match) — or use a custom filter type.

2. **Where do the "specific child terms" come from in `apply_smart_pathway_filtering` for non-CK-omics datasets?**
   - What we know: CK-omics ships a fixed list `['actin', 'vesicle', 'endocytosis', 'cytoskeleton']`. This is hardcoded to muscle/myo biology (DMD/WT context).
   - What's unclear: whether the list generalizes to other datasets the package will encounter.
   - RESOLVED: ship the verbatim list per D-13 (it's a locked client contract). Document in a code comment that the list is dataset-shaped. Future-phase follow-up: expose the specific-child terms as a configurable parameter (already in 14-CONTEXT.md "Deferred Ideas: Configurable per-source pathway cap").

3. **For G3 perf, is the cache hit path measurable as sub-second OR does the bulk-init step itself take time?**
   - What we know: `ensureLocationColumn` in `volcano.ts:120-168` already short-circuits when the column exists with the right hash (`LOCATION_HASH_TAG`). The cold path is the slow one.
   - What's unclear: whether the warm path (cache hit, but column needs re-creating) measures fast enough.
   - RESOLVED: planner adds a perf-bench test that runs `ensureLocationColumn` against a 5000-protein fixture twice — first cold, then warm — and asserts the warm path is <1s. If not, add memoization on the `accForRow` array within the function.

4. **Should the `Display Name` column replace `Gene name` in the grid layout by default?**
   - What we know: CONTEXT.md Claude's Discretion explicitly defers this to the planner.
   - What's unclear: client preference.
   - RESOLVED: keep BOTH columns. `Display Name` becomes the visible-by-default label; `Gene name` is hidden in default grid layout (via `df.col('Gene name').setTag('.row-height', null)` or a layout) but accessible via column chooser. The Volcano binds to `Display Name`; the grid shows `Display Name` as the primary; analysts who want raw IDs can show `Gene name` and `Source ID`.

5. **Banner rendering location for R5 — above the enrichment grid or inside the existing enrichment dot/bar chart layout?**
   - What we know: UI-SPEC.md "Banner rendering" says "docked ABOVE the enrichment grid in the same table view".
   - What's unclear: whether docking lands in the existing table view (which already has dot/bar charts via `openEnrichmentVisualization`) or a separate widget panel.
   - RESOLVED: `tv.dockManager.dock(ui.divText(banner), DG.DOCK_TYPE.TOP, null, 'Smart Filter Info', 0.1)` on the enrichment table view, BEFORE the dot/bar dock calls in `enrichment-viewers.ts:189-196`. Persists for the life of the view.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Ensembl REST API (`https://rest.ensembl.org`) | R1/D-09 | live external service | n/a | Graceful degrade: render raw IDs + single warning toast; do not block import (per CONTEXT.md and UI-SPEC.md error states) |
| UniProt REST API (Phase 13 carry-forward) | R3 panel | live external service | n/a | already in place; this phase doesn't change the dependency |
| g:Profiler API (Phase 8 carry-forward) | R5 wraps existing call | live external service | n/a | already in place |
| Datagrok platform `userDataStorage` | R1 cache | platform-bundled | n/a (host platform owns version) | none — required; if absent, R1 falls back to in-session-only resolution (slower on re-import, no data lost) |
| `grok.dapi.fetchProxy` | All external REST | platform-bundled | n/a | none |
| Node + npm + grok CLI (build) | dev/test | dev environment | host-dependent | `grok publish localhost` for fast iteration (skip webpack — won't work here since src/* changes) |
| R compute environment | none in Phase 14 | n/a — Phase 14 does NOT call any R scripts | n/a | n/a |

**Missing dependencies with no fallback:**
- None. Phase 14 graceful-degrades on every external service.

**Missing dependencies with fallback:**
- Ensembl REST unavailable → analyst sees raw IDs + one toast; phase still ships per the locked R1 fallback chain (`Ensembl name → Ensembl description → raw ID`).

## Sources

### Primary (HIGH confidence)
- `js-api/src/viewer.ts:158` — `Viewer.getOptions()` returns `{id, type, look: {...}}`
- `js-api/src/viewer.ts:235` — `Viewer.filters(t, options?)` static factory; takes `Partial<IFiltersSettings>`
- `js-api/src/interfaces/d4.ts:1674-1710` — `IFiltersSettings` shape (columnNames, filters, showBoolCombinedFilter, showSearchBox)
- `js-api/src/interfaces/d4.ts:322-331` — `RowSet` enum (Selected, Filtered, MouseOverRow, etc.)
- `js-api/src/interfaces/d4.ts:552-555` — `labelColumnNames` + `showLabelsFor`
- `js-api/src/dataframe/data-frame.js:387,389` — `onSelectionChanged`/`onFilterChanged` observables
- `libraries/statistics/src/correlation-coefficient.ts:68` — only `kendallsTau` exported (Pearson/Spearman absent)
- `~/Downloads/ck/CKomics_tool2.py:756-855` — `get_ensembl_annotations` (R1 contract source)
- `~/Downloads/ck/CKomics_tool2.py:895-1059` — `improve_gene_labels_with_ensrnog_marking` (R1 contract source)
- `~/Downloads/ck/CKomics_tool2.py:1062-1100` — `extract_readable_description` (R1 cleanup contract)
- `~/Downloads/ck/CKomics_tool2.py:4685-4736` — `apply_smart_pathway_filtering` (R5 contract source)
- `packages/Proteomics/src/analysis/subcellular-location.ts:165-393` — Phase 13 cache + worker pool pattern (mirror for R1)
- `packages/Proteomics/src/viewers/volcano.ts:1-250` — `createVolcanoPlot` / `recomputeVolcano` / `ensureDirectionColumn` (extend in place for G1/D-03/D-04/D-06)
- `packages/Proteomics/src/package.ts:117-131` — `dockComparisonFilterIfMultiContrast` (extend for D-05/G4)
- `packages/Proteomics/src/panels/uniprot-panel.ts:109-176` — `renderUniProtWidget` (extend for R3/D-11)
- `packages/Proteomics/src/analysis/enrichment.ts:90-130, 250-367, 371-467` — `gGOSt`, `runEnrichmentPipeline`, `showEnrichmentDialog` (extend for R5/D-13/D-14)
- `packages/Proteomics/src/viewers/enrichment-viewers.ts:9` — `activeSubscriptions` cleanup pattern
- `packages/Proteomics/src/viewers/qc-computations.ts:28-32` — `ensureFreshFloat` re-run-safe column pattern
- `packages/Proteomics/src/parsers/spectronaut-candidates-parser.ts:1-269` — full Candidates parser contract (Phase 13)
- `rest.ensembl.org/documentation/info/lookup_post` (WebFetch verified) — Ensembl /lookup/id POST contract: max 1000 IDs, unknown IDs return "no data", optional species/expand/db_type/object_type/format params
- `github.com/Ensembl/ensembl-rest/wiki/Rate-Limits` (WebFetch verified) — 55,000 req/hour ≈ 15 req/sec; X-RateLimit-{Limit,Period,Remaining,Reset} response headers; 429 includes `Retry-After` in sub-second float
- `github.com/Ensembl/ensembl-rest/wiki/POST-Requests` (WebFetch verified) — POST batching guidance, Content-Type: application/json required
- `datagrok.ai/help/visualize/viewers/filters` (WebFetch verified) — filter type names (`categorical`, `free-text`, `histogram`, `text`, `expression`, `hierarchical`, `combined-boolean`); per-column object shape `{type, column}`

### Secondary (MEDIUM confidence)
- `~/Downloads/ck/DMD_vs_WT/volcano_plots/Subcellular_Location_Classification_README.txt` — locked Phase 13 contract; relevant to G3 (color path the speed-fix targets)
- Project memory `feedback_dapi_fresh_instance_patch` — userDataStorage returns fresh instance per access; do not patch the instance
- Project memory `feedback_dg_column_bulk_init` — bulk `col.init` over per-row `col.set` is 137–255× faster on volcano/DE
- Project memory `feedback_dg_column_init_null_sentinel` — `col.init(()=>null)` on numeric col leaves FLOAT_NULL sentinel unsynced; use `isNone()`/`set(i, null)` for null writes
- Project memory `feedback_proteomics_execute_phase_inline` — execute-phase subagents lack Bash; commit plans before execution, run inline sequential

### Tertiary (LOW confidence)
- None — every claim in this document is backed by an in-tree file or an authoritative external doc.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — every API in the recommended stack is verified in `js-api` source code
- Architecture: HIGH — every pattern has an in-tree precedent (Phase 13 cache, Phase 9 enrichment cross-link, Phase 11 bulk-init)
- Pitfalls: HIGH (Pitfalls 1, 3, 4, 5, 6, 9 backed by code or docs) / MEDIUM (Pitfalls 7, 8, 10 are forward-looking risks based on the wiring we're about to do)
- R1 CK-omics port: HIGH — read the source function in full; documented every contract element line-by-line in Pattern 8
- R5 CK-omics port: HIGH — read the source function in full; verified `.head(maxPerSource)` is a combined cap (line 4731) not per-source

**Research date:** 2026-05-28
**Valid until:** 2026-06-27 (30 days — Datagrok platform, CK-omics contracts, and Ensembl REST are all stable surfaces; Phase 14 should ship inside this window per ROADMAP cadence)
