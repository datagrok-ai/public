# Phase 13: CK-omics Volcano and Enrichment Parity - Context

**Gathered:** 2026-05-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Bring the Datagrok Proteomics volcano and enrichment outputs to feature parity
with the client's CK-omics tool deliverables for the BP DMD/WT engagement, so an
analyst gets the same volcano figure and the same enrichment scientific result
without leaving Datagrok. Six requirements (R1–R6, ROADMAP.md): subcellular-location
coloring on the volcano, up/down-split enrichment, Spectronaut Candidates contrast
normalization, multi-contrast selection, WikiPathways enrichment source, and a
Q-value/P-value significance toggle. Scope is fixed by ROADMAP.md — this
discussion clarifies HOW, not WHETHER.

</domain>

<decisions>
## Implementation Decisions

### Subcellular-location fetch & cache (R1)
- **D-01:** Fetch via the UniProt **stream** endpoint through `grok.dapi.fetchProxy()`
  — `/uniprotkb/stream?query=accession:(A OR B …)&fields=accession,cc_subcellular_location,go_c`,
  chunked ~100 accessions per request. Replaces the per-row raw `fetch()` shape;
  CORS-safe (fixes the documented raw-fetch concern in CONCERNS.md/INTEGRATIONS.md).
- **D-02:** **Persistent cross-session cache** keyed by accession (IndexedDB or
  `grok.userSettings` — Claude's discretion on mechanism). Lives in a shared
  subcellular-location module so the existing `uniprot-panel.ts` can be refactored
  to reuse it. This **fully closes** the cache-uniprot todo (see Folded Todos).
- **D-03:** Keep CK-omics' **reviewed-entry-by-gene-name fallback**: when an
  unreviewed entry has an empty subcellular field, look up the reviewed entry by
  gene name (one extra batched query for the misses) and use its location.
- **D-04:** Classify into the **11 categories + Unknown** using the keyword map,
  ordering rule (UniProt order-of-importance; GO-CC fallback), and **exact hex
  palette** from the CK-omics `Subcellular_Location_Classification_README.txt`.
  This README is a **locked client contract** — port verbatim, do not re-derive.

### Volcano UI: location colour + Q/P toggle (R1, R6)
- **D-05:** **One volcano**, not two. Add a `Subcellular Location` categorical
  column with the locked 11-colour map; the volcano switches its colour between
  significance (up/down/NS) and Subcellular Location via the viewer colour
  property plus a menu/dialog switch. (CK-omics ships two static HTMLs only
  because it has no interactive viewer.)
- **D-06:** R6 Q-value vs P-value is a **live viewer-property toggle** — switched
  from the property panel after creation. It recomputes the Y axis **and** the
  up/down/NS classification **and** the threshold lines together, kept in sync on
  every toggle so dots and colouring never disagree. Default = `adj.p-value`
  (= the Q-value the client volcano and the concordance work use).

### Spectronaut Candidates contrast (R3, R4)
- **D-07:** R4 — **import all comparisons** (no modal numerator/denominator
  picker). Retain the `Comparison` column and expose it through a **native
  Datagrok Filter viewer** docked with the volcano so the analyst narrows to one
  contrast interactively. Single-comparison files need no filter. **Tradeoff
  explicitly accepted:** no explicit "view the mirror of a declared contrast"
  control — filter-widget simplicity over the CK-omics modal.
- **D-08:** R3 — canonical convention: **positive log2FC = group1-enriched**,
  where group1 is the declared numerator in each `A / B` comparison string. The
  parser normalizes each row's sign to be consistent with its own declared
  `group1/group2` label and swaps/relabels the AVG Group Quantity Numerator/
  Denominator columns to match (CK-omics `create_subset_data`). The declared
  label IS the orientation — no per-run orientation prompt.
- **D-09:** R3 **also unparks the report-import DE direction fix** in
  `src/analysis/differential-expression.ts`: honor the intended/declared contrast
  instead of the alphabetical-condition-order default, and retire the manual
  Comparison-dropdown workaround. This phase is the explicit prompt that
  overrides the "do not implement unprompted" note in memory
  `project_proteomics_spectronaut_de_direction_default`.

### Up/down enrichment presentation (R2)
- **D-10:** Two g:Profiler queries (up-regulated genes, down-regulated genes)
  merged into **one enrichment DataFrame with a `Direction` column** (Up/Down).
  Reuses the Phase-9 dot-plot/bar-chart and volcano cross-link unchanged.
  Background stays **all detected proteins** (existing custom-background
  behaviour) for both directional queries.
- **D-11:** Present the two directions as **split Up/Down dot + bar charts
  side-by-side** (mirrors CK-omics directional deliverables; source filter still
  applies to both).

### WikiPathways (R5)
- **D-12:** Add `'WP'` to the g:Profiler source list, **default-on**. The
  `/api/gost/profile/` endpoint already supports WikiPathways
  (INTEGRATIONS.md) — minimal change, matches CK-omics' default-on WP.

### Claude's Discretion
- Accession chunk size; cache storage mechanism (IndexedDB vs `grok.userSettings`)
  and invalidation policy; progress-indicator cadence for the location fetch.
- Exact wording/placement of the colour-toggle switch and how the live Q/P
  toggle is wired into the `ScatterPlotViewer` property panel.
- Dock arrangement/sizing for the side-by-side Up/Down enrichment charts and for
  the `Comparison` Filter viewer relative to the volcano.
- Whether the directional gene lists are also surfaced (CK-omics writes
  `*_up_genes_all*.txt` / `*_down_genes_all*.txt`) — not required, planner's call.

### Folded Todos
- **`2026-03-03-cache-uniprot-protein-content-to-avoid-repeated-api-calls.md`**
  — original problem: `uniprot-panel.ts` does one uncached raw `fetch()` per row
  navigation. Folded into **D-02**: the persistent cross-session cache in a
  shared module that the panel is refactored to reuse closes this fully (not just
  for the panel — also for R1's bulk location fetch). Mark the todo done when
  Phase 13 ships D-02.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope
- `.planning/ROADMAP.md` § "Phase 13: CK-omics Volcano and Enrichment Parity" — the six locked requirements R1–R6 and their one-line specs

### Client reference tool (parity source — the deliverable being matched)
- `~/Downloads/ck/CKomics_tool2.py` — the client's CK-omics tool. Key functions:
  `get_subcellular_locations_uniprot` (~:1148), `parse_subcellular_location`
  (:1384), `get_location_colors` (:1608), `create_subset_data` (:1625, the R3
  flip reference), `run_gprofiler_analysis` (:4450, the R2 up/down split),
  `get_selected_databases` (:3850, R5 sources)
- `~/Downloads/ck/DMD_vs_WT/volcano_plots/Subcellular_Location_Classification_README.txt`
  — **LOCKED client contract** for D-04: the 11-category keyword map, ordering
  rule, and exact hex palette. Port verbatim.
- `~/Downloads/ck/Spectronaut-Concordance-Note-2026-05-15.md` — why default
  significance is adj.p (= Storey-q the client uses) at Q<0.05 / |log2FC|>0.58

### Existing code (the integration surface)
- `packages/Proteomics/src/panels/uniprot-panel.ts` — existing UniProt fetcher
  (`fetchUniProtData` ~:74-83); refactor to share the D-02 cache + D-01 transport
- `packages/Proteomics/src/viewers/volcano.ts` — significance hardwired to adj.p
  (~:42); D-05 colour toggle and D-06 live Q/P toggle land here
- `packages/Proteomics/src/analysis/enrichment.ts` — `gConvert` (~:77), `gGOSt`
  (~:99), custom background; D-10/D-11 up/down split and D-12 WP source
- `packages/Proteomics/src/viewers/enrichment-viewers.ts` — Phase-9 dot/bar +
  `onCurrentRowChanged` volcano cross-link to extend for the Direction dimension
- `packages/Proteomics/src/parsers/spectronaut-candidates-parser.ts` — keeps
  `Comparison` verbatim (~:85-89); D-07/D-08 contrast normalization land here
- `packages/Proteomics/src/analysis/differential-expression.ts` — D-09 report-
  import direction fix (alphabetical default → declared/intended contrast)
- `packages/Proteomics/src/utils/proteomics-types.ts` — add
  `SEMTYPE.SUBCELLULAR_LOCATION`; mirror in `packages/Proteomics/detectors.js`

### Codebase maps & prior decisions
- `.planning/codebase/INTEGRATIONS.md` — UniProt/g:Profiler transport, raw-fetch
  CORS concern, g:Profiler WP support
- `.planning/codebase/STRUCTURE.md` — where parsers/viewers/analysis/tests live
- `.planning/codebase/CONVENTIONS.md` — package coding conventions
- `.planning/milestones/v1.2-phases/09-enrichment-visualization-volcano-integration/09-CONTEXT.md`
  — Phase-9 enrichment-viewer + volcano cross-link pattern R2 must extend
- `~/.claude/projects/-Users-edjaeger-datagrok-src-public/memory/project_proteomics_spectronaut_de_direction_default.md`
  — the parked report-DE direction defect that D-09 now unparks

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `uniprot-panel.ts` `fetchUniProtData` — extend/refactor for D-01 batched
  transport + D-02 shared cache; the UniProt JSON `comments[].subcellularLocations`
  + `go` (cellular-component) fields back D-03/D-04
- `enrichment.ts` custom-background + `gGOSt` — call twice (up/down) for D-10;
  add `'WP'` to sources for D-12
- `enrichment-viewers.ts` Phase-9 dot/bar + `onCurrentRowChanged` cross-link —
  extend with a `Direction` column rather than rebuild (D-10/D-11)
- `volcano.ts` significance/threshold-line logic — refactor so the metric is a
  parameter recomputed on toggle (D-06)
- `findColumn`/`findProteomicsColumns` (`utils/column-detection.ts`),
  `SEMTYPE.*` (`utils/proteomics-types.ts`) — per package conventions

### Established Patterns
- `parseXText(text) → DG.DataFrame` parser purity (Phase 12): the Candidates
  parser stays pure; the Filter-viewer wiring (D-07) is orchestrated by the
  `package.ts` import handler / viewer layer, not the parser
- Bulk `Column.init`/`getRawData` not per-row `col.set` (memory
  `feedback_dg_column_bulk_init`) for the location column and the R3 sign flip
- New `SEMTYPE.*` requires a mirrored `detectors.js` entry (package contract)
- Native Datagrok viewer property panel for live options (Phase 9 precedent for
  D-05/D-06); native Filter viewer for D-07

### Integration Points
- Volcano colour + Q/P toggle: `volcano.ts` + its menu handler in `package.ts`
- Candidates contrast: `spectronaut-candidates-parser.ts` (sign normalization) +
  import handler in `package.ts` (Filter viewer docking)
- Report-DE direction: `differential-expression.ts` group-ordering + the DE
  dialog that currently carries the manual Comparison-dropdown workaround
- Up/down enrichment: `enrichment.ts` (two queries) + `enrichment-viewers.ts`
  (Direction-aware split charts)

</code_context>

<specifics>
## Specific Ideas

- Default significance metric = `adj.p-value` at the client's locked cut
  (Q<0.05, |log2FC|>0.58) so the Datagrok volcano matches the client's CK-omics
  figure out of the box; raw `p-value` is the alternate toggle state (D-06).
- The 11-category palette and keyword map are reproduced exactly from the
  CK-omics README (D-04) — it is shipped to the client as a contract; any
  divergence is a parity regression.
- R3 flip mechanics mirror CK-omics `create_subset_data`: invert `log2FC` sign
  and swap `AVG Group Quantity Numerator`/`Denominator` for rows whose declared
  comparison is the reverse of the canonical group1/group2 orientation.
- UniProt fields for the stream query: `accession,cc_subcellular_location,go_c`
  (cellular-component GO as the D-03 fallback signal).

</specifics>

<deferred>
## Deferred Ideas

- **Explicit "view the mirror" contrast control** — the CK-omics modal lets the
  analyst request the inverse of a declared contrast. Traded away in D-07 for
  the Filter-viewer UX. Revisit only if analysts hit it in practice.
- **Surfacing directional gene-list exports** (`*_up_genes_*.txt` /
  `*_down_genes_*.txt`) — CK-omics writes these; not required for parity of the
  in-product result. Planner's discretion, otherwise a future nicety.
- Note: the report-import DE direction fix is **NOT deferred** — it is in scope
  via D-09.

### Reviewed Todos (not folded)
- **`2026-03-03-expand-de-dialog-with-method-selection-comparison-picker-and-options.md`**
  — related (a "comparison picker") but distinct: that todo is about the **DE
  dialog**; R3/R4 here are about the **Candidates import path** and the
  **report-DE direction default**. D-09 addresses the orientation defect but not
  the broader DE-dialog expansion — that todo stays open for its own phase.
- The other ~19 phase-13 todo matches (titles/axis labels, imputation/
  normalization dialogs, import refactor, QC dashboard, FragPipe MaxLFQ, SPC,
  read-only sharing, streaming text parity) matched only on generic keywords
  (`proteomics`, `src`, `viewers`, `analysis`) at the matcher's noise floor —
  out of Phase-13 scope. Reviewed and not folded.

</deferred>

---

*Phase: 13-ck-omics-volcano-and-enrichment-parity*
*Context gathered: 2026-05-16*
