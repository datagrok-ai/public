# Phase 14: CK-omics Analyst-Experience Enhancements - Context

**Gathered:** 2026-05-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Analyst-facing UX and interpretation enhancements that match CK-omics' working experience, layered on Phase 13's parity foundation. Single phase covers two input streams:

**Stream A — ROADMAP R1–R5 (cluster B of the CK-omics gap analysis):**
- R1: predicted/uncharacterized gene-label improvement — detect ENSRNOG/ENSMUSG/LOC/ENSG/ENSDARG/MGP_/RGD/AABR IDs, resolve via Ensembl, mark provenance, keep original for hover, warn on duplicate descriptions (ref CK-omics `improve_gene_labels_with_ensrnog_marking`)
- R2: live filter-aware counters on the volcano (up/down/NS + by-location counts)
- R3: click-point → per-protein group-quantity detail in the UniProt context panel
- R4: group-mean correlation scatter (numerator-mean vs denominator-mean abundance, Pearson/Spearman, colored by significance class — distinct from QC sample×sample heatmap)
- R5: smart hierarchical pathway filtering — drop generic GO parents when specific children present, cap per source (ref CK-omics `apply_smart_pathway_filtering`)

**Stream B — 4 forwarded gaps from Phase 13 UAT (see `14-INBOX.md`):**
- G1 (major): volcano visual parity gaps — generic title, raw axis labels (negLog10P/log2FC), no default labels, "direction" generic semantics with red/blue, no default legend with counts, no protein search affordance
- G2: Volcano Options dialog shows defaults on re-open instead of current viewer state
- G3: Color → Location switch is slow and shows no progress indicator
- G4: Multi-contrast Filters viewer auto-includes the Spectronaut `Flags` column despite `columnNames: ['Comparison']` scoping

Scope is fixed: this discussion clarifies HOW, not WHETHER. New capabilities outside R1–R5 + G1–G4 belong in their own phase.

</domain>

<decisions>
## Implementation Decisions

### Scope packaging & sequencing
- **D-01:** **One phase, all in.** R1–R5 + G1–G4 land together in Phase 14. Volcano polish (G1+G2+G3) folds naturally into R2 (live counters on volcano) since both touch the same viewer; G4 (Filters scoping) is a 1-task fix.
- **D-02:** **MUST-have for BP DMD/WT client deliverable:** G1 (volcano visual parity), R1 (gene-label resolution), R2 (live counters). Planner should sequence these first. Nice-to-haves: R3, R4, R5, G2, G3, G4.

### Volcano polish (G1 + G2 + G3 + R2)
- **D-03:** **Default point labels:** auto-label top ~15–20 proteins by the active significance metric (adj.p-value default; p-value when toggled). Recompute on metric/filter/selection change. Matches CK-omics's visible-by-default top-hits labeling. Analyst can clear/edit via the search affordance (D-05).
- **D-04:** **Direction-label semantics from `proteomics.groups` tag.** Read group1/group2 names from the existing GroupAssignment JSON. Render legend categories as `"Enriched in <group1>"` / `"Enriched in <group2>"` / `"Not significant"`. **Color map is LOCKED:** magenta = group1 (numerator), cyan = group2 (denominator), gray = NS — mirrors CK-omics with no per-dataset config. Single code path works for both Spectronaut Candidates and non-Candidates (MaxQuant/text) flows.
- **D-05:** **Protein search via native Datagrok Filters viewer**, NOT a custom textarea. Dock a Filters viewer scoped to `Gene name` + `Protein ID` (string columns get a built-in search box for free). Matches wire to `df.selection.set(...)` — NOT `df.filter` — so highlighted points stay visible in the cloud and the volcano renders them with selection color in place, matching the CK-omics "Highlighted proteins" semantic. **This unifies with G4:** add `Gene name` and `Protein ID` to the same Filters viewer that already holds `Comparison (group1/group2)`, and explicitly exclude `Flags` from `columnNames`. Investigate whether the Filters viewer respects strict `columnNames` (D-07 of Phase 13 had it added Flags anyway); if not, set `columnNames` post-create via `viewer.setOptions(...)` or use `DG.Viewer.fromType('Filters', df, {...})` with the explicit columns + `meta`/exclusion options.
- **D-06:** **Live counter overlay** anchored bottom-right of the volcano (floating, not docked). Shows "Visible Proteins: Total: N / Enriched in <g1>: N / Enriched in <g2>: N / Not Significant: N" — and the by-location breakdown when Color = Subcellular Location. Recomputes on `df.filter` changes, `df.selection` changes, and viewer-property changes (metric/color toggle). Counts respect filtered-in rows only.

### Volcano polish — Claude's discretion (planner picks)
- Contrast-aware title source: synthesize from `proteomics.de_method` + group names → `"Volcano Plot: <Numerator> vs <Denominator>"`. Use Plotly title or a separate label widget.
- Axis label rewriting: `negLog10P` → `"-Log10(p-value)"` or `"-Log10(Q-value)"` (depending on active metric); `log2FC` → `"Log2 Fold Change"`. Live-rewrite on toggle.
- G2 (dialog state preload): source-of-truth between `viewer.getOptions()`, `df.tags`, and last-input cache — planner picks. Goal: re-opening Volcano Options shows the values currently applied.
- G3 (Color → Location speed + progress): wrap location-color path in `DG.TaskBarProgressIndicator` (or viewer-host `ui.setUpdateIndicator`). Profile the cache hit path — confirm we hit the cross-session `userDataStorage __schema_v` cache from Phase 13 D-02 rather than refetching; batch+memoize per-protein lookups; render incrementally so first paint isn't blocked by full-table classification.

### Gene-label resolution (R1)
- **D-07:** **Eager resolution at parse time, all rows.** New `Display Name` column populated by the parser's post-parse step before any viewer renders. Analyst never sees raw ENSRNOG IDs in the volcano. Matches CK-omics's pre-processed CSV behavior.
- **D-08:** **Inline provenance markers in `Display Name` string.** Examples: `Myh7*` for grouped, `Predicted-LOC123†` for predicted/reclassified. Marker travels with the gene name everywhere (volcano labels, heatmap rows, enrichment table, tooltips). Raw ID kept in `Source ID` column for hover/disambiguation. Single source-of-truth — no separate Provenance column.
- **D-09:** **Ensembl REST POST `/lookup/id` batched by species (~1000 IDs per request).** Group IDs by detected prefix → single POST per species → cross-session cache via `grok.userSettings` (mirror Phase-13 D-02 UniProt cache pattern, including `__schema_v` invalidation). **Three-level fallback:** Ensembl name → Ensembl description → raw ID if both empty. All traffic via `grok.dapi.fetchProxy()` (CORS-safe; matches Phase 13 D-01 UniProt stream pattern).
- **D-10:** **Duplicate descriptions: disambiguate by appending source ID.** If `Myh7` resolves from both `ENSRNOG00000001` and `ENSRNOG00000099`, render as `Myh7 (ENSRNOG00000001)` / `Myh7 (ENSRNOG00000099)`. Names stay unique in viewers; no silent merging. Emit `grok.shell.warning('N duplicate gene names disambiguated')` once after import.

### Per-protein panel + correlation scatter + smart pathway filter
- **D-11:** **Per-group quantities in UniProt panel: compact bar chart with mean ± SD whiskers.** Two/N bars (one per experimental group from `proteomics.groups`). Numeric mean+SD as small text below the bars. Extends `src/panels/uniprot-panel.ts`. Reuses canvas or a tiny svg helper — no new viewer.
- **D-12:** **Group-mean correlation scatter at `Proteomics | Visualize | Group-Mean Correlation…`** — sibling of Volcano Plot / Heatmap / PCA. Native Datagrok scatter viewer with two derived columns added to the protein DataFrame (`Numerator Mean`, `Denominator Mean`). Color by existing significance category column (matches D-04 magenta/cyan/gray). Pearson + Spearman correlations as title annotation (e.g. `"Group-Mean Correlation — r=0.93 (Pearson), ρ=0.91 (Spearman)"`). `createGroupMeanCorrelation()` factory in `viewers/` — same shape as `createVolcanoPlot`. Distinct from QC sample×sample heatmap (per ROADMAP).
- **D-13:** **Smart pathway filtering: port CK-omics `apply_smart_pathway_filtering` verbatim.** Treat `~/Downloads/ck/CKomics_tool2.py:apply_smart_pathway_filtering` as a locked client contract — same parent-detection heuristic (term is in ancestor list of any child in the result set), same per-source cap (top-N by FDR). Match-for-match port. Researcher reads the function during research; planner translates to TypeScript. Sources covered: GO (BP/CC/MF), KEGG, Reactome, WP — same rules apply per source.
- **D-14:** **Smart pathway filtering default-on with dialog toggle to disable.** Mirrors Phase 13 D-12 (WikiPathways default-on) pattern. Adds one checkbox to the Enrichment dialog. Most analysts benefit from the cleaner result set; power users doing methodology checks can disable.

### Claude's Discretion
- Exact wording/placement of axis label rewrites and the counter overlay.
- Storage mechanism for G2 dialog state preload (viewer.getOptions vs df.tags vs input cache).
- Progress indicator widget choice for G3 (TaskBarProgressIndicator vs viewer overlay).
- Ensembl batch size cap if 1000 turns out too aggressive for some species endpoints.
- Whether Display Name column replaces or augments the existing Gene name column (Claude reads the parsers and decides).
- How `*` vs `†` marker rules map to specific CK-omics cases — planner reads the reference function.
- Per-bar styling in the per-group quantity panel (height units, colors, error-bar variant).

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Locked client contracts (port verbatim)
- `~/Downloads/ck/CKomics_tool2.py` — source for R1 `improve_gene_labels_with_ensrnog_marking` (D-07/D-08/D-10) and R5 `apply_smart_pathway_filtering` (D-13). MUST READ before implementing R1 or R5.
- `~/Downloads/ck/Subcellular_Location_Classification_README.txt` — already locked by Phase 13 D-04. Relevant to G3 (Color → Location performance) since the 11-category palette + classification rules drive the path being optimized.
- `~/Downloads/ck/DMD_vs_WT/volcano_plots/` — target visual reference for G1 (volcano visual parity). Compare side-by-side after each volcano polish task.

### Phase 13 carry-forward (do NOT regress)
- `.planning/phases/13-ck-omics-volcano-and-enrichment-parity/13-CONTEXT.md` — D-01..D-12 of Phase 13 are LOCKED. Particularly: D-04 (11-color subcellular palette), D-05 (ONE volcano with property-toggle for color mode), D-07 (multi-contrast Filters viewer — D-05 here extends it), D-10/D-11 (split Up/Down enrichment — R5 here modifies the rows that flow into this layout).
- `.planning/phases/13-ck-omics-volcano-and-enrichment-parity/13-HUMAN-UAT.md` — full structured YAML for G1–G4 in the `## Gaps` section.
- `.planning/phases/14-ck-omics-analyst-experience-enhancements/14-INBOX.md` — forwarded-gaps summary table.

### Project-level
- `.planning/ROADMAP.md` Phase 14 entry — R1–R5 requirements (source of truth for Stream A scope).
- `packages/Proteomics/CLAUDE.md` — pipeline-stage tag conventions, parser/analysis/viewer file-layout contract, `findColumn`/`SEMTYPE` usage requirements.
- `packages/Proteomics/.planning/codebase/STRUCTURE.md` — "Where to Add New Code" contract for new parsers/analysis/viewers (D-12 adds a new viewer; R1 modifies parsers).

### Platform docs
- `grok.dapi.fetchProxy(...)` — CORS-safe external fetch (Ensembl REST in D-09).
- `grok.userSettings` — cross-session cache mechanism (D-09 cache pattern mirrors Phase 13 D-02).
- `DG.Viewer.filters(df, opts)` and `DG.Viewer.fromType('Filters', df, opts)` — for the unified Filter viewer in D-05/G4.
- `df.selection.set(idx, true)` — selection-driven highlight (D-05 wires search matches here, not `df.filter`).
- `DG.TaskBarProgressIndicator` / `ui.setUpdateIndicator` — progress UX for G3.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`createVolcanoPlot()` factory in `src/viewers/volcano.ts`** — extend in place for G1/D-03/D-04/D-06. Don't fork.
- **`dockComparisonFilterIfMultiContrast()` in `src/package.ts:117`** — extend `columnNames` to include `Gene name` + `Protein ID` (D-05); investigate strict-scoping fix to drop `Flags` (G4).
- **Phase 13 UniProt subcellular cache (`grok.userSettings` + `__schema_v`)** — D-09 mirrors this exact pattern for Ensembl labels.
- **`src/panels/uniprot-panel.ts`** — extend with the per-group bar chart (D-11). Already wired via `@grok.decorators.panel` with `semType` filter.
- **`src/utils/column-detection.ts:findColumn`** — use this for resolving `Display Name` / `Source ID` / group-quantity columns (per CLAUDE.md convention — never use raw `df.col('Gene names')`).
- **`src/utils/proteomics-types.ts:SEMTYPE`** — add new semantic types for `Display Name` / `Source ID` / `Numerator Mean` / `Denominator Mean` here AND mirror them in root-level `detectors.js` (single source of truth).
- **`createVolcanoPlot` color toggle pattern (Phase 13 D-05)** — D-04 direction semantics + D-06 counter overlay must layer onto this without forking.

### Established Patterns
- **Function-naming prefixes** (CLAUDE.md): `parseX` (parsers, return df), `showX` (dialogs), `openX` (table views), `createX` (viewer factories, no side effects). D-11 panel extension stays in `panels/`; D-12 follows `createX` for the new correlation viewer.
- **Pipeline tag convention:** `proteomics.*` tags coordinate workflow state. D-07 (eager Ensembl resolution at parse-time) does NOT need a new tag — `Display Name` column existence is the signal. D-13/D-14 (smart pathway filter) MAY set `proteomics.enrichment_smart_filtered` or similar — planner decides.
- **R fallback chain** (Phase 13 carry-forward, CLAUDE.md): always provide a client-side fallback when calling external services. R1's Ensembl resolution should degrade gracefully — if Ensembl is unreachable, render raw IDs with a single warning toast, do not block the import.
- **Heatmap clone-for-filter-isolation pattern** (CLAUDE.md): if any new viewer (D-12 correlation scatter) needs its own filter scope, use `df.clone(filter)` so it doesn't leak to volcano sharing the source DataFrame.

### Integration Points
- **Spectronaut Candidates parser** (`src/parsers/spectronaut-candidates-parser.ts`) — D-07 eager Ensembl resolution hooks in here for Candidates flow; same for `maxquant-parser.ts`, `spectronaut-parser.ts`, `fragpipe-parser.ts`, `generic-parser.ts`. All five parsers need the same post-parse Ensembl-resolution call (extract to `src/utils/gene-label-resolver.ts` or similar — shared helper, not duplicated).
- **`uniprot-panel.ts`** — currently fetches per-protein metadata on cell click. R3 (D-11) adds the per-group bar chart to the same panel; reuse the existing fetch trigger.
- **Enrichment dialog** (`src/analysis/enrichment.ts` — or wherever the dialog lives) — D-14 adds one checkbox ("Apply smart pathway filter (default on)"); D-13 logic runs post-g:Profiler-response before the DataFrame is returned.
- **Volcano + Filters viewer dock arrangement** — D-05 extends the existing Phase 13 D-07 Filter viewer rather than docking a second one. Result: ONE Filters viewer with `[Comparison (group1/group2), Gene name, Protein ID]` (minus `Flags`).

</code_context>

<specifics>
## Specific Ideas

- **CK-omics figure as visual ground-truth:** the volcano in `~/Downloads/ck/DMD_vs_WT/volcano_plots/` is the comparison standard for G1. Phase 13 confirmed distribution/orientation match; Phase 14 closes the surface-finish gap (title, axes, labels, semantics, legend, search, counters).
- **User explicitly preferred Datagrok-native Filters viewer over a custom search textarea** (D-05). Rationale: reuses platform muscle memory, gets the search box for free, wires to `df.selection` for in-place highlight without adding a new component.
- **Color lock is non-negotiable** (D-04): magenta = numerator-enriched, cyan = denominator-enriched, gray = NS — same hex values across all volcanos in the package. Mirrors how the subcellular palette is treated in Phase 13 D-04.

</specifics>

<deferred>
## Deferred Ideas

- **Custom volcano JsViewer** — already in PROJECT.md "Out of Scope". Phase 14 keeps using Datagrok's native scatter with formula lines (Phase 13 D-05/D-06 lock-in).
- **Per-protein waterfall chart in the UniProt panel** (vs the simpler bar chart in D-11) — captured for a future enhancement if D-11's bars prove insufficient.
- **Configurable per-source pathway cap exposed in the Enrichment dialog** (alternative to D-14's binary toggle) — captured if analysts later want fine-grained control.
- **Group-mean correlation viewer fold into QC dashboard** (alternative to D-12 standalone menu entry) — preserved if the dashboard ever gets a "diagnostic correlations" tab.
- **Provenance-driven color or filter in viewers** (alternative to D-08's inline-marker-only approach) — captured if the inline `*`/`†` proves too subtle.
- **Methodology-check / power-user mode for raw enrichment** (alternative to D-14 default-on toggle) — captured if a separate "research" enrichment flow is later needed.

</deferred>

---

*Phase: 14-ck-omics-analyst-experience-enhancements*
*Context gathered: 2026-05-28*
