---
phase: 15-read-only-publishing-foundation
plan: 02
type: execute
wave: 2
depends_on: ["15-00", "15-01"]
files_modified:
  - src/publishing/trim-dataframe.ts
autonomous: true
requirements: [PUB-01, PUB-02, PUB-03, PUB-11]
must_haves:
  truths:
    - "Calling `trimForPublish(df)` returns a deep clone of `df` containing ONLY the 7-column allowlist (Protein ID, Gene, log2FC, p-value, adj.p-value, significant, direction) — raw intensities and peptide counts are absent"
    - "Mutating the source `df` AFTER calling `trimForPublish` does NOT change the returned clone (deep-clone, Pitfall 1)"
    - "The returned clone has all 13 `proteomics.published*` tags re-set explicitly post-clone (Pitfall 3 mitigation)"
    - "The returned clone has 13 single-row `_meta_published_*` columns matching the tags (belt-and-braces, PUB-11)"
    - "Column lookup goes through `findColumn(df, semType, nameHints)` — NEVER `df.col('Gene names')` raw (CLAUDE.md convention)"
  artifacts:
    - path: "src/publishing/trim-dataframe.ts"
      provides: "trimForPublish + trimEnrichmentForPublish"
      exports: ["trimForPublish", "trimEnrichmentForPublish"]
  key_links:
    - from: "src/publishing/trim-dataframe.ts"
      to: "src/utils/column-detection.ts"
      via: "findColumn import"
      pattern: "import.*findColumn.*column-detection"
    - from: "src/publishing/trim-dataframe.ts"
      to: "src/publishing/publish-state.ts"
      via: "setPublishedTags + META_COLUMNS + PUBLISHED_TAGS imports"
      pattern: "import.*publish-state"
---

<objective>
Implement `trimForPublish(df, meta)` and `trimEnrichmentForPublish(df, meta)`: the deep-clone-first primitives that produce frozen DataFrame snapshots ready for `DG.Project.save`. Combines three v1.3 patterns — `heatmap.ts` clone-for-isolation + `column-detection.ts` semType-first lookup + Plan 01's tag namespace — into the single, idempotent contract the orchestrator (Plan 04) calls.

Purpose: Pitfall 1 (stale-snapshot leak) AND Pitfall 3 (tag-stripping serializer) are mitigated HERE, not in `publishAnalysis`. The orchestrator never touches the source DF directly — it hands the source to `trimForPublish` and operates on the returned clone. This is the only file where `df.clone(...)` lives for publishing.

Output: `src/publishing/trim-dataframe.ts` exporting two functions: protein-DF trim (always called) and enrichment-DF trim (called opportunistically when source has `proteomics.enrichment === 'true'` and enrichment DF is in `grok.shell.tables` per D-05).
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md
@.planning/phases/15-read-only-publishing-foundation/15-01-SUMMARY.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/utils/proteomics-types.ts
@packages/Proteomics/src/utils/column-detection.ts
@packages/Proteomics/src/viewers/heatmap.ts
@packages/Proteomics/src/viewers/qc-computations.ts
@packages/Proteomics/src/publishing/publish-state.ts

<interfaces>
<!-- Pulled from src/utils/column-detection.ts -->
From src/utils/column-detection.ts:
  export function findColumn(df: DG.DataFrame, semType: string, nameHints: string[]): DG.Column | null
  export function findProteomicsColumns(df: DG.DataFrame): {
    proteinId: DG.Column | null;
    geneSymbol: DG.Column | null;
    log2FC: DG.Column | null;
    pValue: DG.Column | null;
    adjPValue: DG.Column | null;
    significant: DG.Column | null;
  }

From src/utils/proteomics-types.ts:
  SEMTYPE.PROTEIN_ID = 'Proteomics-ProteinId'
  SEMTYPE.GENE_SYMBOL = 'Proteomics-GeneSymbol'
  SEMTYPE.LOG2FC = 'Proteomics-Log2FC'
  SEMTYPE.P_VALUE = 'Proteomics-PValue'
  SEMTYPE.ADJ_P_VALUE = 'Proteomics-AdjPValue'
  SEMTYPE.SIGNIFICANT = 'Proteomics-Significant'
  (Plan 07 may add SEMTYPE.DIRECTION conditionally)

From datagrok-api:
  df.clone(filter?: DG.BitSet, columns?: string[] | DG.Column[]): DG.DataFrame
  df.columns.addNewString(name: string): DG.Column
  df.columns.addNewDateTime(name: string): DG.Column
  df.columns.addNewFloat(name: string): DG.Column
  df.columns.addNewInt(name: string): DG.Column
  col.init(callback: (rowIdx: number) => any): DG.Column
  col.setTag(key: string, value: string): void
  // saved memory: col.init(()=>null) on numeric leaves FLOAT_NULL sentinel — use isNone() check, NOT raw get()
  // saved memory: DG.Column bulk fills via init() are 137x-255x faster than per-row col.set
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Implement trimForPublish — clone + allowlist + tags + metadata columns</name>
  <files>src/publishing/trim-dataframe.ts</files>
  <read_first>
    - @packages/Proteomics/src/viewers/heatmap.ts (lines 56-60 — `df.clone(filter)` precedent)
    - @packages/Proteomics/src/viewers/qc-computations.ts (lines 27-32 — `ensureFreshFloat` idempotency)
    - @packages/Proteomics/src/utils/column-detection.ts (full file — `findColumn` + `findProteomicsColumns` API)
    - @packages/Proteomics/src/utils/proteomics-types.ts (full file — SEMTYPE constants)
    - @packages/Proteomics/CLAUDE.md ("Use findColumn / findProteomicsColumns, never by raw name")
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 2 — exact patterns + notes)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"DG.DataFrame.clone with column allowlist + tag reset" + §"Pattern 3")
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (Phase 15-specific pitfall — volcano color binding columns must be in allowlist)
    - @.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md (spike output — confirms which tags need explicit re-set post-clone)
  </read_first>
  <action>
Create `src/publishing/trim-dataframe.ts`. Standard imports (`* as grok`, `* as DG`, `SEMTYPE` from `../utils/proteomics-types`, `findColumn` + `findProteomicsColumns` from `../utils/column-detection`, `PUBLISHED_TAGS`, `META_COLUMNS`, `setPublishedTags`, `PublishedMetadata` from `./publish-state`).

Export `trimForPublish(source: DG.DataFrame, meta: PublishedMetadata): DG.DataFrame`.

**Step A — resolve allowlist columns via findColumn (NEVER `df.col('Gene names')`):**
Use `findProteomicsColumns(source)` to get `proteinId`, `geneSymbol`, `log2FC`, `pValue`, `adjPValue`, `significant`. For `direction`, fall back to `findColumn(source, SEMTYPE.DIRECTION ?? 'unused', ['direction', 'regulation', 'up_down'])` (SEMTYPE.DIRECTION may not exist — pass `''` or `undefined` to skip semType check and rely on name hints; coordinate with Plan 07 conditional SEMTYPE add).

If any of `proteinId` / `log2FC` / `pValue` / `adjPValue` / `significant` is null, throw `Error('Cannot publish: missing required column [<name>]. Run Differential Expression to completion before publishing.')` — defensive, since `requireDifferentialExpression` in Plan 07 menu handler should have already gated. `geneSymbol` and `direction` are SOFT — if null, do NOT include them in allowlist; downstream volcano falls back gracefully (matches Phase 14 D-04 behavior).

Build the allowlist as an array of column NAMES (not Column objects): `[proteinId.name, geneSymbol?.name, log2FC.name, pValue.name, adjPValue.name, significant.name, direction?.name].filter(Boolean)`.

**Step B — deep clone via the allowlist form of df.clone:**
`const frozen = source.clone(undefined, allowlist);` (filter=undefined keeps every row; columns=allowlist is the allowlist form per RESEARCH §"DG.DataFrame.clone with column allowlist + tag reset"). This is the Pitfall 1 mitigation point — the clone is fully independent of the source thereafter.

**Step C — re-set every required `proteomics.*` tag explicitly (Pitfall 3 mitigation):**
Per spike output in `15-00-SUMMARY.md`: if tags DID survive clone, this is defensive. If tags did NOT survive, this is required. Either way, re-set them — there is no cost to re-setting:
  - Preserve source-pipeline tags: `proteomics.source`, `proteomics.de_method`, `proteomics.groups`, `proteomics.de_complete` (copy from source to frozen)
  - Write Phase 15 published tags via `setPublishedTags(frozen, meta)` (Plan 01 helper)
Use `frozen.setTag(key, source.getTag(key) ?? '')` for the carry-forward source tags; skip nulls.

**Step D — set the frozen DataFrame's name (per RESEARCH Assumption A5 + Pattern 1):**
`frozen.name = \`${source.name || 'analysis'}_published_${meta.publishedAt.toISOString().slice(0,10)}\`;`

**Step E — re-assign semTypes on the cloned columns (per spike output A6):**
If spike confirmed semTypes are stripped on clone (likely), re-assign:
```
frozen.col(proteinId.name)!.semType = SEMTYPE.PROTEIN_ID;
frozen.col(geneSymbol.name)!.semType = SEMTYPE.GENE_SYMBOL; // if present
frozen.col(log2FC.name)!.semType = SEMTYPE.LOG2FC;
frozen.col(pValue.name)!.semType = SEMTYPE.P_VALUE;
frozen.col(adjPValue.name)!.semType = SEMTYPE.ADJ_P_VALUE;
frozen.col(significant.name)!.semType = SEMTYPE.SIGNIFICANT;
```

**Step F — write 13 belt-and-braces metadata columns (PUB-11 + RESEARCH §"Pattern 3"):**
Per saved memory `feedback_dg_column_bulk_init.md`: use `col.init(() => constant)` form, NOT per-row set. Per saved memory `feedback_dg_column_init_null_sentinel.md`: when a numeric column's init callback returns null, the underlying buffer holds the FLOAT_NULL sentinel `2.6789e-34` — reader must use `col.isNone(rowIdx)` not raw `col.get(rowIdx)`. The metadata fields are ALL non-null by construction (PublishedMetadata is fully populated by the orchestrator), so this trap does not bite us — but document the pattern in JSDoc anyway.

Write the columns using helper `addMetadataColumn(frozen, name, value, type)`:
  - `_meta_published_target` (string): `meta.target` (RAW user input, not slug)
  - `_meta_published_at` (datetime): `meta.publishedAt`
  - `_meta_published_by` (string): `meta.publishedBy`
  - `_meta_published_by_email` (string): `meta.publishedByEmail ?? ''`
  - `_meta_published_de_method` (string): `meta.deMethod`
  - `_meta_published_fc_threshold` (float): `meta.fcThreshold`
  - `_meta_published_p_threshold` (float): `meta.pThreshold`
  - `_meta_published_version` (int): `meta.version`
  - `_meta_published_id` (string): `meta.publishId`
  - `_meta_published_includes_enrichment` (string): `meta.includesEnrichment ? 'true' : 'false'`
  - `_meta_supersedes` (string): `meta.supersedes ?? ''`
  - `_meta_superseded_by` (string): `meta.supersededBy ?? ''`

Use the `ensureFreshFloat`-style idempotency: `if (frozen.columns.contains(name)) frozen.columns.remove(name);` before `addNew...` (cheap insurance since clone shouldn't contain `_meta_*` columns but defensive).

Mark each metadata column hidden from the default grid: `col.setTag('.hidden', 'true')` (per RESEARCH §"Pattern 3" note 5 — reviewer doesn't see metadata as data).

Return `frozen`.

JSDoc the function with: "Returns a deep clone of `source` with only the 7-column allowlist + 13 belt-and-braces metadata columns + every required `proteomics.published*` tag re-set explicitly. Source DataFrame is NOT mutated."
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/trim-dataframe.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export function trimForPublish" src/publishing/trim-dataframe.ts | grep -q '^1$' &amp;&amp; grep -c "findProteomicsColumns" src/publishing/trim-dataframe.ts | grep -qv '^0$' &amp;&amp; grep -c "source\.clone(undefined" src/publishing/trim-dataframe.ts | grep -qv '^0$' &amp;&amp; grep -c "setPublishedTags" src/publishing/trim-dataframe.ts | grep -qv '^0$' &amp;&amp; { ! grep -E "df\\.col\\('[A-Za-z]" src/publishing/trim-dataframe.ts; }</automated>
  </verify>
  <done>`trimForPublish` exported; uses `findProteomicsColumns` for allowlist resolution; calls `source.clone(undefined, allowlist)`; calls `setPublishedTags(frozen, meta)`; writes 13 `_meta_*` columns; assigns SEMTYPEs on cloned columns; sets `frozen.name`; never references raw column names via `df.col('SomeName')`; TypeScript strict compiles.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 2: Implement trimEnrichmentForPublish (D-05 opportunistic carry)</name>
  <files>src/publishing/trim-dataframe.ts</files>
  <read_first>
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-05 enrichment carry — opportunistic, no-enrichment source still publishes cleanly)
    - @packages/Proteomics/src/viewers/enrichment-viewers.ts (enrichment-DF column conventions)
    - @packages/Proteomics/src/analysis/enrichment.ts (enrichment-DF column names produced by the analysis step)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 2 — clone+idempotency pattern reused)
  </read_first>
  <action>
Add to `src/publishing/trim-dataframe.ts`:

`trimEnrichmentForPublish(enrichSource: DG.DataFrame, meta: PublishedMetadata): DG.DataFrame` — same shape as `trimForPublish` but with the enrichment DataFrame's allowlist.

Look up the enrichment DataFrame's columns:
  - Term Name (string)
  - Source (`'GO:BP'` / `'KEGG'` / `'REAC'` etc; string)
  - p-value (float)
  - adj.p-value (FDR) (float)
  - Intersection / intersection-genes (string)
  - Direction (if Phase 13 split present — Up/Down/Mixed) (string, OPTIONAL — soft-fail if absent)

Use `findColumn(enrichSource, '', ['term name', 'term'])` etc. (no SEMTYPE for enrichment columns — name hints only). For required columns (Term Name, Source, p-value, adj.p-value, Intersection) throw with clear error if absent. For Direction, omit from allowlist if null.

Clone via `enrichSource.clone(undefined, allowlist)`.

Carry-forward source enrichment tag: `enrichClone.setTag('proteomics.enrichment', 'true');` (single tag — enrichment DataFrames only have one workflow tag per CLAUDE.md tag table).

Set `enrichClone.name = \`enrichment_published_${meta.publishedAt.toISOString().slice(0,10)}\``;

Write 2 belt-and-braces metadata columns (minimal — enrichment has its own metadata column from the protein DF in same Project):
  - `_meta_published_id` (string) — link back to publish id
  - `_meta_published_includes_enrichment` (string) — sentinel `'true'`

JSDoc: "Returns a deep clone of `enrichSource` containing only the term/source/p/adj.p/intersection allowlist. Called by orchestrator only when source DataFrame has `proteomics.enrichment === 'true'` AND an enrichment DF is in `grok.shell.tables` (D-05 opportunistic carry)."

(Discovery of the enrichment DF in the shell happens in Plan 04 orchestrator — this helper is pure: given the enrichment DF, trim it.)
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/trim-dataframe.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export function trimEnrichmentForPublish" src/publishing/trim-dataframe.ts | grep -q '^1$' &amp;&amp; grep -c "enrichSource\.clone(undefined" src/publishing/trim-dataframe.ts | grep -qv '^0$'</automated>
  </verify>
  <done>`trimEnrichmentForPublish` exported; uses `findColumn` (never raw `df.col`); clones with allowlist; sets `proteomics.enrichment` tag; sets `name`; TypeScript strict compiles.</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| live source DF -> trimForPublish (clone boundary) | Source DF MUST NOT be mutated; clone is the new trust domain |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-02 | Information disclosure (stale-snapshot leak) | trimForPublish | mitigate | `source.clone(undefined, allowlist)` is the FIRST mutation operation in the function; no `source.setTag` / `source.columns.remove` / `source.col(...).set` calls anywhere in this file. Plan 08 round-trip test asserts source-unchanged after `trimForPublish`. |
| T-15-03 | Tampering (target value injection into column data) | metadata column writers | accept | Raw target string preserved in `_meta_published_target` column for biologist audit context; column is HIDDEN from default grid (`.hidden=true` tag); never used in URL/SQL/shell context |
</threat_model>

<verification>
- TypeScript strict-mode compiles
- `grep -E "df\.col\('[A-Z]"` returns no hits (no hard-coded column names)
- `grep "source\.clone(undefined"` present
- `grep "setPublishedTags"` present
- (Smoke test, optional) Build a minimal DF in REPL, call `trimForPublish`, verify column count + tag set
</verification>

<success_criteria>
- `src/publishing/trim-dataframe.ts` exports `trimForPublish` and `trimEnrichmentForPublish`
- Source DF is never mutated (every operation goes through the clone)
- 13 `_meta_*` columns written via `init(() => constant)` (per saved memory bulk-fill pattern)
- Every column lookup uses `findColumn` / `findProteomicsColumns` per CLAUDE.md
- Required tag namespace re-set post-clone via `setPublishedTags(frozen, meta)`
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-02-SUMMARY.md` when done with:
- Allowlist columns resolved (which 7 / 6 if no direction / etc.)
- Whether Spike output confirmed that source clone preserves tags (defensive re-set is no-cost either way)
- Whether semType re-assignment was needed post-clone (per spike A6 output)
</output>
