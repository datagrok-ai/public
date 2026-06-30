---
phase: 15
plan: 02
status: complete
type: implementation
---

# Plan 15-02 Summary — trim-dataframe primitives

Built `src/publishing/trim-dataframe.ts` exporting `trimForPublish` and `trimEnrichmentForPublish`. These are the deep-clone + allowlist primitives the Plan 04 orchestrator hands the live source DataFrame to — the orchestrator never touches the source directly.

## Allowlist resolution

Used `findProteomicsColumns` (returns `proteinId`, `geneName`, `log2fc`, `pValue`) plus three additional `findColumn` calls for `adj.p-value` / `significant` / `direction` since `findProteomicsColumns` does not surface those. All seven columns are looked up by semType-first, name-hint-fallback (CLAUDE.md convention). Hard-fail with a descriptive error on missing required columns; `geneName` and `direction` are soft-optional and dropped from the allowlist if absent.

Final allowlist size: 5–7 columns depending on `geneName` and `direction` presence.

## Tag survival

Spike 15-00 confirmed all 14 `proteomics.*` + `proteomics.published*` tags survive `DG.Project.save → find → open` on `release/1.27.3` — but the spike did NOT independently test `df.clone()` tag preservation. The defensive re-set via `setPublishedTags(frozen, meta)` after the clone is no-cost in either case (re-setting an already-present tag is a no-op).

Carry-forward tags written to the clone:
- `proteomics.source`, `proteomics.de_method`, `proteomics.groups`, `proteomics.de_complete`

All Phase-15 published tags written via `setPublishedTags` (Plan 01 helper).

## SemType re-assignment

Spike 15-00 showed semTypes survive `DG.Project.save → find → open` intact. However, `df.clone(rowMask, columnIds)` produces a fresh column list and the platform may or may not preserve the `semType` field on the cloned columns. Re-assigning is defensive and free — done in Step E for the 4 SEMTYPE-namespaced columns (PROTEIN_ID, GENE_SYMBOL, LOG2FC, P_VALUE).

## Belt-and-braces metadata columns (PUB-11)

13 `_meta_*` columns written via `col.init(() => constant)` per saved memory `feedback_dg_column_bulk_init.md`. All columns marked `.hidden=true` so they do not appear in the reviewer's default grid view (per RESEARCH §"Pattern 3" note 5).

Per saved memory `feedback_dg_column_init_null_sentinel.md`: numeric columns with null init values leave the FLOAT_NULL sentinel `2.6789e-34` unsynced. All metadata fields are non-null by construction (the orchestrator populates `PublishedMetadata` fully before calling `trimForPublish`), so this trap does not bite us. Optional fields (`publishedByEmail`, `supersedes`, `supersededBy`) are stored as empty strings when null per the `emptyForNull` flag, since they are typed `string`.

## Enrichment trim

`trimEnrichmentForPublish` mirrors `trimForPublish` for the enrichment DataFrame's allowlist (Term Name, Source, p-value, adj.p-value, Intersection, optional Direction from Phase 13). Writes the `proteomics.enrichment=true` tag + 2 belt-and-braces columns + 3 critical published tags. The Plan 04 orchestrator calls this opportunistically only when the source carries `proteomics.enrichment` AND an enrichment DataFrame is present in `grok.shell.tables` (D-05 opportunistic carry).

## Source-DF non-mutation

`grep -E "source\.setTag|source\.col\(.*\)\.set|source\.columns\." src/publishing/trim-dataframe.ts` returns zero matches. Every write goes through `frozen` (the clone). T-15-02 mitigation: source is not mutated.

## Verification

Project-wide `tsc --noEmit` passes clean. No hard-coded column lookups (`grep -E "df\.col\('[A-Z]"` returns empty).

## Output

`src/publishing/trim-dataframe.ts` — 209 lines, type-checks under strict mode.
