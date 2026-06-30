# Phase 12: Spectronaut Input Coverage - Context

**Gathered:** 2026-05-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Stream a multi-GB Spectronaut **precursor-level long report** (one row per peptide × precursor × fragment × run) into the existing pivot/semType/log2/auto-group logic, entirely client-side with bounded memory, so it imports through the existing `Proteomics | Import | Spectronaut Report` menu item and drives the existing pipeline (Annotate → Normalize → Impute → DE → Volcano/Heatmap → Enrichment) — with no manual pre-processing step.

</domain>

<spec_lock>
## Requirements (locked via SPEC.md)

**6 requirements are locked.** See `12-SPEC.md` for full requirements, boundaries, and acceptance criteria.

Downstream agents MUST read `12-SPEC.md` before planning or implementing. Requirements are not duplicated here.

**In scope (from SPEC.md):**
- Streaming precursor-level long-report import through the existing `Proteomics | Import | Spectronaut Report` menu item
- Client-side TypeScript aggregation using browser-native `File.stream()` / `TextDecoderStream` (no new npm deps)
- Filter + aggregate parity with `/tmp/spectronaut-aggregate.sql` (EG.Qvalue ≤ 0.01, CON__/REV__ drop, max PG.Quantity, min EG.Qvalue per group)
- `TaskBarProgressIndicator` and event-loop yielding during import
- Committing the duckdb script into the package as a documented fallback
- A small derived precursor fixture + automated test

**Out of scope (from SPEC.md):**
- Pivot/Matrix Report support — deferred to a later phase; distinct schema problem, no fixture yet
- Peptide-level quantification/analysis — prior project decision; aggregate to protein-group level only
- Server-side aggregation script (Python/R) — client-side streaming dominates on deps/UX/env
- Web Worker offload — resolved as a HOW decision below (not used; see D-02)
- Spectronaut spectral-library exports — not DE inputs
- Changes to the downstream pipeline (Normalize/Impute/DE/viewers) — unchanged; phase only adds an input path

</spec_lock>

<decisions>
## Implementation Decisions

### Format routing (same menu item, two shapes)
- **D-01:** Header-sniff then branch. Read only the header line first; if precursor-level signature columns are present (any of `EG.ModifiedPeptide`, `FG.Charge`, `PEP.StrippedSequence` — columns a PG-level long report does not carry) → take the new streaming-aggregate path; otherwise → keep the existing `file.text()` in-memory path unchanged. Rationale: keeps the proven small-file path fast and behaviorally identical; size-threshold rejected because file size does not reliably indicate report shape.

### Execution model
- **D-02:** Main thread + await-yield. Process the `File.stream()` reader in an async chunk loop; the reader's `await` boundaries yield the event loop between chunks. Drive a `DG.TaskBarProgressIndicator` by bytes-read. No Web Worker — avoids webpack worker bundling against Datagrok's constrained externals and message-passing complexity for a ~30 s job. This is the same responsiveness pattern that resolved the DE Page-Unresponsive issue earlier.

### Test fixture
- **D-03:** Synthetic generator. Extend the existing `makeLongFormatTsv` helper in `src/tests/spectronaut-parser.ts` to emit precursor-level rows: multiple precursor/fragment rows per (protein × sample), `PEP.*`/`EG.*`/`FG.*` columns, plus deliberate `CON__`/`REV__` decoy rows and sub-threshold `EG.Qvalue` rows so every filter branch is exercised. Fully reproducible, no client-data licensing concerns.
- **D-04:** Equivalence oracle = committed golden file. Run the duckdb script (`/tmp/spectronaut-aggregate.sql`, to be moved per D-05) ONCE over the synthetic fixture; commit its output as a golden `.tsv` beside the fixture. The streaming test asserts streaming output == golden. Provide a documented regeneration command (e.g. a make/npm target) for when the fixture changes. Pins behavior to the *real* duckdb logic, not a hand-derived or re-ported approximation.

### Duckdb fallback location
- **D-05:** New `packages/Proteomics/tools/` directory holding the `.sql` + `.sh`. Documented in the package README. The streaming-import failure path surfaces a hint pointing at it (toast text and/or README link). Rejected `files/` (runtime-served, semantically wrong for a dev/CLI tool) and docs-only (no in-product pointer).

### Claude's Discretion
- Chunk size for the stream reader, line-buffer carry-over handling across chunk boundaries, exact progress-update cadence, and the precise wording of the failure-path hint — standard implementation choices for research/planning.

### Reviewed Todos (not folded)
- "Add titles to various plots" (2026-03-03), "Auto-zoom volcano plot…" (2026-03-05), "Improve chart titles and axis labels across all viewers" (2026-03-10) — all matched only on `area: ui` (score 0.3); unrelated to a streaming TSV importer. Reviewed and deferred.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Locked requirements
- `.planning/phases/12-spectronaut-input-coverage/12-SPEC.md` — Locked requirements, boundaries, acceptance criteria — MUST read before planning

### Existing import / parser code (the integration surface)
- `packages/Proteomics/src/package.ts` §`importSpectronaut` (~lines 99-115) — the menu handler to modify; currently `file.text()` → `parseSpectronautText`
- `packages/Proteomics/src/parsers/spectronaut-parser.ts` — existing PG-level long-form parser: `REQUIRED_COLUMNS`, `QUANTITY_COLUMNS` (PG.IBAQ→PG.Quantity), `pivotSpectronaut`, `buildWideDataFrame`, `autoPopulateGroups`; the streaming path must feed this same pivot/semType/log2/group logic
- `packages/Proteomics/src/parsers/shared-utils.ts` — log2 transform / log2-status detection / delimiter detection / primary-id extraction reused by the parser
- `packages/Proteomics/src/tests/spectronaut-parser.ts` §`makeLongFormatTsv` — synthetic TSV builder to extend per D-03

### Aggregation reference logic
- `/tmp/spectronaut-aggregate.sql` and `/tmp/spectronaut-aggregate.sh` — proven duckdb aggregation to be ported to streaming TS and relocated per D-05; defines the exact filters/aggregates the TS path must match (EG.Qvalue ≤ 0.01 non-numeric-passes, CON__/REV__ drop, max PG.Quantity & min EG.Qvalue per protein×condition×replicate)

### Codebase maps
- `.planning/codebase/STRUCTURE.md` — directory contract (parsers in `src/parsers/`, tests in `src/tests/`, fixtures in `files/demo/`)
- `.planning/codebase/TESTING.md` — `grok test` conventions and fixture patterns
- `.planning/codebase/CONVENTIONS.md` — package coding conventions

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `pivotSpectronaut` / `buildWideDataFrame` / `autoPopulateGroups` in `spectronaut-parser.ts` — the streaming path produces the same aggregated long-form rows these already consume; ideally the streaming aggregator emits a structure the existing pivot can ingest unchanged (or refactor `parseSpectronautText` to accept a pre-aggregated row source vs. raw text)
- `makeLongFormatTsv` test helper — extend rather than write a new fixture builder (D-03)
- `DG.TaskBarProgressIndicator` — used elsewhere in the package (DE dialog) for progress; reuse for byte-read progress (D-02)
- `QUANTITY_COLUMNS` already accepts `PG.Quantity` (shipped commit `fc139e3af3` this milestone) — precursor reports carry `PG.Quantity`, so no quantity-column work needed

### Established Patterns
- `parseXText(text) → DG.DataFrame` parser contract (caller sets `df.name`, adds to shell) — the streaming path breaks the `text` assumption; planner must decide whether to add a sibling `parseSpectronautStream(file)` or restructure around a row-iterator the existing pivot consumes
- Bridge-hop avoidance (`Column.init` / `getRawData`, memory `feedback_dg_column_bulk_init`) — the final wide-DataFrame construction must use bulk init, not per-row `col.set`
- Tag contract: streaming path must set the same `proteomics.source=spectronaut`, `proteomics.preNormalized=true`, and auto-populate `proteomics.groups` as the existing parser

### Integration Points
- `Proteomics | Import | Spectronaut Report` menu handler in `package.ts` — single entry point; header-sniff branch (D-01) lives here
- Existing small-file path must remain byte-for-byte behaviorally unchanged (regression risk) — the branch must be additive

</code_context>

<specifics>
## Specific Ideas

- Precursor-signature detection columns (D-01): `EG.ModifiedPeptide`, `FG.Charge`, `PEP.StrippedSequence` — presence of any one indicates a precursor-level report. The reference file `2026-05-13 BP DMD WT.tsv` header (76 cols) contains all three; a PG-level long report (the existing HYE-mix demo) contains none.
- Reference (manual release-gate) file: `~/Downloads/2026-05-13 BP DMD WT.tsv` (2.6 GB / 2.88 M rows → 8,328 proteins × 24 samples). Not committable; CI uses the synthetic fixture (D-03).
- Golden-file regeneration must be a documented one-liner so the fixture↔golden pair stays in sync (D-04).

</specifics>

<deferred>
## Deferred Ideas

- **Pivot/Matrix Report support** — already-wide protein × sample Spectronaut export; rejected by `parseSpectronautText` (no `R.Condition`/`R.Replicate`). Distinct schema problem, no fixture yet. Own future phase.
- **Peptide-level quantification** — standing project decision to defer (doubles parser complexity for marginal gain); this phase stays at protein-group aggregation.

</deferred>

---

*Phase: 12-spectronaut-input-coverage*
*Context gathered: 2026-05-15*
