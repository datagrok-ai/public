# Phase 12: Spectronaut Input Coverage — Specification

**Created:** 2026-05-15
**Ambiguity score:** 0.134 (gate: ≤ 0.20)
**Requirements:** 6 locked

## Goal

A user can select a raw multi-GB Spectronaut **precursor-level long report** (one row per peptide × precursor × fragment × run) via `Proteomics | Import | Spectronaut Report` and the package streams + aggregates it client-side into the same wide protein × sample DataFrame the existing pipeline consumes — without loading the whole file into memory and without any manual pre-processing step.

## Background

`importSpectronaut` (`src/package.ts:99-115`) reads the file with `file.text()`, which materializes the entire file as one JS string — fatal on a 2.6 GB input (V8 max string length ≈ 512 MB). `parseSpectronautText` (`src/parsers/spectronaut-parser.ts`) expects an **already-PG-level** long report (`R.Condition`, `R.Replicate`, `PG.ProteinGroups`, `PG.IBAQ`/`PG.Quantity`) and pivots it to wide.

The precursor-level report repeats every protein-group metric across millions of precursor/fragment rows. The proven workaround is an external duckdb script (`/tmp/spectronaut-aggregate.sql` + `.sh`) that collapses the reference file `2026-05-13 BP DMD WT.tsv` (2.6 GB / 2.88 M rows) to ~198 K long-form rows (8,328 proteins × 24 samples, ~21 MB), applying an `EG.Qvalue ≤ 0.01` precursor filter and `CON__/REV__` decoy drop, taking `max(PG.Quantity)` and `min(EG.Qvalue)` per (protein × condition × replicate). After this manual step, the existing pipeline (Annotate → Normalize → Impute → DE → Volcano/Heatmap → Enrichment) runs end-to-end and was verified manually.

Key insight: the aggregation output is bounded by **proteins × samples**, not by input row count, because `PG.Quantity` is constant per (protein, run). The in-memory aggregation state is therefore tens of MB regardless of input size — making a browser-native streaming TypeScript implementation feasible with zero new dependencies.

## Requirements

1. **Streaming precursor importer**: The Spectronaut Report import path ingests a precursor-level long report via incremental streaming, never materializing the whole file.
   - Current: `importSpectronaut` calls `file.text()`; fails on 2.6 GB with a V8 string-length / OOM error; `parseSpectronautText` only accepts PG-level long form
   - Target: the import handler reads via `File.stream()` → `TextDecoderStream` → incremental line processing, aggregating into a `Map` keyed by (protein × condition × replicate), then feeds the existing pivot/semType/log2/auto-group logic
   - Acceptance: selecting the 2.6 GB reference TSV completes import with no string-length or out-of-memory error

2. **Aggregation equivalence with the duckdb script**: Streaming aggregation reproduces the duckdb script's filtering and per-group aggregates.
   - Current: aggregation exists only in `/tmp/spectronaut-aggregate.sql` — `EG.Qvalue ≤ 0.01` (non-numeric passes), `CON__/REV__` row drop, `max(PG.Quantity)` and `min(EG.Qvalue)` per (PG × R.Condition × R.Replicate)
   - Target: the TS streaming aggregation applies identical filters and the same per-group aggregates
   - Acceptance: on a small derived precursor fixture, streaming output is row-equivalent to the duckdb script's output — same set of (protein, condition, replicate) group keys, `PG.Quantity` equal within 1e-3, `EG.Qvalue` equal

3. **Pipeline continuity**: The streamed result is shape-identical to today's `parseSpectronautText` output.
   - Current: `parseSpectronautText` yields a wide protein × sample DataFrame with `SEMTYPE` assignment, raw + `log2(...)` columns, `proteomics.source=spectronaut`, `proteomics.preNormalized=true`, and auto-populated `proteomics.groups` when exactly 2 conditions
   - Target: the streaming path produces the identical DataFrame shape and tag set for the reference file (8,328 rows; raw + log2 sample columns; 2-condition auto-groups)
   - Acceptance: after streamed import of the reference file, Annotate → Normalize → Impute → DE → Volcano runs to completion (parity with the manually-verified duckdb path)

4. **Progress + main-thread responsiveness**: Import shows progress and never trips the browser unresponsive guard.
   - Current: no progress feedback on import; large synchronous work has tripped Chrome's "Page Unresponsive" dialog elsewhere (DE)
   - Target: a `DG.TaskBarProgressIndicator` advances by bytes-read; the stream-reader `await` boundaries yield the event loop between chunks
   - Acceptance: importing the 2.6 GB file shows a monotonically advancing progress indicator, no "Page Unresponsive" dialog appears, and the UI stays interactive (tab switching works) during import

5. **Bundled, documented duckdb fallback**: The duckdb aggregation is preserved in-repo as a documented escape hatch.
   - Current: `spectronaut-aggregate.sql` / `.sh` live untracked in `/tmp`
   - Target: both files are committed under the package (tracked location, e.g. `tools/` or `files/`) with README usage documentation and a pointer from the import-failure path
   - Acceptance: the duckdb script exists in the package's git tree with usage docs; an import failure message or README references it as a manual fallback

6. **Automated regression fixture**: A small derived precursor fixture exercises the streaming path in CI.
   - Current: no precursor-level fixture or test exists; the only reference is the uncommittable 2.6 GB file
   - Target: a small synthetic/derived precursor-level TSV fixture plus a test under `src/tests/` asserting streaming-aggregation equivalence (Requirement 2) and output shape (Requirement 3)
   - Acceptance: `grok test` runs a Spectronaut precursor streaming test that passes against the committed small fixture

## Boundaries

**In scope:**
- Streaming precursor-level long-report import through the existing `Proteomics | Import | Spectronaut Report` menu item
- Client-side TypeScript aggregation using browser-native `File.stream()` / `TextDecoderStream` (no new npm deps)
- Filter + aggregate parity with `/tmp/spectronaut-aggregate.sql` (EG.Qvalue ≤ 0.01, CON__/REV__ drop, max PG.Quantity, min EG.Qvalue per group)
- `TaskBarProgressIndicator` and event-loop yielding during import
- Committing the duckdb script into the package as a documented fallback
- A small derived precursor fixture + automated test

**Out of scope:**
- **Pivot/Matrix Report support** — deferred to a later phase; no fixture exists yet and it is a distinct schema problem (already-wide, no `R.Condition`/`R.Replicate`)
- **Peptide-level quantification/analysis** — prior project decision ("doubles parser complexity for marginal gain"); this phase aggregates to protein-group level only, same as the duckdb script
- **Server-side aggregation script (Python/R)** — client-side streaming dominates on deps, UX, and environment requirements; not needed
- **Web Worker offload** — an implementation choice for discuss-phase, not a locked requirement (event-loop yielding is the requirement)
- **Spectronaut spectral-library exports** — not differential-expression inputs
- **Changes to the downstream pipeline** (Normalize/Impute/DE/viewers) — unchanged; this phase only adds an input path

## Constraints

- **Zero new npm dependencies** (standing v1.x principle) — browser-native streaming APIs only; no duckdb-wasm, no parser libraries
- **Bounded memory** — in-memory aggregation state must scale with (#proteins × #samples), not input row count; the implementation must not buffer the whole file or retain all raw rows
- **Filter preservation** — must apply the existing parser's `EG.Qvalue ≤ 0.01` (non-numeric q-values pass) and `CON__/REV__` decoy drop, matching the duckdb script semantics
- **Quantity column** — must accept either `PG.IBAQ` or `PG.Quantity` (already implemented in `spectronaut-parser.ts`, preference order PG.IBAQ → PG.Quantity)
- **Datagrok package conventions** — TypeScript strict, webpack, `//name:` metadata, kebab-case filenames, CRLF, `SEMTYPE.*` from `src/utils/proteomics-types.ts`
- **Fixture logistics** — the 2.6 GB reference file (`~/Downloads/2026-05-13 BP DMD WT.tsv`) is too large to commit; it is the **manual release-gate** fixture, while CI uses the small derived fixture from Requirement 6

## Acceptance Criteria

- [ ] Selecting the 2.6 GB reference TSV via `Proteomics | Import | Spectronaut Report` completes with no V8 string-length or out-of-memory error
- [ ] Resulting DataFrame has 8,328 protein rows and the same 24 sample columns the duckdb path produced
- [ ] Tags set: `proteomics.source=spectronaut`, `proteomics.preNormalized=true`, `proteomics.groups` auto-populated with exactly 2 conditions
- [ ] Streaming aggregation output is row-equivalent to `/tmp/spectronaut-aggregate.sql` output on the small derived fixture (same group keys; `PG.Quantity` within 1e-3; `EG.Qvalue` equal)
- [ ] `EG.Qvalue ≤ 0.01` filter and `CON__/REV__` decoy drop applied identically to the duckdb script
- [ ] A `TaskBarProgressIndicator` advances monotonically during the 2.6 GB import; no Chrome "Page Unresponsive" dialog; tab switching works mid-import
- [ ] Existing pipeline (Annotate → Normalize → Impute → DE → Volcano) runs to completion on the streamed result
- [ ] The duckdb aggregation script is committed under the package with usage documentation
- [ ] A small derived precursor fixture and an automated `grok test` for the streaming path are added under `src/tests/`

## Ambiguity Report

| Dimension          | Score | Min  | Status | Notes                                                              |
|--------------------|-------|------|--------|--------------------------------------------------------------------|
| Goal Clarity       | 0.90  | 0.75 | ✓      | Scope locked: precursor-level only, streaming client-side path     |
| Boundary Clarity   | 0.90  | 0.70 | ✓      | Pivot/Matrix + peptide-level + server-script explicitly deferred   |
| Constraint Clarity | 0.78  | 0.65 | ✓      | Streaming bounds memory; CPU/event-loop yielding is the live limit  |
| Acceptance Criteria| 0.85  | 0.70 | ✓      | E2E on 2.6 GB + duckdb-equivalence on small fixture, pass/fail      |
| **Ambiguity**      | 0.134 | ≤0.20| ✓      |                                                                    |

Status: ✓ = met minimum, ⚠ = below minimum (planner treats as assumption)

## Interview Log

| Round | Perspective                  | Question summary                                  | Decision locked                                                       |
|-------|------------------------------|---------------------------------------------------|-----------------------------------------------------------------------|
| 1     | Researcher / Boundary / Simplifier | Scope split? Automation bar? Aggregation runtime? | Precursor-level only; pivot deferred; runtime question superseded by streaming-TS |
| 1     | (technical)                  | Can TS stream-aggregate to avoid browser memory?  | Yes — output bounded by proteins×samples; native streams, zero new deps; becomes the recommended path |
| 2     | Failure Analyst / Seed Closer | Pass/fail acceptance bar? Fate of duckdb script?  | End-to-end on the 2.6 GB fixture (manual gate) + small-fixture equivalence (CI); keep duckdb script as a bundled, documented fallback |

---

*Phase: 12-spectronaut-input-coverage*
*Spec created: 2026-05-15*
*Next step: /gsd-discuss-phase 12 — implementation decisions (streaming chunk strategy, worker vs. yield, fixture derivation, etc.)*
