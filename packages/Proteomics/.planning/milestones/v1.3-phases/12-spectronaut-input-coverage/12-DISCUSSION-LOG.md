# Phase 12: Spectronaut Input Coverage - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-05-15
**Phase:** 12-spectronaut-input-coverage
**Areas discussed:** Format routing, Execution model, Test fixture, Duckdb home, Equivalence oracle

---

## Format routing

| Option | Description | Selected |
|--------|-------------|----------|
| Header-sniff, then branch | Read header line; precursor columns → streaming path, else existing path | ✓ |
| Always stream | Route every import through streaming aggregator | |
| File-size threshold | Stream only above a size threshold | |

**User's choice:** Header-sniff, then branch
**Notes:** Keeps proven small-file path unchanged; size is not a reliable proxy for report shape.

---

## Execution model

| Option | Description | Selected |
|--------|-------------|----------|
| Main thread + await-yield | Async chunk loop; reader await yields event loop; bytes-read progress | ✓ |
| Web Worker | Offload parse+aggregate to a Worker | |
| Decide in plan-phase | Defer to planner/spike | |

**User's choice:** Main thread + await-yield
**Notes:** Same responsiveness pattern that fixed the DE Page-Unresponsive issue; avoids webpack worker config.

---

## Test fixture

| Option | Description | Selected |
|--------|-------------|----------|
| Synthetic generator extending makeLongFormatTsv | Emit precursor rows + decoys + sub-threshold qvalues | ✓ |
| Derived slice of the real 2.6GB file | Commit a head -N slice | |
| Both | Synthetic + curated real slice | |

**User's choice:** Synthetic generator extending makeLongFormatTsv
**Notes:** Reproducible, no client-data licensing, exercises every filter branch.

---

## Duckdb home

| Option | Description | Selected |
|--------|-------------|----------|
| tools/ + README + import-error hint | New tools/ dir, README docs, failure-path pointer | ✓ |
| files/ (runtime-served) | Under files/ with demo data | |
| docs/ only | Docs folder, no in-product pointer | |

**User's choice:** tools/ + README + import-error hint
**Notes:** Dev/CLI tool, semantically wrong under runtime-served files/.

---

## Equivalence oracle

| Option | Description | Selected |
|--------|-------------|----------|
| Committed golden file | Run duckdb once over synthetic fixture, commit output; test asserts equality | ✓ |
| Hand-computed expected values | Encode expected rows in the test | |
| Port duckdb logic to a TS test helper | Reimplement SQL as TS reference | |

**User's choice:** Committed golden file
**Notes:** Pins behavior to real duckdb logic; needs a documented regeneration command.

---

## Claude's Discretion

- Stream chunk size, cross-chunk line-buffer carry-over, progress-update cadence, exact failure-hint wording.

## Deferred Ideas

- Pivot/Matrix Report support — own future phase (no fixture, distinct schema).
- Peptide-level quantification — standing deferral (parser-complexity vs. gain).
