# Phase 13: CK-omics Volcano and Enrichment Parity - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-05-16
**Phase:** 13-ck-omics-volcano-and-enrichment-parity
**Areas discussed:** UniProt location fetch + caching, Volcano UI: location colour + Q/P, Candidates contrast picker + flip, Up/down enrichment presentation

---

## UniProt location fetch + caching (R1)

| Option | Description | Selected |
|--------|-------------|----------|
| Batched stream via fetchProxy | Chunked accession stream query through grok.dapi.fetchProxy; few requests, CORS-safe | ✓ |
| Per-accession via fetchProxy | One request per accession through fetchProxy; simplest reuse, slow | |
| Keep raw fetch(), batched | Batched stream but raw fetch(); keeps CORS risk | |

| Option | Description | Selected |
|--------|-------------|----------|
| Session in-memory, shared module | Map keyed by accession, session-only; persistent stays a todo | |
| Persistent cross-session cache | IndexedDB / grok.userSettings; closes the cache-uniprot todo fully | ✓ |
| No cache this phase | Fetch every run | |

| Option | Description | Selected |
|--------|-------------|----------|
| Keep reviewed-entry-by-gene-name fallback | Look up reviewed entry by gene name for unreviewed/empty | ✓ |
| Drop it — mark Unknown | Unreviewed/empty → Unknown directly | |

**User's choice:** Batched stream via fetchProxy; persistent cross-session cache; keep the fallback.
**Notes:** Persistent cache fully folds todo `2026-03-03-cache-uniprot-protein-content`; shared module to be reused by `uniprot-panel.ts`.

---

## Volcano UI: location colour + Q/P (R1, R6)

| Option | Description | Selected |
|--------|-------------|----------|
| One volcano, toggle colour column | Switch colour between significance and Subcellular Location via viewer property | ✓ |
| Separate location-coloured volcano | Second viewer/menu item, CK-omics-faithful, more wiring | |
| Creation-dialog choice only | Pick colouring mode once, no live toggle | |

| Option | Description | Selected |
|--------|-------------|----------|
| Creation-dialog option, drives axis + significance | Metric chosen at creation; default adj.p | |
| Live viewer property toggle | Switch metric post-creation; recomputes axis + significance + threshold lines live | ✓ |
| Y-axis only, significance stays adj.p | Toggle changes only Y scale | |

**User's choice:** One volcano with colour-column toggle; live viewer-property Q/P toggle.
**Notes:** Live toggle must keep classification + threshold lines in sync on every switch. 11-category palette locked from CK-omics README. Default metric adj.p-value (= client Q).

---

## Candidates contrast picker + flip (R3, R4)

| Option | Description | Selected |
|--------|-------------|----------|
| Import-time dialog, pick numerator/denominator | CK-omics-style modal dropdowns | |
| Import all + filter column | Keep all comparisons; selector/filter on the view | ✓ (refined) |
| One volcano per contrast | Auto-split into N views | |

| Option | Description | Selected |
|--------|-------------|----------|
| Numerator-positive, auto-flip reverse rows | positive log2FC = chosen-numerator-enriched; invert sign + swap quantity cols | ✓ |
| Honor file's declared direction as-is | Pick contrast, never flip | |

| Option | Description | Selected |
|--------|-------------|----------|
| Candidates-only this phase | Report-DE alphabetical defect stays parked | |
| Also unpark the report-DE fix | Fix differential-expression.ts to honor declared/intended contrast | ✓ |

**User's choice:** Import all + **native Datagrok Filter viewer** on `Comparison` (free-text refinement of option 2); numerator-positive auto-flip; also unpark the report-DE fix.
**Notes:** No explicit numerator/denominator modal → canonical orientation taken from the declared `group1/group2` label; parser normalizes sign to that label. Tradeoff explicitly confirmed by user: lose the explicit "view the mirror" control in exchange for filter-widget simplicity. D-09 is the explicit prompt that overrides the memory's "do not fix unprompted".

---

## Up/down enrichment presentation (R2)

| Option | Description | Selected |
|--------|-------------|----------|
| One DataFrame + Direction column | Two queries merged; reuse Phase-9 viewers + cross-link | ✓ |
| Two separate enrichment table views | CK-omics-faithful; doubles viewer wiring | |

| Option | Description | Selected |
|--------|-------------|----------|
| Split Up/Down charts side-by-side | Two dot/bar charts, Up vs Down | ✓ |
| Single chart, Direction filter | One chart, toggle direction | |
| Single chart, Direction as colour/facet | One chart, Up/Down encoded | |

**User's choice:** One DataFrame + Direction column; split Up/Down charts side-by-side.
**Notes:** Background stays all-detected-proteins for both directional queries (pre-stated as decided, not contested).

---

## Claude's Discretion

- Accession chunk size; cache storage mechanism (IndexedDB vs grok.userSettings) and invalidation; progress cadence.
- Colour-toggle switch wording/placement; how the live Q/P toggle is wired into the ScatterPlotViewer property panel.
- Dock arrangement/sizing for side-by-side Up/Down enrichment charts and the `Comparison` Filter viewer.
- Whether directional gene lists are also surfaced as artifacts.
- R5/WP: pre-decided (add `'WP'` to g:Profiler sources, default-on) — no discussion turn.

## Deferred Ideas

- Explicit "view the mirror of a declared contrast" control — traded away in D-07.
- Surfacing directional gene-list exports (`*_up_genes_*.txt` / `*_down_genes_*.txt`) — not required for in-product parity.
- Reviewed-not-folded todo: `2026-03-03-expand-de-dialog-with…comparison-picker` (related but DE-dialog scope, distinct from R3/R4); ~19 other matches were keyword noise.
