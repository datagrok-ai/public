# Phase 14: CK-omics Analyst-Experience Enhancements - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-05-28
**Phase:** 14-ck-omics-analyst-experience-enhancements
**Areas discussed:** Scope packaging, Volcano polish (G1/G2/G3/R2), Gene-label resolution (R1), Per-protein panel + correlation scatter + smart pathway filter (R3/R4/R5)

---

## Scope packaging

### Q1.1 — Packaging the 5 R + 4 gaps

| Option | Description | Selected |
|--------|-------------|----------|
| One phase, all in | Phase 14 covers R1–R5 plus all 4 forwarded gaps; volcano polish folds into R2. | ✓ |
| Split into 14a (UX polish) + 14b (interpretation) | 14a = G1+G2+G3+G4+R2; 14b = R1+R3+R4+R5. Two smaller phases. | |
| Split by source | Phase 14 = ROADMAP R1–R5 only; spin a separate Phase 15 for the 4 gaps. | |

**User's choice:** One phase, all in
**Notes:** Risks a fat phase but keeps the work coherent; volcano polish and R2 share the same viewer.

### Q1.2 — Must-haves for BP DMD/WT client deliverable

| Option | Description | Selected |
|--------|-------------|----------|
| G1 — Volcano visual parity | Title, axes, labels, semantics, legend, search | ✓ |
| R1 — Predicted gene-label resolution | ENSRNOG etc. → names + provenance | ✓ |
| R2 — Live filter-aware counters on volcano | Counters update as legend/filter toggles | ✓ |
| R3 — Click-point per-protein group-quantity panel | Adds per-group means to the UniProt panel | |

**User's choice:** G1, R1, R2
**Notes:** Nice-to-haves (R3, R4, R5, G2, G3, G4) layer in time-permitting; planner sequences MUST-haves first.

---

## Volcano polish (G1 + G2 + G3 + R2)

### Q2.1 — Default protein labels

| Option | Description | Selected |
|--------|-------------|----------|
| Top-N by significance | Auto-label top ~15–20 by adj.p-value (or active metric); recompute on metric/filter change. Matches CK-omics. | ✓ |
| Top-N by \|log2FC\| | Top ~15 by absolute fold-change. Catches biological outliers; less aligned with CK-omics. | |
| Top-N split: top of each direction | ~10 up + ~10 down. Balanced labeling. | |
| User must opt in (no defaults) | No default labels; analyst types into search. | |

**User's choice:** Top-N by significance

### Q2.2 — Direction-label semantics source

| Option | Description | Selected |
|--------|-------------|----------|
| From `proteomics.groups` tag | Read group1/group2 names from GroupAssignment JSON. Locked magenta/cyan/gray. Single code path. | ✓ |
| From Comparison column string in Candidates files | Parse 'DMD / WT'. Works for Candidates only — two code paths. | |
| Configurable in Volcano Options dialog | Defaults from group names; analyst can rename/recolor. More flexibility, less parity. | |

**User's choice:** From proteomics.groups tag

### Q2.3 — Protein search/highlight UX

| Option | Description | Selected |
|--------|-------------|----------|
| Inline below volcano (CK-omics-style) | Always-visible textarea, multi-line + comma-separated. | |
| Volcano Options dialog only | Lives inside the dialog. Less discoverable, cleaner volcano area. | |
| Property panel input | Datagrok-native pattern — viewer property. | |
| Floating search box / keyboard shortcut | Press '/' or '⌘F' overlay. | |
| **(Freeform) Native Datagrok Filters viewer with text-search box** | User-proposed: dock a Filters viewer scoped to Gene name + Protein ID. String-column filters have built-in search; wire matches to `df.selection.set(...)` so points highlight in place instead of being filtered out. Unifies with G4 (extend the existing Comparison Filter viewer rather than adding a new one). | ✓ |

**User's choice:** Filters viewer (freeform redirect)
**Notes:** Better than any of my prepared options — reuses platform muscle memory, gets the search box for free, no custom component to build, and unifies with G4 by extending the existing Comparison Filter viewer.

### Q2.4 — Live counter location and recomputation trigger

| Option | Description | Selected |
|--------|-------------|----------|
| Floating overlay on volcano + filter/selection events | Small text overlay anchored bottom-right; recompute on df.filter / df.selection / viewer-property changes. | ✓ |
| Datagrok status bar at table-view bottom | Custom statbar entry; discoverable for Datagrok-natives, less prominent. | |
| Inside legend block | Per-category counts in legend ('Enriched in WT (n=339)'). Matches CK-omics legend format. | |
| Both legend counts + overlay total | Splits info. | |

**User's choice:** Floating overlay on volcano + filter/selection events

### Q2.5 — Continue probing or move on?

| Option | Description | Selected |
|--------|-------------|----------|
| Move on — remaining volcano items are mechanical | G2/G3/title/axis labels become Claude's discretion. | ✓ |
| One more pass — specify G3 progress UX explicitly | Lock progress indicator pattern + caching strategy. | |
| One more pass — specify G2 dialog state strategy | Lock state source-of-truth before planner picks. | |

**User's choice:** Move on

---

## Gene-label resolution (R1)

### Q3.1 — Resolution scope (when and on what scale)

| Option | Description | Selected |
|--------|-------------|----------|
| Eager on import, all rows | Detect-and-resolve runs as part of the parser. Display Name column populated before any viewer renders. | ✓ |
| Lazy on first visualization, all rows | Resolution kicks in on first volcano/heatmap open. | |
| Lazy on visible/highlighted rows only | Only resolve viewport/selected rows. Smallest network footprint. | |

**User's choice:** Eager on import, all rows

### Q3.2 — Provenance marker placement

| Option | Description | Selected |
|--------|-------------|----------|
| Inline in Display Name string | E.g. 'Myh7*' / 'Predicted-LOC123†'. Travels everywhere; original ID kept in Source ID column. | ✓ |
| Separate Provenance column | Categorical column with 'reviewed'/'grouped'/'predicted'. Filter/color affordance. | |
| Both — inline marker + provenance column | Belt-and-suspenders. | |

**User's choice:** Inline in Display Name string

### Q3.3 — Resolver backend

| Option | Description | Selected |
|--------|-------------|----------|
| Ensembl REST POST /lookup/id in batches of ~1000 | Group by species → POST per species → cross-session cache via grok.userSettings. Three-level fallback name→desc→raw. CORS via grok.dapi.fetchProxy(). | ✓ |
| g:Profiler /convert endpoint (reuses existing client) | Single REST stack to maintain. Tradeoff: no descriptions, no provenance signal. | |
| Cached static reference file shipped with package | Zero network latency. Tradeoff: stale, big bundle, no help for novel IDs. | |

**User's choice:** Ensembl REST POST /lookup/id in batches of ~1000

### Q3.4 — Duplicate-description handling

| Option | Description | Selected |
|--------|-------------|----------|
| Disambiguate by appending source ID | 'Myh7 (ENSRNOG00000001)' / 'Myh7 (ENSRNOG00000099)'. Single warning toast. | ✓ |
| Group rows on import (CK-omics-style * grouped marker) | Roll duplicates into one row with aggregation. Matches CK-omics literally; lossy. | |
| Warn only, leave names ambiguous | Log warning; render all as 'Myh7'. Hover shows source ID. | |

**User's choice:** Disambiguate by appending source ID

---

## Per-protein panel + correlation scatter + smart pathway filter

### Q4.1 — Per-group quantity visualization in click-point panel

| Option | Description | Selected |
|--------|-------------|----------|
| Compact bar chart, one bar per group with mean ± SD whiskers | Two/N bars; numeric mean+SD below. Reuses ui.canvas / svg helper. | ✓ |
| Numeric table only | Group / Mean / SD / N. Cheapest, no gut-check. | |
| Sparkline-style strip with per-sample dots | Reveals within-group spread; harder to glance. | |
| Reuse the existing QC waterfall pattern | Plug in qc-computations.ts. Consistent; visually heavy for a panel. | |

**User's choice:** Compact bar chart, mean ± SD whiskers

### Q4.2 — Correlation scatter menu location and viewer type

| Option | Description | Selected |
|--------|-------------|----------|
| Proteomics \| Visualize \| Group-Mean Correlation… — native scatter | Sibling of Volcano/Heatmap/PCA. createGroupMeanCorrelation() factory. Native scatter + derived mean columns. | ✓ |
| Proteomics \| QC Dashboard — add as 8th panel | Fold into the existing dashboard. ROADMAP explicitly says 'distinct from QC sample×sample heatmap'. | |
| Proteomics \| Analyze \| Group-Mean Correlation… | Goes under Analyze. Mirrors where Enrichment lives. | |

**User's choice:** Visualize | Group-Mean Correlation… — native scatter

### Q4.3 — Smart pathway filtering rule source

| Option | Description | Selected |
|--------|-------------|----------|
| Port CK-omics rules verbatim | Locked contract — same parent-detection heuristic, same per-source FDR cap. Researcher reads the function during research. | ✓ |
| Hardcoded simple rule — drop GO parents + cap top 20 per source | Simpler; may diverge from CK-omics result set. | |
| Configurable cap with sensible default | Dialog input. Added dialog complexity for an edge case. | |

**User's choice:** Port CK-omics rules verbatim

### Q4.4 — R5 default behavior

| Option | Description | Selected |
|--------|-------------|----------|
| Default-on, with a dialog toggle to disable | Mirrors D-12 (WikiPathways default-on). One checkbox added. | ✓ |
| Default-off, opt-in toggle | Conservative — don't silently mutate. | |
| Always-on, no toggle | Simplest; methodology check requires code dive. | |

**User's choice:** Default-on with toggle

---

## Claude's Discretion

- Volcano title source synthesis (proteomics.de_method + group names)
- Axis label rewriting on metric toggle (negLog10P → -Log10(p-value)/-Log10(Q-value); log2FC → Log2 Fold Change)
- G2 dialog state preload mechanism (viewer.getOptions vs df.tags vs input cache)
- G3 progress indicator widget (DG.TaskBarProgressIndicator vs viewer host overlay)
- Ensembl batch size cap fine-tuning per species
- Whether Display Name replaces or augments existing Gene name column
- `*` vs `†` specific rules per CK-omics reference function
- Per-bar styling in the per-group quantity panel

## Deferred Ideas

- Custom volcano JsViewer (already PROJECT.md out-of-scope)
- Per-protein waterfall chart in UniProt panel (vs simpler bars in D-11)
- Configurable per-source pathway cap in Enrichment dialog (vs D-14 binary toggle)
- Correlation viewer as QC dashboard panel (vs D-12 standalone menu entry)
- Provenance-driven color/filter in viewers (vs D-08 inline-marker only)
- Methodology-check / power-user raw-enrichment flow (vs D-14 default-on toggle)
