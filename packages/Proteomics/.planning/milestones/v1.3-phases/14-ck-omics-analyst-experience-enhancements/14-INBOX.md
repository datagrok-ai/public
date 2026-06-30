---
source: phase-13 UAT (2026-05-28)
purpose: gaps forwarded from Phase 13 verify-work session — fold into 14-CONTEXT.md / 14-SPEC.md when /gsd:discuss-phase 14 runs
---

# Forwarded Gaps from Phase 13 UAT

Phase 13 (`13-ck-omics-volcano-and-enrichment-parity`) UAT completed 2026-05-28 with 4/5 tests pass and 1 test issue + 3 adjacent issues on passing tests. User chose to route all 4 gaps into Phase 14 rather than open a Phase 13 gap-closure plan.

See `.planning/phases/13-ck-omics-volcano-and-enrichment-parity/13-HUMAN-UAT.md` `## Gaps` section for full structured YAML (truth/status/reason/severity/test/missing).

## Summary table

| # | Gap | Severity | From Test |
|---|---|---|---|
| 1 | Volcano lacks CK-omics visual parity: generic "Volcano" title, raw axis labels (negLog10P/log2FC), no default labels for top hits, generic "direction" coloring instead of "Enriched in WT/DMD" with magenta/cyan, no default legend with counts, no protein search input | major | 1 |
| 2 | Volcano Options dialog shows defaults on re-open instead of current viewer state | minor | 3 |
| 3 | Switching Color → Location is slow and shows no progress indicator | minor | 3 |
| 4 | Multi-contrast Filters viewer auto-includes the Spectronaut `Flags` column (Valid/significant) despite `DG.Viewer.filters(df, {columnNames: ['Comparison (group1/group2)']})` scoping | minor | 4 |

## What passed in Phase 13 (do NOT regress when fixing the above)

- Volcano point distribution and orientation match CK-omics with no manual Comparison flip needed (Test 1, partial)
- Directional enrichment layout: 4 viewers (Up/Down × dot/bar), Direction column, WP source, cross-link highlighting on volcano (Test 2)
- Significance-metric and color toggles stay synchronized in the live volcano (Test 3)
- Multi-contrast Candidates auto-docks a Filters viewer (Test 4 — gap is which columns it contains, not that it's docked)
- UniProt subcellular-location fetch + cross-session cache works and looks biologically sane (Test 5)
