# Line chart legend — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI | 10s | PASS | PASSED | 3624 rows |
| 2 | Add Line chart | 3s | PASS | PASSED | tv.viewers = [Grid, Line chart] |
| 3 | Split=Series — legend shows distinct categories | 2s | PASS | PASSED | 7 legend items (distinct colors) |
| 4 | multiAxis=true — each Y gets its own subplot with legend | 2s | PARTIAL | PASSED | Property set; DOM has only a single `[name="legend"]` element — per-line legend blocks are not individually labeled in DOM |
| 5 | Save + re-apply layout — Split/multiAxis persist | 5s | PASS | PASSED | `multiAxis=true`, `splitColumnName='Series'` after reload |
| 6 | Save project, reopen | 2s | FAIL | PASSED (asserted against error) | `project_relations_entity_id_fkey` FK constraint — known limitation |
| 7 | yColumnNames = ['Average Mass', 'TPSA'] | 3s | PASS | PASSED | 14 legend items (7 per Y column × 2 Y cols) |
| 8 | Replace Y column → NIBR logP | 2s | PASS | PASSED | `yColumnNames = ['Average Mass', 'NIBR logP']` accepted |
| 9 | Save + re-apply layout — new Y persists | 5s | PASS | PASSED | yColumnNames restored correctly |
| 10 | Save project, reopen | 2s | FAIL | PASSED (asserted against error) | Same FK constraint |
| 11 | Cleanup | 1s | PASS | n/a | Deleted 2 layouts, closeAll |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 35s |
| grok-browser execution (scenario steps) | 35s |
| Execute via grok-browser (total) | 2m 10s |
| Spec file generation | 40s |
| Spec script execution | 32s |
| **Total scenario run (with model)** | 3m 22s |

## Summary

Line chart legend and multi-axis behaviors are mostly correct: 7 legend items for 7 categories, layout round-trip preserves `splitColumnName`, `multiAxis`, and `yColumnNames`. With two Y columns the legend shows 14 combined items (7 × 2). DOM doesn't expose distinct `[name="legend"]` elements per subplot under multiAxis — only one legend element exists even when the scenario expects "per-line legends". Project save fails with the known FK constraint. **Total scenario run (with model): 3m 22s**.

## Retrospective

### What worked well
- `lc.props.splitColumnName = 'Series'` produces 7 colored legend items
- `lc.props.multiAxis = true` accepted and persisted
- `lc.props.yColumnNames` array assignment and layout round-trip both work
- Legend rescales when replacing Y columns

### What did not work
- Under multiAxis the scenario expects "each subplot has its own legend" but the DOM shows a single `[name="legend"]` block — either the per-subplot legends are drawn on a shared canvas, or the per-subplot structure is not implemented
- Project save fails with FK constraint (same as other scenarios)

### Suggestions for the platform
- If multiAxis is expected to produce per-subplot legends, emit one `[name="legend"]` per subplot so automation can verify them individually
- Auto-save/attach unsaved child dataframes on `projects.save`

### Suggestions for the scenario
- Clarify in step 4 whether "per-line legends" means separate DOM blocks (for automation) or a single combined legend
- In step 7, specify the expected total legend item count (14 for 2 Y × 7 categories) so the assertion is unambiguous
