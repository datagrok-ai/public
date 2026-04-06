# Similarity Search — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open smiles.csv | PASS | Opened via JS API, 1000 rows, 20 cols |
| 2 | Initiate Similarity Search | PASS | Added Chem Similarity Search viewer via JS API. Panel shows "Most similar structures" with Tanimoto/Morgan, reference molecule and 11 similar structures with scores (0.63 to 0.29). |
| 3 | Access properties (gear icon) | PASS | Viewer properties available: fingerprint, limit, distanceMetric, size, moleculeProperties, cutoff, followCurrentRow |
| 4 | Test property modifications | PASS | Changed fingerprint Morgan→MACCS (no error), limit 12→5 (no error), cutoff set to 1 — only reference molecule remained (correct behavior). All changes applied without errors. |

## Summary

All 4 steps passed. Similarity Search viewer displays correctly with molecule similarity scores. Property modifications (fingerprint, limit, cutoff) all work without errors. Cutoff=1 correctly shows only the reference molecule.

## Retrospective

### What worked well
- Similarity search computed and displayed results immediately
- All property changes applied without errors or crashes
- Cutoff=1 correctly filtered to only the reference molecule
- Multiple fingerprint types (Morgan, MACCS) work

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- None

### Suggestions for the scenario
- Add expected similarity score ranges for verification
- Specify which molecule to use as reference for reproducible results
- Add step to verify "size" property visual changes (small/normal/large)
