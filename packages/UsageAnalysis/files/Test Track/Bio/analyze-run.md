# Bio Analyze — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open Bio > Analyze > Sequence Space... on FASTA | PASS | Dialog opened; Table=sample_FASTA, Column=Sequence, Method=UMAP, Similarity=Hamming |
| 2 | Click OK — viewer opens | PASS | Scatter plot "Sequence space" added; X_1/Embed_Y_1/Cluster(DBSCAN) columns created |
| 3 | Open Bio > Analyze > Activity Cliffs... on FASTA | PARTIAL | Menu click threw NullError; opened via `func.prepare().edit()` workaround; dialog shows correctly |
| 4 | Click OK — viewer opens | PASS | "Activity cliffs" viewer appeared: 2 CLIFFS detected |
| 5 | Open Bio > Analyze > Composition on FASTA | PASS | WebLogo viewer added to view |
| 6 | Click Gear on Composition viewer — check settings | PASS | Context panel shows Data section (Sequence, Value, Aggr Type, Skip Empty, Shrink Tail), Layout, Behavior, Style |
| 7 | Run Sequence Space again with changed params (t-SNE) | PASS | Method changed to t-SNE; new Embed_X_3/Embed_Y_3 columns added |
| 8 | Repeat all 3 analyses on sample_HELM.csv (540 rows) | PASS | Sequence Space ✓, Activity Cliffs ✓ (14529 cliffs), Composition/WebLogo ✓ |
| 9 | Repeat all 3 analyses on sample_MSA.csv (540 rows) | PASS | Sequence Space ✓, Activity Cliffs ✓, Composition/WebLogo ✓ |

## Summary

All three Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three sample datasets (FASTA, HELM, MSA). The only issue was that Bio > Analyze > Activity Cliffs fails silently when triggered via the menu on FASTA — it throws "molecules: Value not defined / NullError: method not found 'name' on null" and no dialog appears. The workaround was to call `func.prepare().edit()` which correctly showed the dialog with the Sequence column pre-selected.

## Retrospective

### What worked well
- Sequence Space dialog opens correctly from the menu with sensible defaults
- All 3 sample files have Macromolecule semantic type auto-detected correctly
- Composition (WebLogo) viewer opens directly without a dialog, adds immediately
- Gear icon on WebLogo opens full settings panel with Data/Layout/Behavior/Style sections
- Activity Cliffs produces correct output once dialog is shown

### What did not work
- Bio > Analyze > Activity Cliffs menu item fails silently on first call — NullError when trying to detect the molecule column from context; workaround: call via `grok.functions.eval('Bio:activityCliffs').prepare(...).edit()`
- `NullError: method not found: 'name' on null` appears in console for each Activity Cliffs menu click attempt

### Suggestions for the platform
- Activity Cliffs menu handler should detect the Macromolecule column from the current table and pre-populate the dialog instead of failing silently
- Console errors from Bio package functions should show user-visible toasts rather than silent null errors

### Suggestions for the scenario
- Clarify that Activity Cliffs requires an activity (numeric) column — the default picked was "Length" not "Activity"; scenario should note to set Activities=Activity
- Step 3 "Run once more with arbitrary changed parameters" should specify one function to rerun, not all three
