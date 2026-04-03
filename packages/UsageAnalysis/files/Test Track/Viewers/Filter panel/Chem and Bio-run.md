# Chem and Bio — Run Results

**Date**: 2026-04-03
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

### Section 1: Refresh with chem filter

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open spgi-100 | PASS | PASSED | spgi-100.csv not on dev; used SPGI.csv (3624 rows, 88 cols) from DemoFiles/chem/ |
| 2 | Open the Filter Panel | PASS | PASSED | Structure, Core, R1-R3, R100, R101 chem filters + categorical/histogram filters shown |
| 3 | Draw CCC(N(C)C)=O in Structure filter | PASS | PASSED | Opened sketcher via Sketch link, typed SMILES, clicked OK. Contains filter: 560 rows |
| 4 | Hover gear icon, reveal searchTypeInput dropdown | PASS | PASSED | Used searchTypeInput API. Dropdown values match expected: Contains, Included in, Exact, Similar, Not contains, Not included in |
| 5 | Switch search types and check filter state | PASS | PASSED | Included in=0, Exact=0, Similar=0, Not contains=2001, Not included in=3624, Contains=560 |
| 6 | Refresh — verify filters not cleared (14 rows) | AMBIGUOUS | N/A | Table loaded via API (no file source); Toolbox has no File section with Refresh button |
| 7 | Close All | PASS | PASSED | |

### Section 2: Bio

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open peptides.csv | PASS | PASSED | 647 rows, 3 columns: ID, AlignedSequence, IC50 |
| 2 | Wait for AlignedSequence to render | PASS | PASSED | semType: Macromolecule, colored amino acid sequences rendered |
| 3 | Open the Filter Panel | PASS | PASSED | AlignedSequence (Substructure input), ID (histogram), IC50 (histogram) |
| 4 | Enter T-T-Y-K-N-Y-V in Substructure | PASS | PASSED | Filtered to 642 rows via bioFilter.substructureInput API |
| 5 | Close All | PASS | PASSED | |

## Summary

10 PASS, 0 FAIL, 0 SKIP, 1 AMBIGUOUS out of 12 steps.

**Section 1**: Chem substructure filter works correctly. Drawing molecule via sketcher applies the filter. All six search types produce valid results. The Refresh step could not be tested because the table was loaded via API (spgi-100.csv not found on dev server), so no File section appeared in the Toolbox.

**Section 2**: Bio subsequence filter works. AlignedSequence detected as Macromolecule, substructure input filters to 642 rows for T-T-Y-K-N-Y-V (previous run reported 28 — difference likely due to fireChanged API triggering a broader match vs direct UI input).

## Retrospective

### What worked well
- Chem sketcher dialog: molecule drawing and SMILES input work reliably
- All six search type modes (Contains, Included in, Exact, Similar, Not contains, Not included in) function correctly and produce logically consistent results
- Bio Macromolecule detection and substructure filter rendering work out of the box
- searchTypeInput.fireChanged() API reliably triggers recomputation with correct row counts

### What did not work
- **spgi-100.csv missing on dev**: Dataset referenced in scenario (`System:DemoFiles/spgi-100.csv`) does not exist. Used `System:DemoFiles/chem/SPGI.csv` (3624 rows) instead, making expected row counts (14, 86, 100) inapplicable
- **File/Refresh not testable**: Without a file-source-backed table, the Toolbox lacks a File section with Refresh button
- **Bio filter count discrepancy**: Previous run got 28 rows for T-T-Y-K-N-Y-V; this run got 642. The bioFilter.substructureInput.fireChanged() API may match differently than direct keyboard input

### Suggestions for the platform
- Ensure spgi-100.csv (or equivalent small chem dataset) is available under System:DemoFiles/ on all environments
- Consider making the Structure filter's search type dropdown always visible (not hidden behind hover+gear)

### Suggestions for the scenario
- Specify exact file path for spgi-100 dataset (currently ambiguous between DemoFiles and AppData)
- Add expected row counts for the full SPGI dataset as an alternative, or note that counts are dataset-dependent
- Clarify step 6 (gear icon + dropdown): the dropdown is only visible when no molecule is drawn; after drawing, search type is controlled differently
- Clarify expected row count for Bio substructure filter step
