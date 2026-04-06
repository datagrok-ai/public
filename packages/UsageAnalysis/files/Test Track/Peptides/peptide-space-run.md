# Peptide Space — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Launch SAR from Bio->Analyze->SAR | PASS | PASSED | Opened Bio ribbon → Analyze submenu → SAR...; "Analyze Peptides" dialog appeared with Sequence=AlignedSequence, Activity=IC50, Scaling=none |
| 2 | Wait for calculation results | PASS | PASSED | Clicked OK; SAR view loaded with MCL network graph (orange nodes, 1 cluster of 646), Logo Summary Table, Sequence Variability Map, Most Potent Residues; Columns: 26 |
| 3 | Open settings using wrench button | PASS | PASSED | Clicked `.fa-wrench` icon; "Peptides settings" dialog opened showing General (Activity scaling), MCL (Similarity Threshold=70, Min Cluster Size=5), Viewers section |
| 4 | Adjust arbitrary parameters and click "ok" | PASS | PASSED | Changed Activity scaling: none→lg, Similarity Threshold: 70→55, Min Cluster Size: 5→3 via native setters; clicked OK |
| 5 | MCL Viewer should output different results | PASS | PASSED | MCL graph recalculated: nodes changed from orange to blue, different oval/disc layout; new Cluster size (MCL) and Connectivity (MCL) columns visible in selection panel |

## Summary

All 5 steps passed. SAR was successfully launched via Bio→Analyze→SAR (alternative to the "Launch SAR" button in the context panel). The MCL viewer correctly reflected the new parameters after clicking OK — node color and layout both changed, confirming recalculation.

## Retrospective

### What worked well
- Bio→Analyze→SAR is a reliable alternative path to launch SAR analysis
- The "Analyze Peptides" dialog auto-populated Sequence and Activity fields from the open dataset
- Wrench button (`.fa-wrench`) opened Peptides settings reliably after SAR was launched
- Native setter pattern for inputs (Similarity Threshold, Min Cluster Size) worked correctly
- MCL visual output clearly changed color (orange→blue) confirming parameter effect

### What did not work
- MCL computation took ~20–25 seconds total; no progress indicator visible during recalculation

### Suggestions for the platform
- Add a visible progress bar or spinner during MCL recalculation (currently the graph area is blank until complete)
- Expose `data-testid` attributes on the MCL graph nodes for automated cluster count verification

### Suggestions for the scenario
- Step 1 should clarify that the dataset (peptides.csv) must be open with AlignedSequence column selected before using Bio→Analyze→SAR
- Step 4 should specify which parameters to change for repeatable testing (e.g. "change Similarity Threshold to 55")
- Step 5 should clarify what "different results" means — e.g., node count, cluster count, or node color
