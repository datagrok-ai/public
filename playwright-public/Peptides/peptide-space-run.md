# Peptide Space — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Launch SAR from Bio > Analyze > SAR | PASS | 3s | PASSED | Dialog "Analyze Peptides" appeared with AlignedSequence and IC50 pre-selected |
| 2 | Wait for calculation results | PASS | 15s | PASSED | MCL, Most Potent Residues, Sequence Variability Map viewers all created. Cluster (MCL) column added. |
| 3 | Open settings using wrench button | AMBIGUOUS | - | SKIP | No wrench button found in MCL viewer. Tried gear icon (opens standard properties), icon-edit and icon-plus (both open "Add New Column" dialog). |
| 4 | Adjust arbitrary parameters and click OK | SKIP | - | SKIP | Cannot adjust MCL parameters — no settings dialog accessible |
| 5 | MCL Viewer should output different results | SKIP | - | SKIP | Prerequisite (step 3-4) not completed. MCL viewer shows 2 clusters with initial parameters. |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 25s |
| Spec file generation | 3s |
| Spec script execution | 1.3m |

## Summary

SAR analysis launches correctly via Bio > Analyze > SAR and produces MCL, Most Potent Residues, and Sequence Variability Map viewers. MCL clustering generates 2 clusters with EmbedX/Y, Cluster, Cluster size, and Connectivity columns. However, the "wrench button" for adjusting MCL parameters is not present in the current UI — no dialog for changing MCL clustering parameters could be found.

## Retrospective

### What worked well
- SAR analysis launches correctly from Bio > Analyze > SAR menu
- "Analyze Peptides" dialog pre-selects AlignedSequence and IC50 columns
- MCL viewer, Most Potent Residues, and Sequence Variability Map all render correctly
- MCL clustering produces expected columns (Cluster, EmbedX/Y, Connectivity, Cluster size)
- Playwright spec passes with waitForFunction for async MCL column creation

### What did not work
- **No wrench button** in the MCL viewer — the scenario mentions "Open settings using wrench button" but no such UI element exists. The MCL viewer contains only standard scatter plot icons (icon-plus, icon-edit for "Add New Column")
- MCL clustering takes significant time (~60-80s in Playwright) — spec needs long timeout

### Suggestions for the platform
- Add a settings/wrench button to the MCL viewer for adjusting clustering parameters (inflation, expansion, etc.)
- Or add MCL parameter controls to the viewer's properties panel

### Suggestions for the scenario
- Clarify which button opens MCL settings — "wrench button" does not exist in current UI
- Consider adding expected cluster count or parameter values to verify
- Note that MCL clustering is async and may take 60+ seconds
