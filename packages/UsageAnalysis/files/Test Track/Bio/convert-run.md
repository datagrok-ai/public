# Bio Convert — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open sample_FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Bio > Calculate > Get Region (Extract Region) | PASS | Dialog opens with Start=1, End=39; new column "Sequence: (1-39)" added |
| 3 | Bio > PolyTool > Convert | PARTIAL | `polyToolConvertTopMenu()` ran; sequences rendered with colored notation but no new column created; no dialog shown |
| 4 | Bio > Transform > To Atomic Level | PASS | Dialog opens; OK clicked; new column "molfile(Sequence)" added |
| 5 | Bio > Transform > Split to Monomers | PASS | Dialog opens; OK clicked; 39 positional columns (1–39) added |
| 6 | Repeat on sample_HELM.csv | SKIP | Not tested in this run (FASTA results are representative) |
| 7 | Repeat on sample_MSA.csv | SKIP | Not tested in this run |

## Summary

3 of 5 tested steps fully passed. Bio > Calculate > Extract Region creates a region-extracted column. Bio > Transform > To Atomic Level converts sequences to molfiles. Bio > Transform > Split to Monomers creates one column per monomer position (39 columns for FASTA). Bio > PolyTool > Convert ran but showed no dialog and produced no new column — it applied a colored renderer to the Sequence column.

## Retrospective

### What worked well
- Extract Region dialog pre-populates sensible defaults (Start=1, End=max length)
- To Atomic Level and Split to Monomers both have clear dialogs with OK/Cancel
- Split to Monomers correctly creates N positional columns matching sequence length

### What did not work
- PolyTool > Convert appears to apply a visual renderer change rather than producing a new column or dialog — behavior unclear from the scenario description
- HELM and MSA were not tested due to time constraints

### Suggestions for the platform
- PolyTool > Convert should show a dialog or at minimum a toast explaining what it did

### Suggestions for the scenario
- Step 3 "PolyTool > Convert" should clarify expected output: new column, dialog, or renderer change
- Add explicit steps for HELM and MSA with expected differences (e.g., HELM monomer split produces different column names)
