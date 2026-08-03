# Chem Info Panels — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open smiles.csv linked dataset | 12s | PASS | PASSED | 1000 rows, canonical_smiles detected as Molecule |
| 2 | Select canonical_smiles column header | 3s | PASS | PASSED | 11 column-level panes: Details, Filter, Colors, Style, Settings, Plots, Advanced, Sticky meta, Chemistry (Rendering, Highlight) + Dev |
| 3 | Expand all column info panels | 7s | PASS | PASSED | No warnings, no error balloons; all panes render content |
| 4 | Switch to chembl-scaffolds.csv → Chemistry → Rendering | 4s | PASS | n/a | Rendering sub-pane exposes Structures, Scaffold column, Highlight scaffold, Regen coords, Filter type inputs |
| 5 | Click first molecule cell → molecule context + expand panes | 10s | PASS | PASSED | 25 molecule-level panes incl. Chemistry, Descriptors, Properties, MPO, Structural Alerts, Toxicity, DrugBank, CDD Vault, ChEMBL, Chemspace, PubChem, SureChEMBL, Structure, 2D/3D Structure; no errors |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m |
| grok-browser execution (scenario steps) | 1m 10s |
| Execute via grok-browser (total) | 3m 10s |
| Spec file generation | 1m |
| Spec script execution | 31s |
| **Total scenario run (with model)** | ~5m |

## Summary

Info panels work correctly on smiles.csv: column-level (Details, Filter, Colors, Style, Chemistry with Rendering/Highlight sub-panes) and molecule-level (Chemistry, Biology, Databases, Structure groups with many child panels) both populate without warnings or error balloons on dev.datagrok.ai. Automated Playwright replay via `auth` cookie storageState passes in ~31s.

## Retrospective

### What worked well
- 1000-row smiles dataset loads cleanly; Molecule semtype detected on canonical_smiles
- Context panel showContextPanel=true + assigning `grok.shell.o` to `DG.Column` or `DG.SemanticValue.fromTableCell` correctly rebuilds the pane tree
- All 25 molecule-level panes expand without producing a single warning or error balloon
- Chemistry > Rendering input host (`input-host-Scaffold-column`, `input-host-Highlight-scaffold`) is reliably reachable by name

### What did not work
- `grok.dapi.files.readCsv(...)` on a `.sdf` path parses the binary as a single "col 1" CSV (91k rows, no Molecule column) — SDF formats cannot be exercised through this code path
- Spec first failed with `ENOENT` when storageState was pointed inside `test-results/` (Playwright wipes that directory pre-run)

### Suggestions for the platform
- `grok.dapi.files.readCsv` should fail loudly (or auto-dispatch to `readSdf`) when the extension is `.sdf` rather than silently producing a degenerate single-column frame
- Consider a public `grok.shell.context.showColumn(col)` / `showCell(cell)` helper so automation doesn't need `grok.shell.o = DG.SemanticValue.fromTableCell(...)` gymnastics

### Suggestions for the scenario
- Call out that the Context Panel must be enabled first (`grok.shell.windows.showContextPanel = true` / Alt+P) — in Tabs mode it is hidden by default
- Drop the "use chembl-scaffolds.cvs" typo (`.cvs` → `.csv`)
- Enumerate the expected pane names for both the column- and molecule-level contexts so QA has an explicit checklist
- Split the SDF/MOL coverage into a separate scenario that uses a real SDF loader
