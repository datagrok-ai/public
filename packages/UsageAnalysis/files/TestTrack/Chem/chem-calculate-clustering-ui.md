---
feature: chem
realizes_atlas: []   # untyped: no matching atlas type
priority: p2
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  Butina Cluster (Chem | Analyze | Butina Cluster...) and the BitBIRCH-vs-Butina
  contrast are executed by hand. Butina is a server-side Python script; live MCP
  recon on dev 2026-08-19 opened the dialog and dispatched OK cleanly, but the job
  never joined its cluster (Butina) column within ~90s even on a 50-row extract —
  the dialog stayed open with no error, balloon, or column. The server script
  environment does not complete the job in a bounded, headless-automatable window
  on dev, so asserting the column would be a false-RED (the trigger fires but the
  observable never lands). BitBIRCH, Cluster MCS, and Similarity Matrix (Scenarios
  1, 3, 4) are automated in the paired spec.
---

# Chem — Butina Cluster (manual) and BitBIRCH↔Butina contrast

Companion to `chem-calculate-clustering.md` / `chem-calculate-clustering-spec.ts`.

## Setup

1. Open Datagrok, sign in, and open `System:DemoFiles/chem/smiles.csv`. Wait until the
   Molecule column renders structures. Record the baseline row count from the status bar.
2. Run BitBIRCH Clustering first (Chem | Calculate | BitBIRCH Clustering..., accept defaults,
   OK) so a `Cluster (BitBIRCH)` column is present for the contrast in step below.

## Scenario 2 (manual): Butina Cluster — independent clustering engine

1. From the top menu, select Chem | Analyze | Butina Cluster.... The Butina Molecules
   Clustering dialog opens (a Python-script function-call dialog).
2. Confirm the molecule column is auto-detected. Accept the default distance cutoff. Click OK.
3. Wait for the server-side script to finish. A `cluster (Butina)` column is joined to the table.
   Verify:
   - The number of distinct values in `cluster (Butina)` is greater than 1 and strictly less
     than the row count.
   - The baseline row count is unchanged.
   - No error balloon appears and no console errors fire.

## BitBIRCH ↔ Butina contrast (manual)

1. With both `Cluster (BitBIRCH)` and `cluster (Butina)` columns present on the same table,
   compare them row by row.
2. At least one row must have a different cluster assignment between the two columns —
   an exact match on every row would indicate one engine re-used the other's cached result.
   This confirms the two engines ran independently.

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/chem/smiles.csv"]
}
