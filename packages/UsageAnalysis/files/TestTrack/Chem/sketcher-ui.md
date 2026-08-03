---
feature: chem
target_layer: manual
coverage_type: smoke
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher.md
migration_date: 2026-05-11
related_bugs: []
---

# Chem | Sketcher — Stereochemistry preserved on SMILES Highlight (ui-only)

UI-only manual smoke for the **#2448 invariant**: substructure highlights applied to a
SMILES-format molecule column must preserve the stereochemistry rendering of the underlying
structures. The depicted molecules in the dataframe grid (with stereo bond markings) must
render identically whether or not a Highlights structure is active — the highlight overlay
must not flatten / strip / re-canonicalize the molecule's stereo display.

This scenario is split out from `sketcher.md`'s original Check-block and kept as a manual-only
visual-fidelity check, because verifying "stereochemistry kept" would require a structural
render-comparison library the test harness doesn't currently have — this stays a human visual
comparison for now.

## Setup

1. **Provision linked dataset.** `SMILES_highlighted.csv` — a SMILES-notation dataset
   curated to expose stereochemistry rendering (chiral centers, double-bond E/Z markings,
   wedge bonds) such that visual stereo features are inspectable in the rendered grid
   cells.
   - **Preferred path** (try first): `System:AppData/Chem/tests/SMILES_highlighted.csv`.
   - **Fallback (provisioning requirement)**: if the file is absent from the platform's
     System file shares, document the provisioning step as a manual-setup blocker — the
     QA owner uploads `SMILES_highlighted.csv` to a Files share before running this
     scenario. The exact selection of the curated SMILES rows is unspecified by the
     original body; a small set (~5-10 rows) covering chiral C atoms, E/Z double bonds,
     and explicit wedge bonds is sufficient.
2. **Confirm Chem package is loaded** so that:
   - The Cheminformatics > Highlights > Sketch entry is available in the Context panel
     for a molecule column.
   - The Substructure Filter on `isosmiles` column / scaffold tree is available as an
     alternate "make these structures highlighted" entry point.
3. **No fixture consumed.** Independent scenario.

## Scenario

### #2448 — Stereochemistry preserved on SMILES highlight (ui-only)

1. Open `SMILES_highlighted.csv` from its provisioned path.
2. Open the **Context panel** (Window > Context panel, or equivalent shortcut).
3. Click the column header of the molecule column. The Context panel populates with the
   per-column action panels for a molecule column.
4. In the Context panel, scroll down to **Cheminformatics > Highlights > Sketch**.
5. In the Sketch input, enter `C1CCCCC1` (cyclohexane) and press **Enter**. The structure
   should be highlighted in the dataframe (cells whose molecules contain the cyclohexane
   substructure get the visual highlight overlay).
6. Make the (already-cyclohexane-highlighted) structures additionally highlighted via one
   of the alternate entry points — EITHER:
   - apply the **scaffold tree** to highlight the same substructure across the dataset;
     OR
   - apply the **structure filter** on the `isosmiles` column to filter / highlight rows
     by substructure match.

**Expected result (visual, human review):** the highlighted structures are displayed in
the dataframe grid in the same way as their non-highlighted counterparts, **with
stereochemistry kept**. Specifically:
- Chiral wedge bonds (solid wedge / dashed wedge) remain visible on the affected atoms
  after the highlight overlay is applied.
- E/Z double-bond markings (cis/trans depiction) remain visible after the highlight
  overlay is applied.
- The molecule's overall depiction matches its non-highlighted form aside from the
  highlight color overlay — no atoms are repositioned, no bonds are reshaped, no
  stereo features are stripped.

If any structure's stereo features change appearance when a Highlights structure is
active (vs. when no Highlights is active), the scenario fails — this is the #2448
invariant.

## Notes

- **Disjoint from `sketcher.md`.** The main `sketcher.md` (playwright-automated) covers the
  canonical sketcher cell-editor walk (Favorites / Recent / Copy / Paste / backend
  enumeration). This scenario carries only the stereo-preservation visual invariant on the
  Highlights surface — no sketcher modal interaction.
