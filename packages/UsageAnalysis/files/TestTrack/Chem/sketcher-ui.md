---
feature: chem
sub_features_covered: [chem.sketcher, chem.sketcher.cell-editor]
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

Per chain YAML (`scenario-chains/chem.yaml` rev 2 footer note (d)): split out from the
original `sketcher.md` Check-block as `target_layer: manual`, `pyramid_layer: ui-smoke` —
this is a visual-fidelity invariant that requires human visual comparison of stereo
rendering before/after a Highlights structure is applied. The Migrator carries the
invariant verbatim; the Automator does **not** turn this into a playwright spec (any
attempt to verify "stereochemistry kept" via DOM assertions would require a structural
render-comparison library outside the qa-pw harness's standard helper surface — flagged
as a future-work atlas-curator candidate, not in scope here).

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

- **`coverage_type: smoke`** — single visual-fidelity invariant that asserts the
  Highlights overlay does not corrupt stereochemistry depiction in a SMILES dataset.
  Manual-layer ui-smoke; no DOM assertion library can reasonably automate "wedge bond
  remains visible" or "E/Z marking remains visible" without a pixel-diff harness or a
  vector-render comparison library that the section does not currently use. Per chain
  rev 2 footer note (d): `target_layer: manual`, `pyramid_layer: ui-smoke`. The chain-
  level A-STRUCT-02 (`coverage_type: edge` or `perf` at chain level) is satisfied by
  the `bug_focused_candidates[]` + r-group-analysis.md edge scenario elsewhere in the
  chain; no SR carryforward needed on this manual-layer scenario.
- **No playwright automation expected.** This scenario lives at `target_layer: manual`.
  The Migrator surfaces the invariant verbatim; QA runs it as a manual visual check
  per release.
- **#2448 origin.** Originally a "Check:" sub-block at the trailing end of the
  pre-rev-2 `sketcher.md` body; lifted out as a separate ui-only scenario per chain
  rev 2 footer note (d) (Olena 2026-05-11) for `target_layer: manual`,
  `pyramid_layer: ui-smoke` placement.
- **Disjoint from `sketcher.md`.** The main `sketcher.md` (now `coverage_type:
  regression`, `target_layer: playwright`) covers the canonical sketcher cell-editor
  walk (Favorites / Recent / Copy / Paste / backend enumeration). This scenario carries
  only the stereo-preservation visual invariant on the Highlights surface — no
  sketcher modal interaction.
- **Dataset provisioning note.** The `SMILES_highlighted.csv` dataset's exact path
  inside the platform's file shares is unverified at migration time (the original
  body had a TODO: "Add to linked datasets"). Preferred path
  `System:AppData/Chem/tests/SMILES_highlighted.csv` is the section convention for
  Chem test datasets; if absent the QA owner provisions the file as a Setup-time
  prerequisite. Flagged in migration report Unresolved ambiguities.
- **WideSmokeTest sibling.** The section's WideSmokeTest folder
  (`public/packages/UsageAnalysis/files/TestTrack/WideSmokeTest/Chem/`) contains
  ui-only smoke notes for Activity Cliffs / Info Panels / R-Groups Analysis. This
  scenario sits at `public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher-ui.md`
  (alongside its parent `sketcher.md` per the decomposition idiom) rather than under
  WideSmokeTest, mirroring the chain rev 2 output_plan target placement.
