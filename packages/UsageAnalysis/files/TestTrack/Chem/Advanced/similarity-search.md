---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [chem.cp.similarity-search-viewer]
realizes: [chem.search.similarity-search, chem.chem-similarity-search]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/similarity-search.md
migration_date: 2026-05-11
migration_report: similarity-search-migration-report.md
related_bugs: []
---

# Similarity Search

End-to-end exercise of the Chem Similarity Search viewer property-panel
surface on `smiles.csv`. Single phase, four beats:

- **Add** the Chem Similarity Search viewer via **Top menu > Chem >
  Search > Similarity Search...** — search results display with the
  most-similar molecules listed.
- **Open** the viewer's property panel via the gear icon on the Chem
  Similarity Search panel header.
- **Exercise** the six property-panel knobs in turn — **Fingerprint**,
  **Limit**, **Distance Metric**, **Size**, **Molecule Properties**,
  **Cutoff** — verifying each modification re-queries the viewer
  without errors or crashes.

Exercises the viewer end-to-end across all six property-axis dimensions.

## Setup

1. **Provision linked dataset.** The scenario consumes a single bundled
   Datagrok file (no external provisioning required):
   - `System:DemoFiles/chem/smiles.csv` — primary SMILES dataset
     (`canonical_smiles` molecule column, the input to the Chem
     Similarity Search viewer).
2. **Confirm Chem package is loaded** so that **Chem > Search >
   Similarity Search...** is registered as a top-menu entry (package
   source `Chem/src/package.ts#L543`) and the Chem Similarity Search
   viewer is registered.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. The
   dataset is opened fresh inside the scenario.

## Scenarios

### Similarity Search viewer add + property-panel six-knob walk

End-to-end walk: open dataset (via TestTrack 'star' icon "Open test
data") → add viewer via top menu → open gear / property panel → exercise
six property knobs (Fingerprint, Limit, Distance Metric, Size, Molecule
Properties, Cutoff) in turn.

1. Open the `smiles.csv` dataset by pressing the **star** icon in
   TestTrack (tooltip: **Open test data**). The
   `System:DemoFiles/chem/smiles.csv` dataset opens into a table view.
   Verify the `canonical_smiles` column auto-detects as a Molecule
   column (RDKit renderer renders cells; SMILES auto-detection
   populates units / semType / cell renderer).
2. From the top menu, run **Chem > Search > Similarity Search...**.
3. Verify: the **Chem Similarity Search** viewer is added to the
   active table view. The search results display the most-similar
   molecules from `smiles.csv` (top-N hits ordered by similarity
   score; default Fingerprint Morgan + default Distance Metric
   Tanimoto + default Limit). No console errors fire during viewer
   add or initial query.
4. Click the **gear** icon on the Chem Similarity Search panel header
   to open the viewer's property settings. Verify: the property panel
   renders on the Context Panel (viewer-properties mode) without
   errors; the six property controls — **Fingerprint**, **Limit**,
   **Distance Metric**, **Size**, **Molecule Properties**, **Cutoff**
   — are visible and interactable.
5. **Fingerprint knob.** Change the **Fingerprint** property and
   re-test the viewer for each available fingerprint type (Morgan,
   RDKit, MACCS, AtomPair, TopologicalTorsion, Pattern). Verify: each
   fingerprint change re-queries the viewer and renders the result set
   ("Most similar structures" content updates); no console errors
   appear and the viewer does not crash for any fingerprint type.
6. **Limit knob.** Increase and decrease the **Limit** property
   (number of hits returned). Verify: the rendered hit count updates
   to reflect the new limit; the viewer re-queries cleanly; no
   console errors.
7. **Distance Metric knob.** Change the **Distance Metric** property
   and re-test the viewer for each available metric (Tanimoto,
   Asymmetric, Cosine, Sokal). Verify: each metric change re-queries
   the viewer and renders updated similarity scores; no console
   errors; no viewer crashes for any metric.
8. **Size knob.** Test the **Size** property across all available
   values (**small**, **normal**, **large**). Verify: the rendered
   hit-tile size updates accordingly; viewer layout adapts without
   errors; no console errors for any size value.
9. **Molecule Properties knob.** Add a few entries from the
   **Molecule Properties** multi-select list. Verify: the selected
   properties are appended to the "Most similar structures"
   information for each hit tile (per-hit property strip rendered
   beneath the molecule); no console errors.
10. **Cutoff knob.** Set the **Cutoff** property to **1**. Verify:
    only one molecule remains in the Chem Similarity Search viewer
    when the cutoff equals 1 (the exact-match self-hit at similarity
    1.0); the viewer re-queries cleanly; no console errors.
11. **Final assertion.** Across all six property modifications
    (steps 5-10), verify that no errors or crashes occurred — the
    viewer remained interactable throughout the property-panel walk
    and re-queried correctly for each knob change.

## Notes

- **Scope note.** Step 9 intentionally selects a few Molecule Properties entries (e.g. MW, LogP, HBA)
  rather than every available property — the assertion is that selected properties appear on the hit
  tiles, not that every available property has been exercised.
- **Cutoff = 1 assumption.** Step 10's "only one molecule should remain" assertion relies on
  `smiles.csv` containing no duplicate SMILES strings — the query molecule matches itself at similarity
  1.0, and all other rows are assumed below 1.0. If the dataset ever gains exact-duplicate rows, the hit
  count at cutoff = 1 may exceed 1.
